"""
soy_globin_core — Modal-free implementation of the soybean globin-family pipeline.

Every stage here is a plain function with no Modal import, so the identical code
runs (a) inside the Modal app in ``soy_globin_modal.py`` and (b) locally for
validation. Stages:

    1. fetch_inputs        Wm82.gnm4.ann1.T8TQ primary proteome + gene-model GFF3
                           (SoyBase) and the Pfam PF00042 (Globin) HMM (InterPro).
    2. select_family       hmmsearch --cut_ga, union with the four focal Lb genes.
    3. align_and_tree      MAFFT L-INS-i  ->  IQ-TREE 2 (ModelFinder + UFBoot + SH-aLRT).
    4. pairwise_identity   %ID from the MAFFT alignment.
    5. esm2_embeddings     ESM2-650M mean-pooled residue embeddings -> cosine distance.
    6. classify_pairs      tandem / proximal / dispersed from GFF3 gene order.

Reference build: Wm82.gnm4.ann1.T8TQ  ==  Phytozome Wm82.a4.v1.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

SOYBASE_BASE = (
    "https://data.soybase.org/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"
)
PROTEOME_URL = f"{SOYBASE_BASE}/glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz"
GFF3_URL = f"{SOYBASE_BASE}/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_exons.gff3.gz"
PFAM_HMM_URL = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00042?annotation=hmm"

#: Focal leghemoglobins (Wm82.a4.v1 gene IDs) -> symbol.
FOCAL_GENES = {
    "Glyma.10G199100": "Lba",
    "Glyma.10G199000": "Lbc1",
    "Glyma.20G191200": "Lbc2",
    "Glyma.10G198800": "Lbc3",
}

ESM2_MODEL = "facebook/esm2_t33_650M_UR50D"

#: A same-chromosome pair with at most this many intervening protein-coding genes
#: is called a tandem duplicate (the convention used by MCScanX / PGDD).
TANDEM_MAX_INTERVENING = 10
#: Same chromosome, more than TANDEM_MAX_INTERVENING apart but within this bp
#: window -> "proximal"; beyond it -> "dispersed".
PROXIMAL_MAX_BP = 1_000_000

# Strip the assembly/annotation prefix that SoyBase puts on every identifier.
_ID_PREFIX = re.compile(r"^glyma\.Wm82\.gnm4\.ann1\.")
_CHR_RE = re.compile(r"Gm(\d+)$")


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #


def sha256(path: str | Path, buf: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(buf), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd: list[str], log: str | Path | None = None, **kw) -> subprocess.CompletedProcess:
    """Run a subprocess, raising with captured stderr on failure."""
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if log:
        Path(log).write_text(
            f"$ {' '.join(cmd)}\n[exit {proc.returncode}  {time.time()-t0:.1f}s]\n"
            f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}\n"
        )
    if proc.returncode != 0:
        raise RuntimeError(
            f"command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr[-4000:]}"
        )
    return proc


_VERSION_RE = re.compile(r"\d+\.\d+")


def tool_version(cmd: list[str]) -> str:
    """First line of a tool's banner that actually carries a version number.

    hmmsearch's ``-h`` opens with a decoration line and puts the version on
    the next one, so taking line 0 unconditionally records nothing useful.
    """
    try:
        p = subprocess.run(cmd, capture_output=True, text=True)
        lines = [l.strip() for l in (p.stdout or p.stderr).strip().splitlines()]
        for line in lines[:6]:
            if _VERSION_RE.search(line):
                return line.lstrip("# ").strip()
        return lines[0] if lines else "<no output>"
    except Exception as exc:  # pragma: no cover - diagnostic only
        return f"<unavailable: {exc}>"


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read (optionally gzipped) FASTA into {header_id: sequence}, '*' stripped."""
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs: dict[str, str] = {}
    name, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).replace("*", "").upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).replace("*", "").upper()
    return seqs


def write_fasta(seqs: dict[str, str], path: str | Path, width: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def gene_of(protein_id: str) -> str:
    """'glyma.Wm82.gnm4.ann1.Glyma.10G199100.1' -> 'Glyma.10G199100'."""
    core = _ID_PREFIX.sub("", protein_id)
    return core.rsplit(".", 1)[0] if re.search(r"\.\d+$", core) else core


def label_of(gene_id: str) -> str:
    """Tree/matrix label: gene ID, suffixed with the symbol for focal genes."""
    sym = FOCAL_GENES.get(gene_id)
    return f"{gene_id}_{sym}" if sym else gene_id


# --------------------------------------------------------------------------- #
# Stage 1 — inputs
# --------------------------------------------------------------------------- #


def fetch_inputs(datadir: str | Path, force: bool = False) -> dict[str, str]:
    """Download proteome, GFF3 and the PF00042 HMM. Cached by presence on disk."""
    datadir = Path(datadir)
    datadir.mkdir(parents=True, exist_ok=True)
    targets = {
        "proteome": (PROTEOME_URL, datadir / "proteins_primary.faa.gz"),
        "gff3": (GFF3_URL, datadir / "gene_models_exons.gff3.gz"),
        "pfam_hmm_gz": (PFAM_HMM_URL, datadir / "PF00042.hmm.gz"),
    }
    out: dict[str, str] = {}
    for key, (url, dest) in targets.items():
        if force or not dest.exists() or dest.stat().st_size == 0:
            print(f"[fetch] {url}", flush=True)
            with urllib.request.urlopen(url, timeout=300) as r, open(dest, "wb") as fh:
                shutil.copyfileobj(r, fh)
        out[key] = str(dest)

    hmm = datadir / "PF00042.hmm"
    if force or not hmm.exists():
        with gzip.open(out["pfam_hmm_gz"], "rb") as src, open(hmm, "wb") as dst:
            shutil.copyfileobj(src, dst)
    out["pfam_hmm"] = str(hmm)
    return out


def hmm_metadata(hmm_path: str | Path) -> dict[str, str]:
    meta: dict[str, str] = {}
    with open(hmm_path) as fh:
        for line in fh:
            if line.startswith("//"):
                break
            for key in ("NAME", "ACC", "DESC", "LENG", "GA"):
                if line.startswith(key):
                    meta[key.lower()] = line[len(key) :].strip()
    return meta


# --------------------------------------------------------------------------- #
# Stage 2 — family selection
# --------------------------------------------------------------------------- #


def parse_tblout(path: str | Path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            rows.append(
                {
                    "protein_id": f[0],
                    "query": f[2],
                    "evalue_full": float(f[4]),
                    "score_full": float(f[5]),
                    "evalue_best_dom": float(f[7]),
                    "score_best_dom": float(f[8]),
                }
            )
    return pd.DataFrame(rows)


def parse_domtblout(path: str | Path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            rows.append(
                {
                    "protein_id": f[0],
                    "prot_len": int(f[2]),
                    "hmm_len": int(f[5]),
                    "i_evalue": float(f[12]),
                    "dom_score": float(f[13]),
                    "hmm_from": int(f[15]),
                    "hmm_to": int(f[16]),
                    "ali_from": int(f[17]),
                    "ali_to": int(f[18]),
                }
            )
    df = pd.DataFrame(rows)
    if not df.empty:
        df["hmm_coverage"] = (df.hmm_to - df.hmm_from + 1) / df.hmm_len
    return df


def select_family(
    proteome: str | Path,
    hmm: str | Path,
    workdir: str | Path,
    cpu: int = 4,
    focal: dict[str, str] | None = None,
) -> tuple[dict[str, str], pd.DataFrame]:
    """hmmsearch --cut_ga for PF00042, then union with the focal Lb genes.

    Returns (labelled sequences, membership table).
    """
    focal = FOCAL_GENES if focal is None else focal
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    tbl, domtbl = workdir / "PF00042.tblout", workdir / "PF00042.domtblout"

    run(
        [
            "hmmsearch", "--cut_ga", "--cpu", str(cpu),
            "--tblout", str(tbl), "--domtblout", str(domtbl),
            "-o", str(workdir / "hmmsearch.out"),
            str(hmm), str(proteome),
        ],
        log=workdir / "hmmsearch.log",
    )

    hits = parse_tblout(tbl)
    dom = parse_domtblout(domtbl)
    best_dom = (
        dom.sort_values("dom_score", ascending=False)
        .groupby("protein_id", as_index=False)
        .first()
        if not dom.empty
        else dom
    )

    all_seqs = read_fasta(proteome)
    prot_by_gene: dict[str, str] = {}
    for pid in all_seqs:
        prot_by_gene.setdefault(gene_of(pid), pid)

    hit_pids = set(hits.protein_id)
    rescued = []
    for gene in focal:
        pid = prot_by_gene.get(gene)
        if pid is None:
            raise KeyError(f"focal gene {gene} not found in proteome {proteome}")
        if pid not in hit_pids:
            hit_pids.add(pid)
            rescued.append(gene)
    if rescued:
        print(f"[select] focal genes added below the PF00042 GA cutoff: {rescued}")

    rows = []
    for pid in sorted(hit_pids):
        gene = gene_of(pid)
        h = hits[hits.protein_id == pid]
        d = best_dom[best_dom.protein_id == pid] if not best_dom.empty else best_dom
        rows.append(
            {
                "gene_id": gene,
                "label": label_of(gene),
                "protein_id": pid,
                "symbol": focal.get(gene, ""),
                "is_focal": gene in focal,
                "prot_len": len(all_seqs[pid]),
                "pfam_hit": pid not in {prot_by_gene[g] for g in rescued},
                "evalue_full": float(h.evalue_full.iloc[0]) if len(h) else np.nan,
                "score_full": float(h.score_full.iloc[0]) if len(h) else np.nan,
                "dom_score": float(d.dom_score.iloc[0]) if len(d) else np.nan,
                "hmm_coverage": float(d.hmm_coverage.iloc[0]) if len(d) else np.nan,
                "ali_from": int(d.ali_from.iloc[0]) if len(d) else -1,
                "ali_to": int(d.ali_to.iloc[0]) if len(d) else -1,
            }
        )
    members = pd.DataFrame(rows).sort_values(
        ["is_focal", "score_full"], ascending=[False, False]
    ).reset_index(drop=True)

    seqs = {r.label: all_seqs[r.protein_id] for r in members.itertuples()}
    return seqs, members


# --------------------------------------------------------------------------- #
# Stage 3 — alignment and phylogeny
# --------------------------------------------------------------------------- #


def run_mafft(fasta: str | Path, out_aln: str | Path, threads: int = 4) -> str:
    """MAFFT L-INS-i (accurate, suitable for <~200 sequences)."""
    cmd = ["mafft", "--localpair", "--maxiterate", "1000",
           "--anysymbol", "--thread", str(threads), str(fasta)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"mafft failed:\n{proc.stderr[-4000:]}")
    Path(out_aln).write_text(proc.stdout)
    Path(str(out_aln) + ".log").write_text(proc.stderr)
    return str(out_aln)


def iqtree_exe() -> str:
    """Bioconda ships the v2 binary as ``iqtree2`` or ``iqtree`` depending on build."""
    exe = shutil.which("iqtree2") or shutil.which("iqtree")
    if exe is None:
        raise RuntimeError("neither iqtree2 nor iqtree found on PATH")
    return exe


def run_iqtree(
    aln: str | Path, prefix: str | Path, threads: str | int = "AUTO", seed: int = 20240601
) -> dict[str, str]:
    """IQ-TREE 2: ModelFinder + 1000 UFBoot replicates + 1000 SH-aLRT replicates."""
    exe = iqtree_exe()
    run(
        [exe, "-s", str(aln), "-m", "MFP", "-B", "1000", "--alrt", "1000",
         "-T", str(threads), "--seed", str(seed), "--prefix", str(prefix), "-redo"],
        log=str(prefix) + ".cmd.log",
    )
    return {
        "treefile": f"{prefix}.treefile",
        "iqtree_report": f"{prefix}.iqtree",
        "contree": f"{prefix}.contree",
    }


def best_model(iqtree_report: str | Path) -> str:
    for line in Path(iqtree_report).read_text().splitlines():
        if line.startswith("Best-fit model according to"):
            return line.split(":", 1)[1].strip()
    return "unknown"


# --------------------------------------------------------------------------- #
# Stage 4 — pairwise identity
# --------------------------------------------------------------------------- #


def pairwise_identity(aln_path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Percent identity from the MSA.

    ``pid_aligned``  : identities / columns where BOTH sequences have a residue
                       (the usual "pairwise identity" of an alignment).
    ``pid_shorter``  : identities / length of the shorter ungapped sequence
                       (penalises partial/truncated globins).
    """
    aln = read_fasta(aln_path)
    labels = list(aln)
    arr = np.array([list(aln[l]) for l in labels])
    is_res = (arr != "-") & (arr != ".")
    ungapped_len = is_res.sum(axis=1)

    n = len(labels)
    pid_a = np.full((n, n), 100.0)
    pid_s = np.full((n, n), 100.0)
    long_rows = []
    for i, j in combinations(range(n), 2):
        both = is_res[i] & is_res[j]
        n_both = int(both.sum())
        n_id = int((arr[i][both] == arr[j][both]).sum())
        a = 100.0 * n_id / n_both if n_both else 0.0
        s = 100.0 * n_id / int(min(ungapped_len[i], ungapped_len[j]))
        pid_a[i, j] = pid_a[j, i] = a
        pid_s[i, j] = pid_s[j, i] = s
        long_rows.append(
            {
                "label_a": labels[i], "label_b": labels[j],
                "n_aligned_cols": n_both, "n_identical": n_id,
                "pid_aligned": a, "pid_shorter": s,
                "len_a": int(ungapped_len[i]), "len_b": int(ungapped_len[j]),
            }
        )
    mat = pd.DataFrame(pid_a, index=labels, columns=labels)
    return mat, pd.DataFrame(long_rows)


# --------------------------------------------------------------------------- #
# Stage 5 — ESM2 embeddings
# --------------------------------------------------------------------------- #


def esm2_embeddings(
    seqs: dict[str, str],
    model_name: str = ESM2_MODEL,
    device: str | None = None,
    batch_size: int = 8,
) -> tuple[list[str], np.ndarray]:
    """Mean-pooled final-layer ESM2 embeddings, averaged over residues only
    (BOS/EOS/PAD excluded)."""
    import torch
    from transformers import AutoTokenizer, EsmModel

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    dtype = torch.float16 if device == "cuda" else torch.float32

    tok = AutoTokenizer.from_pretrained(model_name)
    # add_pooling_layer=False: we pool ourselves from last_hidden_state, and the
    # checkpoint has no pooler weights (loading it emits a spurious
    # "newly initialized" warning for parameters this code never reads).
    model = (
        EsmModel.from_pretrained(model_name, torch_dtype=dtype, add_pooling_layer=False)
        .to(device)
        .eval()
    )

    labels = list(seqs)
    out = np.zeros((len(labels), model.config.hidden_size), dtype=np.float32)
    with torch.no_grad():
        for start in range(0, len(labels), batch_size):
            chunk = labels[start : start + batch_size]
            enc = tok([seqs[l] for l in chunk], return_tensors="pt", padding=True)
            enc = {k: v.to(device) for k, v in enc.items()}
            hidden = model(**enc).last_hidden_state  # (B, L, H)

            # residue mask = attention mask minus the BOS and EOS positions
            mask = enc["attention_mask"].clone()
            mask[:, 0] = 0
            eos = enc["attention_mask"].sum(dim=1) - 1
            mask[torch.arange(mask.size(0), device=device), eos] = 0
            m = mask.unsqueeze(-1).to(hidden.dtype)

            pooled = (hidden * m).sum(dim=1) / m.sum(dim=1)
            out[start : start + len(chunk)] = pooled.float().cpu().numpy()
    return labels, out


def cosine_distance_matrix(labels: list[str], emb: np.ndarray) -> pd.DataFrame:
    norm = emb / np.linalg.norm(emb, axis=1, keepdims=True)
    dist = 1.0 - norm @ norm.T
    np.fill_diagonal(dist, 0.0)
    return pd.DataFrame(dist, index=labels, columns=labels)


# --------------------------------------------------------------------------- #
# Stage 6 — genomic context and paralog-pair classification
# --------------------------------------------------------------------------- #


def parse_gff_genes(gff: str | Path) -> pd.DataFrame:
    """Gene features from the Wm82.gnm4.ann1 GFF3, with per-seqid rank order."""
    opener = gzip.open if str(gff).endswith(".gz") else open
    rows = []
    with opener(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(
                kv.split("=", 1) for kv in f[8].split(";") if "=" in kv
            )
            gene = attrs.get("Name") or _ID_PREFIX.sub("", attrs.get("ID", ""))
            m = _CHR_RE.search(f[0])
            rows.append(
                {
                    "gene_id": gene,
                    "seqid": f[0],
                    "chrom": f"Gm{int(m.group(1)):02d}" if m else f[0].split(".")[-1],
                    "is_scaffold": m is None,
                    "start": int(f[3]),
                    "end": int(f[4]),
                    "strand": f[6],
                }
            )
    genes = pd.DataFrame(rows)
    genes = genes.sort_values(["seqid", "start", "end"]).reset_index(drop=True)
    genes["rank_on_seqid"] = genes.groupby("seqid").cumcount()
    return genes


def classify_pairs(
    members: pd.DataFrame,
    genes: pd.DataFrame,
    identity_long: pd.DataFrame | None = None,
    cosine: pd.DataFrame | None = None,
    max_intervening: int = TANDEM_MAX_INTERVENING,
    proximal_max_bp: int = PROXIMAL_MAX_BP,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Annotate every within-family pair as tandem / proximal / dispersed.

    tandem    same seqid and <= ``max_intervening`` protein-coding genes between
              them (the chr10 Lb cluster falls here)
    proximal  same seqid, more intervening genes but <= ``proximal_max_bp`` apart
    dispersed different seqid (e.g. the chr20 Lbc2 copy vs the chr10 cluster),
              or same seqid but beyond the proximal window
    """
    ctx = members.merge(genes, on="gene_id", how="left", validate="one_to_one")
    missing = ctx[ctx.seqid.isna()].gene_id.tolist()
    if missing:
        raise KeyError(f"genes absent from GFF3: {missing}")
    ctx["rank_on_seqid"] = ctx["rank_on_seqid"].astype(int)

    by_label = ctx.set_index("label")
    rows = []
    for a, b in combinations(list(by_label.index), 2):
        ga, gb = by_label.loc[a], by_label.loc[b]
        same = ga.seqid == gb.seqid
        if same:
            n_between = int(abs(ga.rank_on_seqid - gb.rank_on_seqid)) - 1
            gap_bp = int(max(ga.start, gb.start) - min(ga.end, gb.end))
            gap_bp = max(gap_bp, 0)
            if n_between <= max_intervening:
                mode = "tandem"
            elif gap_bp <= proximal_max_bp:
                mode = "proximal"
            else:
                mode = "dispersed"
        else:
            n_between, gap_bp, mode = -1, -1, "dispersed"
        rows.append(
            {
                "label_a": a, "label_b": b,
                "gene_a": ga.gene_id, "gene_b": gb.gene_id,
                "symbol_a": ga.symbol, "symbol_b": gb.symbol,
                "chrom_a": ga.chrom, "chrom_b": gb.chrom,
                "same_seqid": bool(same),
                "n_intervening_genes": n_between,
                "intergenic_bp": gap_bp,
                "duplication_mode": mode,
                "both_focal": bool(ga.is_focal and gb.is_focal),
            }
        )
    pairs = pd.DataFrame(rows)

    if identity_long is not None:
        pairs = pairs.merge(
            identity_long[["label_a", "label_b", "pid_aligned", "pid_shorter",
                           "n_aligned_cols"]],
            on=["label_a", "label_b"], how="left",
        )
    if cosine is not None:
        pairs["esm2_cosine_distance"] = [
            float(cosine.loc[r.label_a, r.label_b]) for r in pairs.itertuples()
        ]
    sort_cols = [c for c in ("both_focal", "pid_aligned") if c in pairs.columns]
    if sort_cols:
        pairs = pairs.sort_values(sort_cols, ascending=False).reset_index(drop=True)
    return pairs, ctx


# --------------------------------------------------------------------------- #
# Tree relabelling
# --------------------------------------------------------------------------- #


def annotate_newick(treefile: str | Path, out: str | Path, ctx: pd.DataFrame) -> str:
    """Append chromosome to every tip label: Glyma.10G199100_Lba -> ..._Gm10."""
    nwk = Path(treefile).read_text()
    for r in ctx.itertuples():
        nwk = re.sub(rf"(?<![\w.]){re.escape(r.label)}(?![\w.])",
                     f"{r.label}_{r.chrom}", nwk)
    Path(out).write_text(nwk)
    return str(out)
