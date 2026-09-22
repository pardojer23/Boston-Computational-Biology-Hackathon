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
from collections.abc import Iterable
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# No constants.
#
# Every value that used to live here — the SoyBase URLs, FOCAL_GENES,
# GENE_SYMBOLS, the ESM2 model id, the duplication thresholds, the identifier
# regexes — is now in ``config/config.yaml`` and reaches these functions as an
# argument. The reason is not tidiness: the previous arrangement meant that
# editing FOCAL_GENES and re-running one module produced a clean-looking table
# scored against a different family definition, with nothing to detect it
# (docs/PIPELINE_AUDIT.md finding 1). As config the same edit changes the
# config digest, which is a declared input of every rule.
#
# These functions take a ``cfg`` (a ``pipeline.config.Config``) but deliberately
# do not import it: ``config`` needs PyYAML, and keeping this module free of
# that dependency is what lets the identical code run inside the Modal image.
# The dependency is duck-typed and checked at the call site.
# --------------------------------------------------------------------------- #

if TYPE_CHECKING:  # pragma: no cover - typing only
    from config import Config


def _compiled(pattern: str | re.Pattern) -> re.Pattern:
    return pattern if isinstance(pattern, re.Pattern) else re.compile(pattern)


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


def gene_of(protein_id: str, id_prefix: str | re.Pattern) -> str:
    """``glyma.Wm82.gnm4.ann1.Glyma.10G199100.1`` -> ``Glyma.10G199100``.

    ``id_prefix`` comes from ``family.identifiers.id_prefix_strip``; it was a
    module constant compiled against this one assembly.
    """
    core = _compiled(id_prefix).sub("", protein_id)
    return core.rsplit(".", 1)[0] if re.search(r"\.\d+$", core) else core


def label_of(gene_id: str, symbols: dict[str, str]) -> str:
    """Tree/matrix label: gene ID, suffixed with its symbol if it has one.

    ``symbols`` is ``cfg.gene_symbols`` — focal and outgroup merged. The map was
    a module constant (``GENE_SYMBOLS``) declared in two places at one point,
    and the two copies disagreeing about whether ``Glyma.10G198900`` carries a
    symbol silently cost 6 of 21 rows on any cross-module join (README §6).
    Passing it in means there is one map, held by the config.
    """
    sym = symbols.get(gene_id)
    return f"{gene_id}_{sym}" if sym else gene_id


def gene_from_label(label: str) -> str:
    """Inverse of ``label_of``: ``Glyma.10G199100_Lba`` -> ``Glyma.10G199100``.

    Needs no family knowledge, so it takes no config.
    """
    return label.split("_", 1)[0]


def canonical_pair_order(labels: Iterable[str]) -> list[tuple[str, str]]:
    """Every unordered pair of ``labels``, in one canonical orientation and order.

    Each pair is oriented so ``gene_from_label(a) < gene_from_label(b)``, and the
    pairs themselves are sorted by that same key.

    Every module that emits a pair table iterates this. That is what makes
    ``label_a``/``label_b`` a usable join key across modules: two modules that
    each call ``itertools.combinations`` over their own member ordering will emit
    the same pair in opposite orientations, and a merge on those columns then
    drops the reversed rows without raising. Four of 21 pairs were reversed this
    way before this function existed.

    ``contracts.canonical_pairs`` implements the same ordering and is what the
    table contracts check against; this stays as the name the science modules
    call, and the test suite asserts the two agree.
    """
    pairs = [tuple(sorted((a, b), key=gene_from_label)) for a, b in combinations(labels, 2)]
    return sorted(pairs, key=lambda p: (gene_from_label(p[0]), gene_from_label(p[1])))


# --------------------------------------------------------------------------- #
# Stage 1 — inputs
# --------------------------------------------------------------------------- #


class ChecksumMismatch(RuntimeError):
    """A downloaded or cached reference file does not match its pinned digest."""


def verify_sha256(path: str | Path, expected: str | None) -> str:
    """Return the file's digest, raising if it does not match ``expected``.

    ``expected=None`` records without asserting — used the first time a file is
    pinned, so the digest can be read off the run manifest and pasted into the
    config.
    """
    got = sha256(path)
    if expected and got != expected:
        raise ChecksumMismatch(
            f"{path}\n  expected sha256 {expected}\n  got      sha256 {got}\n"
            f"The pinned reference has changed upstream. SoyBase and InterPro "
            f"both re-release under the same URL, so this is the only thing "
            f"standing between a silent annotation change and a re-scored "
            f"family. Confirm the new release is what you want, then update "
            f"the digest in config/config.yaml."
        )
    return got


def fetch_file(
    url: str,
    dest: str | Path,
    sha256_expected: str | None = None,
    force: bool = False,
    gunzip_to: str | Path | None = None,
) -> dict:
    """Download one file and verify it. Returns a provenance record.

    One file per call rather than the previous all-three-at-once helper: each
    reference is its own workflow rule with its own output, so a changed HMM
    re-runs the family search without re-downloading the proteome.

    Caching is by presence *and digest*: a file already on disk is verified
    rather than trusted, which is what makes a truncated download — previously
    detectable only when a parser failed much later — fail here instead.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    need = force or not dest.exists() or dest.stat().st_size == 0
    if not need and sha256_expected:
        try:
            verify_sha256(dest, sha256_expected)
        except ChecksumMismatch:
            print(f"[fetch] cached {dest.name} fails its digest; re-downloading",
                  flush=True)
            need = True
    if need:
        print(f"[fetch] {url}", flush=True)
        with urllib.request.urlopen(url, timeout=300) as r, open(dest, "wb") as fh:
            shutil.copyfileobj(r, fh)

    got = verify_sha256(dest, sha256_expected)
    rec = {"url": url, "path": str(dest), "sha256": got,
           "bytes": dest.stat().st_size, "downloaded": bool(need)}

    if gunzip_to is not None:
        out = Path(gunzip_to)
        out.parent.mkdir(parents=True, exist_ok=True)
        if force or need or not out.exists():
            with gzip.open(dest, "rb") as src, open(out, "wb") as dst:
                shutil.copyfileobj(src, dst)
        rec["gunzipped_path"] = str(out)
        rec["gunzipped_sha256"] = sha256(out)
    return rec


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


def hmmsearch(
    hmm: str | Path,
    proteome: str | Path,
    workdir: str | Path,
    prefix: str,
    threshold: str = "cut_ga",
    evalue: float | None = None,
    cpu: int = 4,
) -> dict[str, Path]:
    """One hmmsearch invocation. Returns the paths it wrote.

    Split out of ``select_family`` so the relaxed cutoff check can reuse it
    instead of being a hand-run command whose result lived only in prose
    (README §1.2, audit 1.2).
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "tblout": workdir / f"{prefix}.tblout",
        "domtblout": workdir / f"{prefix}.domtblout",
        "out": workdir / f"{prefix}.hmmsearch.out",
    }
    if threshold == "cut_ga":
        thr = ["--cut_ga"]
    elif threshold == "evalue":
        if evalue is None:
            raise ValueError("threshold='evalue' requires an evalue")
        thr = ["-E", str(evalue)]
    else:
        raise ValueError(f"unknown threshold {threshold!r}; use 'cut_ga' or 'evalue'")

    run(
        ["hmmsearch", *thr, "--cpu", str(cpu),
         "--tblout", str(paths["tblout"]), "--domtblout", str(paths["domtblout"]),
         "-o", str(paths["out"]), str(hmm), str(proteome)],
        log=workdir / f"{prefix}.hmmsearch.log",
    )
    return paths


def cutoff_gap_check(
    hmm: str | Path,
    proteome: str | Path,
    workdir: str | Path,
    n_accepted: int,
    relaxed_evalue: float = 10.0,
    cpu: int = 4,
) -> dict:
    """Evidence that the gathering threshold is not truncating the family.

    Runs the search again at a relaxed E-value and measures the bit-score gap
    between the weakest accepted hit and the best rejected one. A large gap
    means the cutoff falls in empty space rather than through a continuum of
    weaker family members.

    For the committed run: the 8th-best protein in the genome scores 13.5
    (E = 0.23) against 63.2 for the weakest real hit — a ~50-bit gap. That
    number was obtained by hand and written into the README; here it is
    computed, recorded, and subject to a configured pass condition.
    """
    paths = hmmsearch(hmm, proteome, workdir, prefix="relaxed",
                      threshold="evalue", evalue=relaxed_evalue, cpu=cpu)
    hits = parse_tblout(paths["tblout"]).sort_values("score_full", ascending=False)
    scores = hits.score_full.tolist()
    accepted = scores[:n_accepted]
    rejected = scores[n_accepted:]
    weakest_accepted = float(accepted[-1]) if accepted else float("nan")
    best_rejected = float(rejected[0]) if rejected else float("-inf")
    gap = weakest_accepted - best_rejected if rejected else float("inf")
    return {
        "relaxed_evalue": relaxed_evalue,
        "n_hits_relaxed": int(len(hits)),
        "n_accepted": int(n_accepted),
        "weakest_accepted_bitscore": weakest_accepted,
        "best_rejected_bitscore": best_rejected if rejected else None,
        "best_rejected_protein": (hits.protein_id.iloc[n_accepted]
                                  if rejected else None),
        "bitscore_gap": gap,
    }


def select_family(
    proteome: str | Path,
    hmm: str | Path,
    workdir: str | Path,
    cfg: "Config",
    seed: pd.DataFrame | None = None,
    cpu: int = 4,
) -> tuple[dict[str, str], pd.DataFrame]:
    """HMM search over the proteome, unioned with the resolved family seed.

    Returns (labelled sequences, membership table). Family membership is
    *discovered*, not declared — the seed genes are only guaranteed to be
    present, and everything else in the table is whatever the HMM found.

    ``seed`` is the resolver's output (``pipeline/family.py``). It carries the
    genes that must be in the family and which of them were declared focal;
    when no gene is declared focal, every discovered member resolves to focal
    and the family has no outgroup. Passing ``seed=None`` falls back to the
    config's own declaration, which is what the ``gene_ids`` mode reduces to.
    """
    import family as fammod

    if seed is None:
        seed = fammod.resolve(cfg)
    seed_focal = set(seed.loc[seed.is_focal_declared.astype(bool), "gene_id"]) \
        if len(seed) else set()
    # Symbols come from the seed first, then the config's outgroup labels.
    symbols = {**cfg.gene_symbols,
               **{r.gene_id: r.symbol for r in seed.itertuples() if r.symbol}}
    # Only rows flagged must_include are force-included. Outgroup symbols are
    # carried in the seed as LABELS and are not flagged, so naming a gene can
    # never by itself put it in the family — for `pfam` the search alone
    # decides membership.
    must_include = (set(seed.loc[seed.must_include.astype(bool), "gene_id"])
                    if len(seed) else set())
    id_prefix = _compiled(cfg["family.identifiers.id_prefix_strip"])
    acc = cfg["sequence.hmm.accession"]

    paths = hmmsearch(
        hmm, proteome, workdir, prefix=acc,
        threshold=cfg["sequence.select.threshold"],
        evalue=cfg.get("sequence.select.evalue"),
        cpu=cpu,
    )
    tbl, domtbl = paths["tblout"], paths["domtblout"]

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
        prot_by_gene.setdefault(gene_of(pid, id_prefix), pid)

    hit_pids = set(hits.protein_id)
    rescued = []
    if cfg.get("sequence.select.rescue_focal_below_cutoff", True):
        for gene in sorted(must_include):
            pid = prot_by_gene.get(gene)
            if pid is None:
                raise KeyError(f"seed gene {gene} not found in proteome {proteome}")
            if pid not in hit_pids:
                hit_pids.add(pid)
                rescued.append(gene)
        if rescued:
            print(f"[select] seed genes added below the {acc} cutoff: {rescued}")
    else:
        absent = [g for g in sorted(must_include)
                  if prot_by_gene.get(g) not in hit_pids]
        if absent:
            raise KeyError(
                f"seed genes below the {acc} cutoff and rescue is disabled: {absent}"
            )

    rows = []
    for pid in sorted(hit_pids):
        gene = gene_of(pid, id_prefix)
        h = hits[hits.protein_id == pid]
        d = best_dom[best_dom.protein_id == pid] if not best_dom.empty else best_dom
        rows.append(
            {
                "gene_id": gene,
                "label": label_of(gene, symbols),
                "protein_id": pid,
                "symbol": symbols.get(gene, ""),
                # Placeholder; resolved against the seed below, because
                # "focal" depends on whether anything was declared focal at
                # all, which is a property of the seed and not of this row.
                "is_focal": gene in seed_focal,
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
    members = pd.DataFrame(rows)
    # No declared focal subset -> the whole discovered set is the family. This
    # is the `focal`-optional path: legal, and materially different from the
    # declared case, since it leaves no outgroup for the separation check.
    members["is_focal"] = fammod.resolve_focal(seed, members)
    if not seed_focal:
        print(f"[select] no focal subset declared; all {len(members)} discovered "
              f"members are focal, so this family has no outgroup")
    members = members.sort_values(
        ["is_focal", "score_full"], ascending=[False, False]
    ).reset_index(drop=True)

    seqs = {r.label: all_seqs[r.protein_id] for r in members.itertuples()}
    return seqs, members


# --------------------------------------------------------------------------- #
# Stage 3 — alignment and phylogeny
# --------------------------------------------------------------------------- #


def run_mafft(
    fasta: str | Path,
    out_aln: str | Path,
    threads: int = 4,
    args: Iterable[str] = ("--localpair", "--maxiterate", "1000"),
) -> str:
    """MAFFT with the configured mode (default L-INS-i).

    ``args`` is ``sequence.alignment.args``. This alignment is the coordinate
    system for pairwise identity, the heme-pocket transfer and the pocket
    identity of every structural pair — the highest fan-out artefact in the
    pipeline — so its mode belongs in the manifest rather than inside this
    function (audit 1.3a). ``--anysymbol`` is not configurable: it guards
    against non-standard residue codes in the proteome and changing it would
    make the run fail rather than differ.
    """
    cmd = ["mafft", *[str(a) for a in args],
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
    aln: str | Path,
    prefix: str | Path,
    threads: str | int = "AUTO",
    seed: int = 20240601,
    model: str = "MFP",
    ufboot: int = 1000,
    alrt: int = 1000,
) -> dict[str, str]:
    """IQ-TREE: ModelFinder plus ultrafast-bootstrap and SH-aLRT replicates.

    The tree is a *terminal* branch of the pipeline: it is reported, relabelled
    with chromosomes, and read by no scoring step (``sequence.phylogeny.
    feeds_score: false``). Internal support within the focal clade is 31-47
    because four sequences at >91% identity over 162 columns do not contain
    enough signal to resolve their branching order — more replicates will not
    change that (README §1.3, §9.3).
    """
    exe = iqtree_exe()
    run(
        [exe, "-s", str(aln), "-m", str(model), "-B", str(ufboot),
         "--alrt", str(alrt), "-T", str(threads), "--seed", str(seed),
         "--prefix", str(prefix), "-redo"],
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
    idx = {lab: i for i, lab in enumerate(labels)}
    long_rows = []
    for la, lb in canonical_pair_order(labels):
        i, j = idx[la], idx[lb]
        both = is_res[i] & is_res[j]
        n_both = int(both.sum())
        n_id = int((arr[i][both] == arr[j][both]).sum())
        a = 100.0 * n_id / n_both if n_both else 0.0
        s = 100.0 * n_id / int(min(ungapped_len[i], ungapped_len[j]))
        pid_a[i, j] = pid_a[j, i] = a
        pid_s[i, j] = pid_s[j, i] = s
        long_rows.append(
            {
                "label_a": la, "label_b": lb,
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
    model_name: str,
    device: str | None = None,
    batch_size: int = 8,
    dtype_name: str | None = None,
) -> tuple[list[str], np.ndarray]:
    """Mean-pooled final-layer ESM2 embeddings, over residue tokens only.

    BOS/EOS/PAD are masked before averaging, which matters because sequences of
    different length are batched together.

    ``dtype_name`` ('float32' | 'float16') is explicit rather than inferred from
    the device. The committed cosine matrix is fp32 on CPU; fp16 on a GPU
    differs in the fourth decimal, and leaving that implicit in the hardware is
    how a provenance difference gets mistaken for a result (README §1.5, §5.2).
    """
    import torch
    from transformers import AutoTokenizer, EsmModel

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    if dtype_name is None:
        dtype_name = "float16" if device == "cuda" else "float32"
    if dtype_name not in ("float16", "float32"):
        raise ValueError(f"dtype_name must be float16 or float32, got {dtype_name!r}")
    dtype = torch.float16 if dtype_name == "float16" else torch.float32

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


def parse_gff_genes(
    gff: str | Path,
    id_prefix: str | re.Pattern,
    chrom_regex: str | re.Pattern,
) -> pd.DataFrame:
    """Gene features from the annotation GFF3, with per-seqid rank order.

    ``id_prefix`` and ``chrom_regex`` come from ``family.identifiers``; they
    were module constants compiled against this one assembly's conventions.
    """
    id_prefix = _compiled(id_prefix)
    chrom_re = _compiled(chrom_regex)
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
            gene = attrs.get("Name") or id_prefix.sub("", attrs.get("ID", ""))
            m = chrom_re.search(f[0])
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
    max_intervening: int,
    proximal_max_bp: int,
    identity_long: pd.DataFrame | None = None,
    cosine: pd.DataFrame | None = None,
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
    for a, b in canonical_pair_order(list(by_label.index)):
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
    # Row order is left as canonical_pair_order emitted it — the same order every
    # other module's pair table uses — rather than re-sorted by interest, so the
    # three tables can be eyeballed side by side as well as merged.
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
