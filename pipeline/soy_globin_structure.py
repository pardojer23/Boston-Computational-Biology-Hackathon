"""Module 2 — structure.

Pairwise structural distance over the same labels the sequence module emits,
so the output joins onto ``paralog_pairs.csv`` on ``label_a``/``label_b``.

Three things are deliberate here:

* **Structures come from AlphaFold DB, not a prediction.** Every member of
  this family has a reviewed or curated UniProt entry and therefore an AFDB
  model. Downloading beats predicting; it also means the structures are
  citable. ``NEXT_RUN.md`` budgeted ESMFold on a GPU for this stage — it is
  not needed. Note the AFDB *file* URL pattern (``AF-<acc>-F1-model_v4.pdb``)
  no longer resolves; go through ``/api/prediction/<acc>`` and read the URL
  it hands back.

* **AFDB models are keyed on the UniProt sequence, which is not always the
  Wm82.a4 primary transcript.** Every model is aligned back to the a4 protein
  and the identity and length delta recorded, so a structure standing in for a
  different isoform is visible rather than silent.

* **The heme pocket is taken from a heme-bound crystal structure, not from
  recalled residue numbers.** AFDB models are apo, so there is no ligand to
  measure contacts against. Residues within ``POCKET_CUTOFF_A`` of the heme in
  a soybean leghemoglobin-a crystal are transferred onto every family member
  through the MAFFT alignment the sequence module already produced. UniProt's
  own heme-binding annotations for the same protein are fetched independently
  and used as a check on that transfer, not as its source.

No Modal import, for the same reason ``soy_globin_core`` has none.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import soy_globin_core as core

# --------------------------------------------------------------------------
# constants

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
RCSB_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"

SOYBEAN_TAXON = 3847

#: Heme-bound soybean leghemoglobin-a crystal structure used to define the
#: pocket. 1BIN is 2.20 A (1FSL, 2.30 A, is the alternative); both are
#: cross-referenced from P02238 and span residues 2-144.
POCKET_TEMPLATE_PDB = "1BIN"
POCKET_LIGAND = "HEM"

#: Any residue with a heavy atom within this distance of any heme heavy atom
#: is called part of the pocket.
POCKET_CUTOFF_A = 5.0

_UA = {"User-Agent": "Python-urllib"}


def _get(url: str, tries: int = 3, pause: float = 2.0) -> bytes:
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=_UA)
            with urllib.request.urlopen(req, timeout=60) as fh:
                return fh.read()
        except Exception as exc:  # noqa: BLE001 - surfaced after retries
            last = exc
            time.sleep(pause * (i + 1))
    raise RuntimeError(f"GET failed after {tries} tries: {url}") from last


# --------------------------------------------------------------------------
# 1. gene id -> UniProt accession


def resolve_uniprot(gene_ids: list[str]) -> pd.DataFrame:
    """Map Wm82.a4 gene IDs to UniProt accessions by gene-name search.

    Resolved at run time rather than hardcoded so that pointing the pipeline
    at another family needs no edit here.
    """
    rows = []
    for gid in gene_ids:
        q = f"(gene:{gid}) AND (organism_id:{SOYBEAN_TAXON})"
        url = (
            f"{UNIPROT_SEARCH}?query={urllib.parse.quote(q)}"
            "&fields=accession,id,protein_name,length,reviewed"
            "&format=tsv&size=5"
        )
        txt = _get(url).decode()
        lines = [ln for ln in txt.splitlines() if ln.strip()]
        if len(lines) < 2:
            rows.append(
                {"gene_id": gid, "uniprot": None, "entry_name": None,
                 "uniprot_protein_name": None, "uniprot_length": np.nan,
                 "uniprot_reviewed": None, "n_uniprot_hits": 0}
            )
            continue
        # Prefer a reviewed (Swiss-Prot) entry when the search returns several.
        hits = [ln.split("\t") for ln in lines[1:]]
        hits.sort(key=lambda h: 0 if h[-1].strip().lower().startswith("reviewed") else 1)
        acc, entry, pname, length, reviewed = hits[0][:5]
        rows.append(
            {"gene_id": gid, "uniprot": acc, "entry_name": entry,
             "uniprot_protein_name": pname, "uniprot_length": int(length),
             "uniprot_reviewed": reviewed.strip(), "n_uniprot_hits": len(hits)}
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 2. AlphaFold DB


def fetch_afdb(acc: str, outdir: str | Path) -> dict:
    """Download the AFDB model for ``acc``. Returns provenance, not just a path."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    meta = json.loads(_get(AFDB_API.format(acc=acc)).decode())
    if not meta:
        raise RuntimeError(f"no AlphaFold DB entry for {acc}")
    m = meta[0]
    url = m.get("pdbUrl")
    if not url:
        raise RuntimeError(f"AFDB entry for {acc} has no pdbUrl")
    dest = outdir / f"AF-{acc}-F1.pdb"
    dest.write_bytes(_get(url))
    return {
        "uniprot": acc,
        "pdb_path": str(dest),
        "afdb_entry": m.get("modelEntityId"),
        "afdb_model_created": m.get("modelCreatedDate", "")[:10],
        "afdb_version": m.get("latestVersion"),
        "mean_plddt": m.get("globalMetricValue"),
        "frac_plddt_very_low": m.get("fractionPlddtVeryLow"),
        "afdb_sequence": m.get("uniprotSequence") or m.get("sequence"),
        "pdb_url": url,
    }


# --------------------------------------------------------------------------
# 3. structure I/O and sequence reconciliation

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def read_pdb_chain(path: str | Path, chain: str | None = None):
    """Minimal PDB reader: CA coordinates + one-letter sequence + residue numbers.

    Returns ``(coords (n,3) float64, seq str, resnums list[int])`` for the
    first chain encountered (or ``chain`` if given). Biopython's parser would
    also work; this keeps the dependency surface at numpy and avoids
    Biopython's warnings on AFDB headers.
    """
    coords, seq, resnums = [], [], []
    seen: set[tuple[int, str]] = set()
    want = chain
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break  # first model only
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[16] not in (" ", "A"):
                continue
            one = _THREE_TO_ONE.get(line[17:20].strip())
            if one is None:
                continue
            ch = line[21]
            if want is None:
                want = ch  # lock onto the first chain that yields a residue
            if ch != want:
                continue
            key = (int(line[22:26]), line[26])
            if key in seen:
                continue
            seen.add(key)
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            seq.append(one)
            resnums.append(key[0])
    return np.asarray(coords, dtype=np.float64), "".join(seq), resnums


def read_pdb_heavy_atoms(path: str | Path):
    """All heavy atoms as ``(coords (n,3), keys list[(chain,resnum,resname,record)])``."""
    coords, keys = [], []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            element = line[76:78].strip() or line[12:16].strip()[:1]
            if element == "H":
                continue
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            keys.append((line[21], int(line[22:26]), line[17:20].strip(), line[:6].strip()))
    return np.asarray(coords, dtype=np.float64), keys


def pairwise_identity_global(a: str, b: str) -> tuple[float, int]:
    """Needleman-Wunsch identity (BLOSUM62) over the aligned region.

    Returns ``(percent identity over aligned columns, n aligned columns)``.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    aln = aligner.align(a, b)[0]
    sa, sb = str(aln[0]), str(aln[1])
    both = [(x, y) for x, y in zip(sa, sb) if x != "-" and y != "-"]
    if not both:
        return 0.0, 0
    ident = sum(1 for x, y in both if x == y)
    return 100.0 * ident / len(both), len(both)


# --------------------------------------------------------------------------
# 4. heme pocket definition


def heme_pocket_residues(
    pdb_id: str = POCKET_TEMPLATE_PDB,
    ligand: str = POCKET_LIGAND,
    cutoff: float = POCKET_CUTOFF_A,
    cache_dir: str | Path = ".",
) -> dict:
    """Residues lining the heme in a heme-bound crystal structure.

    Empirical: any polypeptide residue with a heavy atom within ``cutoff`` of
    any heme heavy atom. Returns the crystal chain used, the residue numbers,
    and the one-letter sequence of the chain so callers can map positions.
    """
    cache = Path(cache_dir) / f"{pdb_id}.pdb"
    cache.parent.mkdir(parents=True, exist_ok=True)
    if not cache.exists():
        cache.write_bytes(_get(RCSB_PDB.format(pdb_id=pdb_id)))

    coords, keys = read_pdb_heavy_atoms(cache)
    is_lig = np.array([k[2] == ligand for k in keys])
    if not is_lig.any():
        raise RuntimeError(f"ligand {ligand} not found in {pdb_id}")

    # Use the heme in the first chain that has one, and that chain's protein.
    lig_chain = keys[int(np.flatnonzero(is_lig)[0])][0]
    lig_mask = np.array([k[2] == ligand and k[0] == lig_chain for k in keys])
    prot_mask = np.array(
        [k[3] == "ATOM" and k[0] == lig_chain and k[2] in _THREE_TO_ONE for k in keys]
    )

    lig_xyz = coords[lig_mask]
    prot_xyz = coords[prot_mask]
    prot_keys = [k for k, m in zip(keys, prot_mask) if m]

    d = np.linalg.norm(prot_xyz[:, None, :] - lig_xyz[None, :, :], axis=2)
    near = d.min(axis=1) <= cutoff
    pocket = sorted({prot_keys[i][1] for i in np.flatnonzero(near)})

    _, chain_seq, chain_resnums = read_pdb_chain(cache, chain=lig_chain)
    return {
        "pdb_id": pdb_id,
        "chain": lig_chain,
        "ligand": ligand,
        "cutoff_a": cutoff,
        "pocket_resnums": pocket,
        "chain_seq": chain_seq,
        "chain_resnums": chain_resnums,
        "pdb_path": str(cache),
    }


def uniprot_ligand_sites(acc: str) -> list[dict]:
    """UniProt ``Binding site`` features, as an independent check on the transfer."""
    d = json.loads(_get(UNIPROT_ENTRY.format(acc=acc)).decode())
    out = []
    for f in d.get("features", []):
        if f.get("type") != "Binding site":
            continue
        out.append(
            {
                "position": f["location"]["start"]["value"],
                "ligand": (f.get("ligand") or {}).get("name"),
                "ligand_part": (f.get("ligandPart") or {}).get("name"),
                "description": f.get("description", ""),
            }
        )
    return sorted(out, key=lambda r: r["position"])


def map_positions_through_alignment(
    ref_seq: str, ref_positions: list[int], aln: dict[str, str], ref_label: str
) -> dict[int, int]:
    """Map 1-based positions in ``ref_seq`` to 0-based column indices of ``aln``.

    ``ref_seq`` is aligned to the reference row of the MSA first, so template
    numbering (a crystal chain, or a UniProt sequence) can differ from the
    a4 protein that is actually in the alignment.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligned_ref = aln[ref_label]
    plain_ref = aligned_ref.replace("-", "")

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    a = aligner.align(ref_seq, plain_ref)[0]

    # ref index (0-based) -> plain_ref index (0-based)
    ref_to_plain: dict[int, int] = {}
    ia = ib = 0
    for x, y in zip(str(a[0]), str(a[1])):
        if x != "-" and y != "-":
            ref_to_plain[ia] = ib
        if x != "-":
            ia += 1
        if y != "-":
            ib += 1

    # plain_ref index -> MSA column
    plain_to_col = [i for i, c in enumerate(aligned_ref) if c != "-"]

    out: dict[int, int] = {}
    for p in ref_positions:
        j = ref_to_plain.get(p - 1)
        if j is not None and j < len(plain_to_col):
            out[p] = plain_to_col[j]
    return out
