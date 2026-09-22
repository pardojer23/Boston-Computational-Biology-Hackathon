"""
Resolve CDS sequences for soybean (Glycine max) paralogs resolved by
soybean_expression.py, prune the group's ML tree to those tips, and run a
codon-aware PRANK alignment + HyPhy aBSREL branch-site selection test on
Modal.

General GeneID/mRNA resolution chain (works for any Glyma ID, not a
hardcoded lookup table):
    Glyma.##G###### --(reformat)--> GLYMA_##G######v4
        --(NCBI esearch, db=gene)--> GeneID
        --(NCBI elink gene_nuccore_refseqrna)--> RefSeq mRNA accession
        --(NCBI efetch fasta_cds_na)--> CDS nucleotide sequence
A GeneID with no gene_nuccore_refseqrna link (e.g. an annotated pseudogene)
is skipped with an explicit message rather than failing.

Usage (as a library -- see pipeline.py for the orchestrated call):
    from soybean_selection import run_soybean_selection_analysis
"""
import copy
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from textwrap import wrap

from Bio import Phylo

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def _fetch(url, retries=4):
    """NCBI E-utils without an API key rate-limits at ~3 req/s; back off and
    retry on 429s rather than failing the whole resolution chain."""
    req = urllib.request.Request(url, headers={"User-Agent": "research-agent"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                time.sleep(0.34)  # stay under ~3 req/s
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise


def translate(seq: str) -> str:
    codons = wrap(seq, 3)
    aa = "".join(CODON_TABLE.get(c, "X") for c in codons if len(c) == 3)
    return aa.rstrip("*")


def resolve_geneid_from_glyma(glyma_id: str, organism="Glycine max"):
    core = glyma_id.replace("Glyma.", "")
    tag = f"GLYMA_{core}v4"
    term = f'{tag}[All Fields] AND "{organism}"[Organism]'
    url = f"{EUTILS}/esearch.fcgi?db=gene&term={urllib.parse.quote(term)}&retmode=json"
    data = json.loads(_fetch(url))
    ids = data["esearchresult"]["idlist"]
    return ids[0] if ids else None


def resolve_mrna_from_geneid(gene_id: str):
    url = f"{EUTILS}/elink.fcgi?dbfrom=gene&db=nuccore&id={gene_id}&linkname=gene_nuccore_refseqrna&retmode=json"
    data = json.loads(_fetch(url))
    linksetdbs = data["linksets"][0].get("linksetdbs")
    if not linksetdbs:
        return None  # no mRNA -- e.g. an annotated pseudogene
    nuccore_uid = linksetdbs[0]["links"][0]

    url2 = f"{EUTILS}/esummary.fcgi?db=nuccore&id={nuccore_uid}&retmode=json"
    data2 = json.loads(_fetch(url2))
    return data2["result"][nuccore_uid]["accessionversion"]


def fetch_cds(mrna_accession: str):
    url = f"{EUTILS}/efetch.fcgi?db=nuccore&id={mrna_accession}&rettype=fasta_cds_na&retmode=text"
    text = _fetch(url).decode()
    lines = text.strip().split("\n")
    return "".join(lines[1:])


def resolve_cds_for_glyma_genes(resolved_records: list) -> list:
    """For each Glyma-resolved record, walk Glyma ID -> GeneID -> mRNA -> CDS,
    validating the CDS translation against the OrthoDB protein sequence.
    Records with no resolvable mRNA (e.g. pseudogenes) get cds_seq=None with
    an explanatory note; the caller is responsible for excluding those from
    the selection test."""
    out = []
    for r in resolved_records:
        glyma_id = r.get("glyma_id")
        if not glyma_id:
            out.append({**r, "cds_seq": None, "cds_note": "no Glyma ID resolved"})
            continue

        gene_id = resolve_geneid_from_glyma(glyma_id)
        if gene_id is None:
            out.append({**r, "cds_seq": None, "cds_note": f"no NCBI GeneID found for {glyma_id}"})
            continue

        mrna_acc = resolve_mrna_from_geneid(gene_id)
        if mrna_acc is None:
            out.append({**r, "ncbi_gene_id": gene_id, "cds_seq": None,
                        "cds_note": "no RefSeq mRNA linked to this GeneID (likely a pseudogene)"})
            continue

        cds_seq = fetch_cds(mrna_acc)
        translated = translate(cds_seq)
        matches_orthodb = translated == r["seq"]
        out.append({**r, "ncbi_gene_id": gene_id, "mrna_accession": mrna_acc,
                    "cds_seq": cds_seq, "cds_len_nt": len(cds_seq),
                    "translation_matches_orthodb": matches_orthodb, "cds_note": None})
    return out


def tree_label(odb_id: str, organism: str) -> str:
    org = (organism or "unknown").replace(" ", "_").replace(".", "")
    return f"{org}__{odb_id.replace(':', '_')}"


def prune_tree_to_labels(full_tree_newick: str, target_labels: list) -> str:
    import io
    full_tree = Phylo.read(io.StringIO(full_tree_newick), "newick")
    all_tips = [t.name for t in full_tree.get_terminals()]
    missing = [t for t in target_labels if t not in all_tips]
    if missing:
        raise ValueError(f"Target tips not found in tree: {missing}")

    pruned = copy.deepcopy(full_tree)
    for tip_name in all_tips:
        if tip_name not in target_labels:
            clade = next(pruned.find_clades(name=tip_name), None)
            if clade is not None:
                pruned.prune(clade)

    buf = io.StringIO()
    Phylo.write(pruned, buf, "newick")
    return buf.getvalue()


def run_prank_hyphy_modal(cds_fasta_text: str, tree_newick_text: str, modal_python: str = None) -> dict:
    """Dispatch the codon-aware alignment + branch-site test to Modal.

    Runs as a SUBPROCESS against a separate Python interpreter that has the
    `modal` SDK installed and authenticated (typically NOT the same
    environment as the scientific stack running this orchestrator -- this
    project's own Modal SDK lives in a dedicated conda env). Pass
    `modal_python` explicitly, or set the MODAL_PYTHON_BIN environment
    variable; defaults to plain `python` on PATH (works if the current
    environment already has `modal` importable)."""
    import json
    import os
    import subprocess
    import tempfile

    modal_python = modal_python or os.environ.get("MODAL_PYTHON_BIN", "python")
    this_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(this_dir, "modal_absrel.py")

    with tempfile.TemporaryDirectory() as tmp:
        cds_path = os.path.join(tmp, "cds.fasta")
        tree_path = os.path.join(tmp, "tree.nwk")
        out_path = os.path.join(tmp, "result.json")
        with open(cds_path, "w") as f:
            f.write(cds_fasta_text)
        with open(tree_path, "w") as f:
            f.write(tree_newick_text)

        proc = subprocess.run([modal_python, script_path, cds_path, tree_path, out_path],
                               capture_output=True, text=True, timeout=1800)
        if proc.returncode != 0 or not os.path.exists(out_path):
            raise RuntimeError(f"Modal dispatch failed (returncode={proc.returncode}):\n"
                                f"STDOUT: {proc.stdout[-2000:]}\nSTDERR: {proc.stderr[-2000:]}")
        with open(out_path) as f:
            return json.load(f)


def run_soybean_selection_analysis(resolved_records: list, full_tree_newick: str, organism="Glycine max") -> dict:
    """Orchestrated entry point. `resolved_records` should be the
    `resolved_records` list from soybean_expression.run_soybean_expression_analysis
    (records with a `glyma_id`, possibly None). Skips cleanly if fewer than 2
    resolvable CDS sequences are available (aBSREL needs >=2 taxa; a
    meaningful branch-site comparison needs >=3)."""
    with_cds = resolve_cds_for_glyma_genes([r for r in resolved_records if r.get("glyma_id")])
    valid = [r for r in with_cds if r.get("cds_seq")]

    if len(valid) < 3:
        return {"skipped": True, "reason": f"only {len(valid)} paralog(s) with a resolvable CDS "
                                            f"(need >=3 for a meaningful branch-site tree)",
                "cds_resolution_detail": with_cds}

    target_labels = [tree_label(r["odb_id"], organism) for r in valid]
    pruned_newick = prune_tree_to_labels(full_tree_newick, target_labels)

    cds_fasta = "\n".join(f">{tree_label(r['odb_id'], organism)}\n{r['cds_seq']}" for r in valid)

    modal_result = run_prank_hyphy_modal(cds_fasta, pruned_newick)

    return {
        "skipped": False,
        "n_taxa": len(valid),
        "cds_resolution_detail": with_cds,
        "pruned_tree_newick": pruned_newick,
        "cds_fasta": cds_fasta,
        "modal_result": modal_result,
    }
