#!/usr/bin/env python
"""One invocable entry point per pipeline step.

    python pipeline/steps.py <step> [--config config/config.yaml] [paths...]

Each subcommand is a thin shim: parse declared inputs, outputs and params, call
the unchanged library function, validate the result against its contract, and
write a provenance sidecar. No science lives here.

This replaces the four hand-ordered runners (``run_local.py``,
``run_structure.py``, ``soy_globin_expression.py --main``,
``run_integration.py``). Those decided the order themselves, in prose in
README §7, with no way to tell whether a downstream output was older than the
input it came from. The order now lives in ``workflow/Snakefile`` as file
dependencies, and each step below is one rule.

Provenance chaining
-------------------
Every step writes ``<output_dir>/.prov/<step>.json`` recording the sha256 of
each input it consumed and each output it produced, plus tool versions and the
config digest. ``run-manifest`` collects them into
``results/run_manifest.json``. That is what makes a partial re-run detectable:
previously there were four independent manifests, none of which recorded the
identity of the upstream outputs behind it, so a module-4 table computed
against a stale module-1 pair table looked exactly like a good one
(docs/PIPELINE_AUDIT.md findings 1 and 2).
"""

from __future__ import annotations

import argparse
import json
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

import config as cfgmod                  # noqa: E402
import contracts                         # noqa: E402
import soy_globin_core as core           # noqa: E402


# --------------------------------------------------------------------------- #
# Provenance
# --------------------------------------------------------------------------- #


def _digest_files(paths) -> dict[str, str]:
    out = {}
    for p in paths:
        p = Path(p)
        if p.exists() and p.is_file():
            out[str(p)] = core.sha256(p)
    return out


def write_prov(
    step: str,
    cfg: cfgmod.Config,
    inputs,
    outputs,
    params: dict | None = None,
    tools: dict | None = None,
    extra: dict | None = None,
    prov_dir: str | Path = "results/.prov",
) -> Path:
    """Record what this step consumed and produced.

    The input digests are the load-bearing part: they are what a later step (or
    a reader) can compare against the *current* digest of the same file to see
    that a re-run is needed. Snakemake's own staleness check is by mtime, which
    a ``touch`` or a checkout defeats.
    """
    prov_dir = Path(prov_dir)
    prov_dir.mkdir(parents=True, exist_ok=True)
    rec = {
        "step": step,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config_digest": cfg.digest(),
        "config_source": str(cfg.source) if cfg.source else None,
        "params": params or {},
        "tools": tools or {},
        "python": platform.python_version(),
        "inputs_sha256": _digest_files(inputs),
        "outputs_sha256": _digest_files(outputs),
    }
    if extra:
        rec.update(extra)
    p = prov_dir / f"{step}.json"
    p.write_text(json.dumps(rec, indent=2, sort_keys=True, default=str))
    return p


def _pkg_versions(*names: str) -> dict[str, str]:
    out = {}
    for n in names:
        try:
            mod = __import__(n)
            out[n] = getattr(mod, "__version__", "?")
        except Exception:
            out[n] = "not installed"
    return out


# --------------------------------------------------------------------------- #
# Stage 0 — reference data
# --------------------------------------------------------------------------- #


def step_fetch_reference(a, cfg):
    """Download and verify one pinned reference file."""
    if a.which == "hmm":
        spec = cfg.section("sequence")["hmm"]
        rec = core.fetch_file(spec["url"], a.out + ".gz", spec["sha256"],
                              force=a.force, gunzip_to=a.out)
    else:
        spec = cfg[f"reference.files.{a.which}"]
        rec = core.fetch_file(spec["url"], a.out, spec["sha256"], force=a.force)

    checks: dict = {}
    if a.which == "proteome":
        n = len(core.read_fasta(a.out))
        cfgmod.evaluate_check(
            cfg, "reference_checksums", True,
            {"file": a.which, "sha256": rec["sha256"]}, sink=checks)
        want = spec.get("n_sequences_expected")
        if want:
            cfgmod.evaluate_check(
                cfg, "reference_checksums", n == int(want),
                {"file": a.which, "n_sequences": n, "expected": int(want)},
                sink=checks)
    if a.which == "hmm":
        meta = core.hmm_metadata(a.out)
        cfgmod.evaluate_check(
            cfg, "hmm_release", meta.get("acc") == spec.get("release_expected"),
            {"acc": meta.get("acc"), "expected": spec.get("release_expected"),
             "leng": meta.get("leng"), "ga": meta.get("ga")}, sink=checks)
        checks["hmm_metadata"] = meta

    write_prov(f"fetch_reference_{a.which}", cfg, [], [a.out],
               params={"which": a.which, "url": spec["url"]},
               extra={"fetch": rec, "checks": checks}, prov_dir=a.prov)
    print(f"[fetch_reference] {a.which} -> {a.out}  sha256={rec['sha256'][:12]}")


# --------------------------------------------------------------------------- #
# Stage 1 — sequence
# --------------------------------------------------------------------------- #


def step_select_family(a, cfg):
    """HMM search over the proteome, unioned with the configured focal genes."""
    seqs, members = core.select_family(
        a.proteome, a.hmm, a.workdir, cfg, cpu=a.threads)
    core.write_fasta(seqs, a.out_faa)
    contracts.write_csv(members, a.out_members,
                        contracts.members_spec(cfg, len(members)))
    n_searched = len(core.read_fasta(a.proteome))
    write_prov("select_family", cfg, [a.proteome, a.hmm],
               [a.out_faa, a.out_members],
               params={"threshold": cfg["sequence.select.threshold"],
                       "accession": cfg["sequence.hmm.accession"],
                       "cfg_section": cfg.section_digest("sequence.select",
                                                         "family.focal_genes")},
               tools={"hmmsearch": core.tool_version(["hmmsearch", "-h"])},
               extra={"n_proteins_searched": n_searched,
                      "n_members": int(len(members)),
                      "n_hits": int(members.pfam_hit.sum()),
                      "n_focal_rescued": int((~members.pfam_hit).sum())},
               prov_dir=a.prov)
    print(f"[select_family] {len(members)} of {n_searched} proteins")


def step_cutoff_check(a, cfg):
    """Relaxed-threshold search: evidence the family is not truncated."""
    members = pd.read_csv(a.members)
    n_hits = int(members.pfam_hit.sum())
    ev = core.cutoff_gap_check(
        a.hmm, a.proteome, a.workdir, n_accepted=n_hits,
        relaxed_evalue=float(cfg["sequence.cutoff_check.relaxed_evalue"]),
        cpu=a.threads)
    checks: dict = {}
    cfgmod.evaluate_check(
        cfg, "family_cutoff_gap",
        ev["bitscore_gap"] >= float(cfg["sequence.cutoff_check.min_bitscore_gap"]),
        {**ev, "min_required": float(cfg["sequence.cutoff_check.min_bitscore_gap"])},
        sink=checks)
    Path(a.out).write_text(json.dumps(checks, indent=2, default=str))
    write_prov("cutoff_check", cfg, [a.hmm, a.proteome, a.members], [a.out],
               params={"relaxed_evalue": cfg["sequence.cutoff_check.relaxed_evalue"]},
               extra={"checks": checks}, prov_dir=a.prov)
    print(f"[cutoff_check] bit-score gap {ev['bitscore_gap']:.1f} "
          f"({ev['weakest_accepted_bitscore']:.1f} accepted vs "
          f"{ev['best_rejected_bitscore']} rejected)")


def step_align(a, cfg):
    """MAFFT with the configured mode."""
    core.run_mafft(a.faa, a.out, threads=a.threads,
                   args=cfg["sequence.alignment.args"])
    aln = core.read_fasta(a.out)
    ncol = len(next(iter(aln.values())))
    checks: dict = {}
    want = cfg.get("sequence.alignment.n_columns_expected")
    if want:
        cfgmod.evaluate_check(cfg, "alignment_columns", ncol == int(want),
                              {"n_columns": ncol, "expected": int(want)},
                              sink=checks)
    write_prov("align", cfg, [a.faa], [a.out],
               params={"args": list(cfg["sequence.alignment.args"]),
                       "cfg_section": cfg.section_digest("sequence.alignment")},
               tools={"mafft": core.tool_version(["mafft", "--version"])},
               extra={"n_sequences": len(aln), "n_columns": ncol,
                      "checks": checks}, prov_dir=a.prov)
    print(f"[align] {len(aln)} sequences, {ncol} columns")


def step_tree(a, cfg):
    """IQ-TREE. Terminal branch: no scoring step reads the tree."""
    p = cfg.section("sequence")["phylogeny"]
    tree = core.run_iqtree(a.aln, a.prefix, threads=str(a.threads),
                           seed=int(p["seed"]), model=p["model"],
                           ufboot=int(p["ufboot"]), alrt=int(p["alrt"]))
    shutil.copy(tree["treefile"], a.out_treefile)
    shutil.copy(tree["iqtree_report"], a.out_report)
    write_prov("tree", cfg, [a.aln], [a.out_treefile, a.out_report],
               params={**p, "cfg_section": cfg.section_digest("sequence.phylogeny")},
               tools={"iqtree": core.tool_version([core.iqtree_exe(), "--version"])},
               extra={"best_model": core.best_model(tree["iqtree_report"]),
                      "feeds_score": bool(p.get("feeds_score", False))},
               prov_dir=a.prov)
    print(f"[tree] best-fit model {core.best_model(tree['iqtree_report'])}")


def step_identity(a, cfg):
    """Pairwise % identity from the MSA, so it is consistent with the tree."""
    mat, long = core.pairwise_identity(a.aln)
    mat.to_csv(a.out_matrix)
    long.to_csv(a.out_pairs, index=False)
    write_prov("identity", cfg, [a.aln], [a.out_matrix, a.out_pairs],
               extra={"n_pairs": int(len(long))}, prov_dir=a.prov)
    print(f"[identity] {len(long)} pairs, "
          f"{long.pid_aligned.min():.1f}-{long.pid_aligned.max():.1f}% identity")


def step_embed(a, cfg):
    """ESM2 embeddings. Runs in the esm2 environment, not the analysis one.

    This is the one genuine environment seam in the pipeline. The old runner
    handled it with three mutually exclusive flags (--skip-esm2 / --reuse-esm2
    / compute), and the reuse path recorded provenance as the sentence "reused
    from a previous run (see that run's manifest)". Here it is a rule with its
    own declared conda environment and a real output digest; a missing input is
    the workflow's problem, not a flag.
    """
    e = cfg.section("sequence")["embedding"]
    seqs = core.read_fasta(a.faa)
    labels, emb = core.esm2_embeddings(
        seqs, model_name=e["model"], batch_size=int(e["batch_size"]),
        device=a.device, dtype_name=e["dtype"])
    cos = core.cosine_distance_matrix(labels, emb)
    np.savez_compressed(a.out_npz, labels=np.array(labels), embeddings=emb)
    cos.to_csv(a.out_matrix)
    write_prov("embed", cfg, [a.faa], [a.out_matrix, a.out_npz],
               params={**e, "device": a.device or "auto",
                       "cfg_section": cfg.section_digest("sequence.embedding")},
               tools=_pkg_versions("torch", "transformers"),
               extra={"embedding_dim": int(emb.shape[1]),
                      "pooling": e["pooling"], "n_sequences": len(labels)},
               prov_dir=a.prov)
    print(f"[embed] {len(labels)} sequences, dim {emb.shape[1]}, dtype {e['dtype']}")


def step_classify_pairs(a, cfg):
    """Genomic context and duplication mode; emits the sequence pair table."""
    d = cfg.section("sequence")["duplication"]
    members = pd.read_csv(a.members)
    genes = core.parse_gff_genes(
        a.gff3, cfg["family.identifiers.id_prefix_strip"],
        cfg["family.identifiers.chromosome_regex"])
    ident_long = pd.read_csv(a.identity_pairs)
    cos = pd.read_csv(a.cosine, index_col=0)
    pairs, ctx = core.classify_pairs(
        members, genes,
        max_intervening=int(d["tandem_max_intervening_genes"]),
        proximal_max_bp=int(d["proximal_max_bp"]),
        identity_long=ident_long, cosine=cos)
    n = len(members)
    contracts.write_csv(pairs, a.out_pairs, contracts.paralog_pairs_spec(cfg, n))
    contracts.write_csv(ctx, a.out_context, contracts.gene_context_spec(cfg, n))
    write_prov("classify_pairs", cfg,
               [a.members, a.gff3, a.identity_pairs, a.cosine],
               [a.out_pairs, a.out_context],
               params={**d, "cfg_section": cfg.section_digest("sequence.duplication")},
               extra={"counts": pairs.duplication_mode.value_counts().to_dict()},
               prov_dir=a.prov)
    print(f"[classify_pairs] {len(pairs)} pairs: "
          f"{pairs.duplication_mode.value_counts().to_dict()}")


def step_annotate_tree(a, cfg):
    """Relabel tree tips with chromosome context. Presentation only."""
    ctx = pd.read_csv(a.context)
    core.annotate_newick(a.treefile, a.out, ctx)
    write_prov("annotate_tree", cfg, [a.treefile, a.context], [a.out],
               prov_dir=a.prov)
    print(f"[annotate_tree] -> {a.out}")


# --------------------------------------------------------------------------- #
# Stage 2 — structure
# --------------------------------------------------------------------------- #


def step_uniprot_map(a, cfg):
    """Persist gene -> UniProt accession with the release stamped on it."""
    import soy_globin_structure as struc

    members = pd.read_csv(a.members)
    up = struc.resolve_uniprot(members.gene_id.tolist(), cfg)
    up = up.merge(members[["gene_id", "label"]], on="gene_id")
    missing = up[up.uniprot.isna()].gene_id.tolist()
    if missing:
        raise SystemExit(f"no UniProt accession for: {missing}")
    contracts.write_csv(up, a.out,
                        contracts.uniprot_accessions_spec(cfg, len(members)))
    write_prov("uniprot_map", cfg, [a.members], [a.out],
               params={"taxon": cfg["family.species.ncbi_taxon"],
                       "pinned": bool(cfg.get("structure.uniprot.pinned"))},
               extra={"uniprot_release": up.uniprot_release.iloc[0],
                      "accessions": dict(zip(up.gene_id, up.uniprot))},
               prov_dir=a.prov)
    print(f"[uniprot_map] release {up.uniprot_release.iloc[0]}: "
          f"{dict(zip(up.label, up.uniprot))}")


def step_fetch_afdb(a, cfg):
    """Download one AlphaFold DB model into the cache (not into results)."""
    import soy_globin_structure as struc

    rec = struc.fetch_afdb(a.acc, Path(a.out).parent,
                           api=cfg["structure.afdb.api"], force=a.force)
    meta = {k: v for k, v in rec.items() if k != "afdb_sequence"}
    Path(a.out_meta).write_text(json.dumps(rec, indent=2, default=str))
    checks: dict = {}
    want = cfg.get("structure.afdb.model_version_expected")
    if want is not None:
        cfgmod.evaluate_check(cfg, "afdb_model_version",
                              str(rec.get("afdb_version")) == str(want),
                              {"acc": a.acc, "version": rec.get("afdb_version"),
                               "expected": want}, sink=checks)
    write_prov(f"fetch_afdb_{a.acc}", cfg, [], [a.out, a.out_meta],
               params={"acc": a.acc, "api": cfg["structure.afdb.api"]},
               extra={"afdb": meta, "checks": checks}, prov_dir=a.prov)
    print(f"[fetch_afdb] {a.acc} {rec['afdb_entry']} pLDDT {rec['mean_plddt']:.1f}")


def step_pocket_template(a, cfg):
    """Ligand-contact residues from the crystal template.

    A DAG root: depends on (template, ligand, cutoff) and nothing else, so the
    workflow caches it on those three values.
    """
    import soy_globin_structure as struc

    p = cfg.section("structure")["pocket"]
    pk = struc.heme_pocket_residues(
        pdb_id=p["template_pdb"], ligand=p["ligand"],
        cutoff=float(p["cutoff_a"]), cache_dir=Path(a.out).parent,
        template_api=cfg["structure.pocket.template_api"])
    checks: dict = {}
    want = p.get("n_residues_expected")
    if want:
        cfgmod.evaluate_check(
            cfg, "pocket_residue_count", len(pk["pocket_resnums"]) == int(want),
            {"n_residues": len(pk["pocket_resnums"]), "expected": int(want),
             "chain": pk["chain"], "resnums": pk["pocket_resnums"]}, sink=checks)
    Path(a.out).write_text(json.dumps({**pk, "checks": checks}, indent=2, default=str))
    write_prov("pocket_template", cfg, [], [a.out],
               params={"template_pdb": p["template_pdb"], "ligand": p["ligand"],
                       "cutoff_a": p["cutoff_a"],
                       "cfg_section": cfg.section_digest("structure.pocket")},
               extra={"chain": pk["chain"],
                      "n_pocket_residues": len(pk["pocket_resnums"])},
               prov_dir=a.prov)
    print(f"[pocket_template] {p['template_pdb']} chain {pk['chain']}: "
          f"{len(pk['pocket_resnums'])} residues")


def step_structures(a, cfg):
    """Reconcile each AFDB model against the reference transcript."""
    import soy_globin_structure as struc

    up = pd.read_csv(a.uniprot)
    prot = core.read_fasta(a.faa)
    metas = [json.loads(Path(p).read_text()) for p in a.afdb_meta]
    prov = pd.DataFrame(metas).merge(up, left_on="uniprot", right_on="uniprot")

    out_pdb = Path(a.out_pdb_dir)
    out_pdb.mkdir(parents=True, exist_ok=True)
    rec, checks = [], {}
    mvr = cfg.section("structure")["model_vs_reference"]
    max_delta = int(mvr["max_len_delta"])
    min_pid = float(mvr["min_identity"])
    known = mvr.get("known_exceptions") or {}
    for r in prov.itertuples():
        dest = out_pdb / Path(r.pdb_path).name
        shutil.copy(r.pdb_path, dest)
        _, struct_seq, _ = struc.read_pdb_chain(dest)
        a4 = prot[r.label]
        pid, ncol = struc.pairwise_identity_global(struct_seq, a4)
        rec.append({"label": r.label, "n_res_model": len(struct_seq),
                    "len_a4_protein": len(a4),
                    "len_delta": len(struct_seq) - len(a4),
                    "pid_model_vs_a4": round(pid, 2), "n_aligned_cols": ncol,
                    "model_is_a4_sequence": bool(struct_seq == a4),
                    "pdb_path": str(dest)})
    rec = pd.DataFrame(rec)
    # A model whose sequence is not the reference transcript is not necessarily
    # wrong — AFDB is keyed on the UniProt sequence — but it has to be declared.
    # Rows named in `known_exceptions` are exempted by gene ID, with the reason
    # carried into the check evidence rather than left in a code comment.
    rec["gene_id"] = [core.gene_from_label(l) for l in rec.label]
    off = rec[(rec.len_delta.abs() > max_delta) | (rec.pid_model_vs_a4 < min_pid)]
    unexplained = off[~off.gene_id.isin(known)]
    cfgmod.evaluate_check(
        cfg, "model_vs_reference_sequence", bool(unexplained.empty),
        {"max_len_delta": max_delta, "min_identity": min_pid,
         "differing": rec.loc[rec.len_delta != 0,
                              ["label", "len_delta", "pid_model_vs_a4"]]
                         .to_dict("records"),
         "exempted": {g: known[g] for g in off.gene_id if g in known},
         "unexplained": unexplained[["label", "len_delta", "pid_model_vs_a4"]]
                        .to_dict("records")},
        sink=checks)

    structures = (prov.merge(rec.drop(columns=["pdb_path", "gene_id"]), on="label")
                      .sort_values("label")
                      .reset_index(drop=True))
    # Sorted on label so row order is a property of the family, not of the
    # order two frames happened to merge in. The pre-refactor table also
    # carried `uniprot_x`/`uniprot_y` merge suffixes because the accession
    # column arrived from both sides; the accession table is now a declared
    # input with one `uniprot` column.
    structures.drop(columns=[c for c in ("afdb_sequence",) if c in structures],
                    errors="ignore").to_csv(a.out, index=False)
    write_prov("structures", cfg, [a.uniprot, a.faa, *a.afdb_meta], [a.out],
               params={"max_len_delta": max_delta, "min_identity": min_pid},
               extra={"checks": checks}, prov_dir=a.prov)
    print(f"[structures] {len(structures)} models; "
          f"{int((rec.len_delta != 0).sum())} differ in length from the reference")


def step_pocket_residues(a, cfg):
    """Transfer crystal pocket numbering onto the MSA columns."""
    import soy_globin_structure as struc

    pk = json.loads(Path(a.pocket).read_text())
    aln = core.read_fasta(a.aln)
    prot = core.read_fasta(a.faa)
    up = pd.read_csv(a.uniprot)
    labels = list(prot)

    # The reference member is configured, not matched by symbol suffix.
    ref_gene = cfg.pocket_reference_gene
    ref_label = cfg.label_of(ref_gene)
    if ref_label not in aln:
        raise SystemExit(
            f"pocket reference member {ref_label!r} (gene {ref_gene}) is not in "
            f"the alignment; members are {labels}")

    pos2col = struc.map_positions_through_alignment(
        pk["chain_seq"], pk["pocket_resnums"], aln, ref_label)
    dropped = [p for p in pk["pocket_resnums"] if p not in pos2col]
    cols = sorted(pos2col.values())

    # Independent cross-check: UniProt's own binding-site annotations.
    ref_acc = up.loc[up.gene_id == ref_gene, "uniprot"].iloc[0]
    sites = struc.uniprot_ligand_sites(ref_acc, cfg["structure.uniprot.entry_api"])
    site_cols = struc.map_positions_through_alignment(
        prot[ref_label], [s["position"] for s in sites], aln, ref_label)
    n_inside = sum(1 for s in sites if site_cols.get(s["position"]) in cols)
    checks: dict = {}
    cfgmod.evaluate_check(
        cfg, "pocket_uniprot_crosscheck",
        n_inside >= int(cfg.get("structure.pocket_crosscheck.min_sites_inside", 0)),
        {"accession": ref_acc, "n_sites": len(sites), "n_inside_pocket": n_inside,
         "sites": sites}, sink=checks)

    rows = []
    for lab in labels:
        row = {"label": lab, "gene_id": core.gene_from_label(lab)}
        for p, col in sorted(pos2col.items()):
            row[f"col{col + 1}_ref{p}"] = aln[lab][col]
        rows.append(row)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    Path(a.out_cols).write_text(json.dumps(
        {"reference_gene": ref_gene, "reference_label": ref_label,
         "alignment_columns": cols, "pos_to_col": {str(k): v for k, v in pos2col.items()},
         "unmappable_resnums": dropped, "checks": checks}, indent=2, default=str))
    write_prov("pocket_residues", cfg, [a.pocket, a.aln, a.faa, a.uniprot],
               [a.out, a.out_cols],
               params={"reference_member": cfg["structure.pocket.reference_member"]},
               extra={"n_alignment_columns": len(cols),
                      "unmappable_resnums": dropped, "checks": checks},
               prov_dir=a.prov)
    print(f"[pocket_residues] {len(cols)} alignment columns; "
          f"{n_inside}/{len(sites)} UniProt sites inside the pocket")


def step_structure_pairs(a, cfg):
    """TM-align all-vs-all plus pocket identity over the mapped columns."""
    import soy_globin_structure as struc
    from tmtools import tm_align

    aln = core.read_fasta(a.aln)
    cols = json.loads(Path(a.pocket_cols).read_text())["alignment_columns"]
    structures = pd.read_csv(a.structures)
    pdb_dir = Path(a.pdb_dir)
    labels = list(structures.label)

    chains = {}
    for r in structures.itertuples():
        path = pdb_dir / f"AF-{r.uniprot}-F1.pdb"
        c, s, _ = struc.read_pdb_chain(path)
        chains[r.label] = (c, s)

    rows = []
    for x, y in core.canonical_pair_order(labels):
        ca, sa = chains[x]
        cb, sb = chains[y]
        res = tm_align(ca, cb, sa, sb)
        # Normalising by the shorter chain is the conservative choice: it cannot
        # be inflated by one protein being a fragment of the other.
        tm1, tm2 = float(res.tm_norm_chain1), float(res.tm_norm_chain2)
        na, nb = len(sa), len(sb)
        pa, pb = [aln[x][c] for c in cols], [aln[y][c] for c in cols]
        both = [(p, q) for p, q in zip(pa, pb) if p != "-" and q != "-"]
        pocket_id = (100.0 * sum(1 for p, q in both if p == q) / len(both)
                     if both else np.nan)
        rows.append({
            "label_a": x, "label_b": y,
            "gene_a": core.gene_from_label(x), "gene_b": core.gene_from_label(y),
            "tm_score": round(tm1 if na <= nb else tm2, 4),
            "tm_norm_a": round(tm1, 4), "tm_norm_b": round(tm2, 4),
            "rmsd": round(float(res.rmsd), 3),
            "pocket_identity": round(pocket_id, 2),
            "n_pocket_cols": len(both), "n_res_a": na, "n_res_b": nb,
        })
    pairs = pd.DataFrame(rows)
    contracts.write_csv(pairs, a.out,
                        contracts.structure_pairs_spec(cfg, len(labels)))
    write_prov("structure_pairs", cfg, [a.aln, a.pocket_cols, a.structures],
               [a.out],
               params={"tm_score_convention": "normalised by the shorter chain",
                       "n_pocket_cols": len(cols)},
               tools=_pkg_versions("tmtools"),
               extra={"tm_score_range": [float(pairs.tm_score.min()),
                                         float(pairs.tm_score.max())]},
               prov_dir=a.prov)
    print(f"[structure_pairs] {len(pairs)} pairs, TM-score "
          f"{pairs.tm_score.min():.4f}-{pairs.tm_score.max():.4f}")


# --------------------------------------------------------------------------- #
# Stage 3 — expression
# --------------------------------------------------------------------------- #


def step_fetch_expression(a, cfg):
    """Download the three CellRanger parts of one library."""
    import soy_globin_expression as expr

    rec = expr.fetch_library(cfg, a.datadir, a.gsm, force=a.force)
    write_prov(f"fetch_expression_{a.gsm}", cfg, [],
               [p["path"] for p in rec["parts"].values()],
               params={"gsm": a.gsm, "stem": rec["stem"], "tissue": rec["tissue"]},
               extra={"library": rec}, prov_dir=a.prov)
    total = sum(p["bytes"] for p in rec["parts"].values())
    print(f"[fetch_expression] {a.gsm} ({rec['tissue']}) {total / 1e6:.0f} MB")


def step_pseudobulk(a, cfg):
    """Sum counts over every barcode in each library. No cell calling."""
    import soy_globin_expression as expr

    counts, cpm = expr.build_pseudobulk(cfg, a.datadir, cache=a.out,
                                        force=a.force)
    write_prov("pseudobulk", cfg, sorted(expr.input_checksums(cfg, a.datadir)),
               [a.out],
               params={"normalisation": cfg["expression.pseudobulk.normalisation"],
                       "cache_key": "sha256 of every matrix and feature file"},
               extra={"library_totals": {c: int(counts[c].sum())
                                         for c in counts.columns},
                      "n_features": int(counts.shape[0])},
               prov_dir=a.prov)
    print(f"[pseudobulk] {counts.shape[0]} features x {counts.shape[1]} libraries")


def step_expression(a, cfg):
    """Per-gene profiles, pair metrics and the marker panel."""
    import soy_globin_expression as expr

    r = expr.run(cfg, a.datadir, a.outdir, members_csv=a.members, fetch=False)
    outs = [Path(a.outdir) / f for f in ("gene_pseudobulk_profiles.csv",
                                         "expression_pairs.csv",
                                         "marker_profiles.csv", "manifest.json")]
    write_prov("expression", cfg, [a.members, a.counts_cache], outs,
               params={"resolution": cfg["expression.resolution"],
                       "cfg_section": cfg.section_digest("expression")},
               extra={"n_pairs": int(len(r["pairs"])),
                      "markers_absent": r["manifest"]["markers_absent_from_series"]},
               prov_dir=a.prov)
    print(f"[expression] {len(r['profiles'])} genes, {len(r['pairs'])} pairs")


# --------------------------------------------------------------------------- #
# Stage 4 — integration
# --------------------------------------------------------------------------- #


def step_integration(a, cfg):
    """Join the three pair tables, score, sweep alpha, run the checks."""
    import soy_globin_integration as integ

    r = integ.run(cfg, a.results_dir, a.outdir, alpha=a.alpha)
    outs = [Path(a.outdir) / f for f in ("redundancy_scores.csv",
                                         "weight_sensitivity.csv", "manifest.json")]
    ins = [Path(a.results_dir) / p for p in (
        "sequence_module/paralog_pairs.csv", "structure_module/structure_pairs.csv",
        "expression_module/expression_pairs.csv",
        "expression_module/gene_pseudobulk_profiles.csv")]
    write_prov("integration", cfg, ins, outs,
               params={"alpha_M": r["scores"].alpha_M.iloc[0],
                       "beta_E": r["scores"].beta_E.iloc[0],
                       "m_axis_weights": cfg.m_axis_weights,
                       "clade_floor": cfg["score.normalisation.clade_floor"],
                       "cfg_section": cfg.section_digest("score")},
               extra={"validation_checks": r["manifest"]["validation_checks"]},
               prov_dir=a.prov)
    top = r["scores"].iloc[0]
    print(f"[integration] {len(r['scores'])} pairs; top "
          f"{top.label_a}|{top.label_b} R_family={top.R_family:.3f}")


def step_run_manifest(a, cfg):
    """Collect every step's provenance sidecar into one run manifest.

    This is what the four independent module manifests could not do: state,
    for each step, the digests of the exact upstream outputs it consumed. A
    partial re-run shows up as an input digest that no longer matches the
    current file.
    """
    prov_dir = Path(a.prov)
    steps = {}
    for p in sorted(prov_dir.glob("*.json")):
        steps[p.stem] = json.loads(p.read_text())

    # Cross-check: every recorded input digest against the file as it is now.
    stale = []
    for name, rec in steps.items():
        for path, digest in rec.get("inputs_sha256", {}).items():
            f = Path(path)
            if f.exists() and core.sha256(f) != digest:
                stale.append({"step": name, "input": path,
                              "recorded": digest[:12], "current": core.sha256(f)[:12]})

    git = "unknown"
    try:
        git = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                             text=True, cwd=Path(a.prov).parent.parent
                             ).stdout.strip() or "unknown"
    except Exception:
        pass

    manifest = {
        "run_id": a.run_id or time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_commit": git,
        "config_digest": cfg.digest(),
        "config_source": str(cfg.source) if cfg.source else None,
        "family": {"name": cfg["family.name"],
                   "species": cfg["family.species.name"],
                   "focal_genes": cfg.focal_genes,
                   "focal_digest": cfg.section_digest("family.focal_genes")},
        "environment": {"python": platform.python_version(),
                        "platform": platform.platform(),
                        **_pkg_versions("numpy", "pandas", "scipy", "Bio")},
        "n_steps_recorded": len(steps),
        "steps": steps,
        "stale_inputs": stale,
    }
    Path(a.out).write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    if stale:
        print(f"[run_manifest] WARNING {len(stale)} recorded input digest(s) no "
              f"longer match the file on disk: {stale}", file=sys.stderr)
    print(f"[run_manifest] {len(steps)} steps -> {a.out}")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

STEPS = {
    "fetch-reference": (step_fetch_reference,
                        ["which", "out", "force"]),
    "select-family": (step_select_family,
                      ["proteome", "hmm", "workdir", "out_faa", "out_members",
                       "threads"]),
    "cutoff-check": (step_cutoff_check,
                     ["hmm", "proteome", "members", "workdir", "out", "threads"]),
    "align": (step_align, ["faa", "out", "threads"]),
    "tree": (step_tree, ["aln", "prefix", "out_treefile", "out_report", "threads"]),
    "identity": (step_identity, ["aln", "out_matrix", "out_pairs"]),
    "embed": (step_embed, ["faa", "out_matrix", "out_npz", "device"]),
    "classify-pairs": (step_classify_pairs,
                       ["members", "gff3", "identity_pairs", "cosine",
                        "out_pairs", "out_context"]),
    "annotate-tree": (step_annotate_tree, ["treefile", "context", "out"]),
    "uniprot-map": (step_uniprot_map, ["members", "out"]),
    "fetch-afdb": (step_fetch_afdb, ["acc", "out", "out_meta", "force"]),
    "pocket-template": (step_pocket_template, ["out"]),
    "structures": (step_structures,
                   ["uniprot", "faa", "afdb_meta", "out_pdb_dir", "out"]),
    "pocket-residues": (step_pocket_residues,
                        ["pocket", "aln", "faa", "uniprot", "out", "out_cols"]),
    "structure-pairs": (step_structure_pairs,
                        ["aln", "pocket_cols", "structures", "pdb_dir", "out"]),
    "fetch-expression": (step_fetch_expression, ["gsm", "datadir", "force"]),
    "pseudobulk": (step_pseudobulk, ["datadir", "out", "force"]),
    "expression": (step_expression,
                   ["datadir", "outdir", "members", "counts_cache"]),
    "integration": (step_integration, ["results_dir", "outdir", "alpha"]),
    "run-manifest": (step_run_manifest, ["out", "run_id"]),
}

_LIST_ARGS = {"afdb_meta"}
_FLAG_ARGS = {"force"}
_INT_ARGS = {"threads"}
_FLOAT_ARGS = {"alpha"}
#: Optional string arguments. ``run_id`` defaults to a UTC timestamp so a run
#: is always identifiable even when nothing supplies one; ``device`` defaults to
#: auto-detect.
_OPTIONAL_ARGS = {"run_id", "device"}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    for name, (_fn, argnames) in STEPS.items():
        p = sub.add_parser(name, help=(_fn.__doc__ or "").strip().split("\n")[0])
        p.add_argument("--config", default=None,
                       help="path to config.yaml (default: config/config.yaml)")
        p.add_argument("--prov", default="results/.prov",
                       help="directory for provenance sidecars")
        for arg in argnames:
            flag = "--" + arg.replace("_", "-")
            if arg in _FLAG_ARGS:
                p.add_argument(flag, action="store_true")
            elif arg in _LIST_ARGS:
                p.add_argument(flag, nargs="+", required=True)
            elif arg in _INT_ARGS:
                p.add_argument(flag, type=int, default=4)
            elif arg in _FLOAT_ARGS:
                p.add_argument(flag, type=float, default=None)
            elif arg in _OPTIONAL_ARGS:
                p.add_argument(flag, default=None)
            else:
                p.add_argument(flag, required=True)
    return ap


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    cfg = cfgmod.load(a.config)
    fn, _ = STEPS[a.step]
    fn(a, cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
