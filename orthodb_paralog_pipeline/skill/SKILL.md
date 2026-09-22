---
name: orthodb-paralog-pipeline
description: "Run an OrthoDB orthologous-group id through a full paralog-redundancy pipeline: fetch group sequences, build a MAFFT+IQ-TREE maximum-likelihood tree, and — assuming the group contains soybean (Glycine max) paralogs — map them onto the GSE270392 soybean expression atlas (tau, expression breadth, Benoit et al. 2025 dosage-balance/dominance/specialized/diverged classification, all-cell-type Pearson/Spearman/log2FC) and run a PRANK-aligned HyPhy aBSREL branch-site selection test on the soybean subtree. Use for questions like 'build a phylogeny and test paralog redundancy for OrthoDB group 706508at2759' or 'is this soybean gene family dosage-balanced or diverged'. Long-running stages (MAFFT/IQ-TREE, PRANK/HyPhy) dispatch to the user's Modal account."
---

# OrthoDB paralog-redundancy pipeline

Generalizes the leghemoglobin (OrthoDB group 706508at2759) analysis built in
this project into a reusable pipeline for **any** OrthoDB group id. The
phylogeny stage is species-agnostic; the expression and selection stages
assume **soybean** (Glycine max) — pass any group id, and the pipeline uses
whichever of its member sequences are soybean, skipping cleanly if none are.

Regression-tested against group 706508at2759 (leghemoglobin A): reproduces
197 sequences / 57 species, the Lbc3/Lbc2 "dosage balanced" Benoit
classification (coexpression 0.945, mean\|log2FC\| 0.76), and the aBSREL
result (0/5 branches significant, omega ranking Lbc2 0.12 < Lbc3 0.29 <
Lbc1 0.37 < Lba 0.99 < Node4 1.20).

## Files in this skill

- `kernel.py` — thin loader; auto-loads into your kernel and re-exports the
  functions below.
- `orthodb_fetch.py` — Stage 1 (species-agnostic).
- `soybean_expression.py` — Stage 3 (Glyma ID resolution + GSE270392 metrics).
- `soybean_selection.py` — Stage 4 prep (CDS resolution + tree pruning).
- `build_tree_modal.py`, `modal_absrel.py` — the two Modal-dispatching CLI
  scripts (Stages 2 and 4b — see the IMPORTANT note below).
- `report.py` — figures (tree colored by omega) + Markdown report.
- `pipeline.py` — a standalone CLI version of the whole flow (useful as a
  reference for the call sequence below, or to run outside a Claude Science
  session entirely).

## IMPORTANT: where each stage runs

Stages 1, 3, and the CDS/tree-prep half of stage 4 are pure Python (urllib +
pandas/numpy/scipy/Biopython) — call them directly via `kernel.py`'s
functions in a `python` cell.

**Stages 2 and 4b (the actual Modal dispatch) must run via the `bash` tool,
not a `python` cell.** Calling `subprocess` against a binary under a
host-granted path (the user's Modal-authenticated Python environment) is
blocked from inside the persistent `python`/`repl` kernels in this sandbox,
but works fine as a `bash` tool call. `kernel.py`'s `prepare_tree_build_command`
/ `prepare_selection_command` write the needed input files and hand back the
exact command string — run that string via `bash` (declaring the `Modal`
credential on the call), then load the result JSON back with
`load_json_result(path)`.

## Workflow

```python
# python cell
import os
os.makedirs("pipeline_out", exist_ok=True)

# Stage 1: fetch + clean (any OrthoDB group id)
fasta_path = "pipeline_out/group_clean.fasta"
group_summary = fetch_and_clean_group("<group_id>", fasta_path)
print(group_summary["n_genes"], group_summary["n_species"], group_summary["group_name"])
```

```bash
# bash tool cell -- declare credentials=["Modal"]; PATH/MODAL_PYTHON_BIN must
# point at the user's Modal-authenticated Python (ask_about_compute if unknown)
export MODAL_PYTHON_BIN=/path/to/modal/env/bin/python
$MODAL_PYTHON_BIN <skill_dir>/build_tree_modal.py pipeline_out/group_clean.fasta <group_id> pipeline_out/tree.json
```

```python
# python cell
tree_result = load_json_result("pipeline_out/tree.json")
full_tree_newick = tree_result["treefile"]

# Stage 3: soybean expression analysis (skips cleanly if no Glycine max in the group)
expression_result = run_soybean_expression_analysis(group_summary["records"], cache_dir="pipeline_out/_cache")
if not expression_result["skipped"]:
    for pc in expression_result["pair_classifications"]:
        if pc["classifiable"]:
            print(pc["gene_a"], pc["gene_b"], pc["benoit_group"], pc["rank_standardized_coexpression"])

# Stage 4a: resolve CDS + prune tree to the soybean paralogs with valid CDS
with_cds = resolve_cds_for_glyma_genes(expression_result["resolved_records"])
valid = [r for r in with_cds if r.get("cds_seq")]
target_labels = [tree_label(r["odb_id"], "Glycine max") for r in valid]
pruned_newick = prune_tree_to_labels(full_tree_newick, target_labels)

cds_fasta = "\n".join(f">{tree_label(r['odb_id'],'Glycine max')}\n{r['cds_seq']}" for r in valid)
with open("pipeline_out/cds.fasta", "w") as f: f.write(cds_fasta)
with open("pipeline_out/pruned_tree.nwk", "w") as f: f.write(pruned_newick)

print(prepare_selection_command("pipeline_out/cds.fasta", "pipeline_out/pruned_tree.nwk",
                                 "pipeline_out/absrel.json", modal_python="$MODAL_PYTHON_BIN"))
# -> run the printed command via the bash tool (Modal credential declared)
```

```python
# python cell, after the bash aBSREL run finishes
import json
modal_result = load_json_result("pipeline_out/absrel.json")
absrel_json = json.loads(modal_result["absrel_result"])
detected_labels = {tree_label(r["odb_id"], "Glycine max") for r in expression_result["resolved_records"]
                    if r.get("glyma_id") in expression_result["detected_glyma_ids"]}
plot_absrel_branches(absrel_json, detected_labels, "pipeline_out/absrel_figure.png")
write_markdown_report(group_summary, expression_result,
                       {"skipped": False}, "pipeline_out/report.md",
                       {"absrel_figure": "absrel_figure.png"})
```

## Known pitfalls (already fixed here, keep them fixed if you edit)

- **aBSREL unroots the tree.** Never assume the tested tree's internal-branch
  topology matches your rooted guide tree — `report.plot_absrel_branches`
  reads topology from `absrel_result["input"]["trees"]`, not from the guide
  tree you supplied. A prior version of this analysis drew the wrong branch
  as the tested internal edge because it assumed the guide tree's shape;
  see this project's `leghemoglobin_absrel_selection/README.md` for the
  incident.
- **Glyma ID resolution** is a two-tier chain: OrthoDB's own `GLYMA_##G######vN`
  tag when present, else BLASTP against SoyBase's public SequenceServer
  (`sequenceserver.soybase.org`) — submit via manual multipart-encoded
  `urllib` POST, not the `requests` library (it fails through this sandbox's
  proxy for this domain).
- **CDS/mRNA resolution** walks `Glyma ID -> GLYMA_##G######vN tag -> NCBI
  gene esearch -> GeneID -> elink gene_nuccore_refseqrna -> mRNA accession ->
  efetch CDS` — no genome-wide GFF download needed. A GeneID with no
  `gene_nuccore_refseqrna` link (e.g. an annotated pseudogene) is skipped, not
  an error.
- **NCBI E-utilities rate limit** at ~3 req/s without an API key; back off
  and retry on HTTP 429 rather than failing (see `soybean_selection._fetch`).
