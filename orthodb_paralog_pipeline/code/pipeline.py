"""
Orchestrator: run the full OrthoDB-group-to-selection-test pipeline for any
OrthoDB group id.

Stages (species-agnostic where noted):
  1. Fetch the group + clean FASTA labels                    [any species]
  2. MAFFT align + IQ-TREE ML tree, on Modal                  [any species]
  3. Soybean expression analysis (GSE270392), if the group    [soybean-specific]
     contains Glycine max records
  4. Soybean branch-site selection (PRANK + HyPhy aBSREL,      [soybean-specific]
     on Modal), on the soybean paralogs with resolvable CDS
  5. Figures + a Markdown report

Stages 3-4 assume soybean per this pipeline's scope (see README); they skip
cleanly (with a stated reason) if the group has no Glycine max records, or
if soybean records exist but none resolve to usable data.

Usage:
    python pipeline.py <orthodb_group_id> <output_dir>
"""
import json
import os
import sys

import orthodb_fetch
import report
import soybean_expression
import soybean_selection


def run_orthodb_pipeline(group_id: str, outdir: str, cache_dir: str = None,
                          run_tree_on_modal: bool = True, full_tree_newick_override: str = None) -> dict:
    os.makedirs(outdir, exist_ok=True)
    cache_dir = cache_dir or os.path.join(outdir, "_cache")

    # --- Stage 1: fetch ---
    fasta_path = os.path.join(outdir, f"{group_id}_clean.fasta")
    group_summary = orthodb_fetch.fetch_and_clean_group(group_id, fasta_path)
    print(f"[1/5] Fetched {group_summary['n_genes']} genes / {group_summary['n_species']} species "
          f"for {group_id} ({group_summary['group_name']!r})")

    # --- Stage 2: tree ---
    tree_result = None
    full_tree_newick = None
    if run_tree_on_modal:
        import build_tree_modal
        tree_result = build_tree_modal.build_tree_via_subprocess(fasta_path, group_id)
        full_tree_newick = tree_result.get("treefile")
        tree_out_path = os.path.join(outdir, f"{group_id}_tree.json")
        with open(tree_out_path, "w") as f:
            json.dump(tree_result, f)
        print(f"[2/5] Tree built: IQ-TREE returncode={tree_result.get('iqtree_returncode')}")
    elif full_tree_newick_override:
        full_tree_newick = full_tree_newick_override
        print("[2/5] Using externally supplied tree (run_tree_on_modal=False)")
    else:
        print("[2/5] Skipped tree build (run_tree_on_modal=False, no override supplied)")

    # --- Stage 3: soybean expression ---
    expression_result = soybean_expression.run_soybean_expression_analysis(
        group_summary["records"], cache_dir=cache_dir)
    if expression_result.get("skipped"):
        print(f"[3/5] Soybean expression analysis skipped: {expression_result['reason']}")
    else:
        print(f"[3/5] Soybean expression analysis: {expression_result['n_gmax_records']} Gmax record(s), "
              f"{len(expression_result['detected_glyma_ids'])} detected in GSE270392")

    # --- Stage 4: soybean selection ---
    selection_result = {"skipped": True, "reason": "no ML tree available"}
    if full_tree_newick and not expression_result.get("skipped"):
        selection_result = soybean_selection.run_soybean_selection_analysis(
            expression_result["resolved_records"], full_tree_newick)
        if selection_result.get("skipped"):
            print(f"[4/5] Selection analysis skipped: {selection_result['reason']}")
        else:
            print(f"[4/5] Selection analysis: {selection_result['n_taxa']} taxa, "
                  f"hyphy_returncode={selection_result['modal_result'].get('hyphy_returncode')}")
    else:
        print(f"[4/5] Selection analysis skipped: {selection_result['reason']}")

    # --- Stage 5: figures + report ---
    figure_paths = {}
    if not selection_result.get("skipped"):
        absrel_json = json.loads(selection_result["modal_result"]["absrel_result"])
        detected_labels = {soybean_selection.tree_label(r["odb_id"], "Glycine max")
                            for r in expression_result["resolved_records"]
                            if r.get("glyma_id") in expression_result["detected_glyma_ids"]}
        absrel_fig_path = os.path.join(outdir, f"{group_id}_absrel.png")
        report.plot_absrel_branches(absrel_json, detected_labels, absrel_fig_path)
        figure_paths["absrel_figure"] = os.path.basename(absrel_fig_path)

    report_path = os.path.join(outdir, f"{group_id}_report.md")
    report.write_markdown_report(group_summary, expression_result, selection_result, report_path, figure_paths)
    print(f"[5/5] Wrote report to {report_path}")

    return {
        "group_summary": {k: v for k, v in group_summary.items() if k != "records"},
        "tree_result_present": tree_result is not None,
        "expression_result_skipped": expression_result.get("skipped"),
        "selection_result_skipped": selection_result.get("skipped"),
        "report_path": report_path,
        "figure_paths": figure_paths,
    }


def main():
    group_id, outdir = sys.argv[1], sys.argv[2]
    run_tree_on_modal = "--no-tree" not in sys.argv
    summary = run_orthodb_pipeline(group_id, outdir, run_tree_on_modal=run_tree_on_modal)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
