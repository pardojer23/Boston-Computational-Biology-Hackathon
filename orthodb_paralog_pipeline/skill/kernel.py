"""
Thin loader for the orthodb-paralog-pipeline skill: locates this skill's own
directory on disk (via co_filename, since __file__ isn't set for sidecars),
adds it to sys.path, and re-exports the pipeline's functions so SKILL.md's
call sites (e.g. `fetch_and_clean_group(...)`) work directly.

The heavier module code (fetch, expression analysis, CDS/tree utilities,
report figures) lives in the sibling .py files in this skill directory, not
inline here -- see SKILL.md for the full workflow and which stages must run
via the `bash` tool instead of directly in this kernel (Modal dispatch).
"""
import sys


def skill_dir_path():
    here = sys._getframe().f_code.co_filename
    import os
    d = os.path.dirname(here)
    if not d:
        raise RuntimeError("skill directory unavailable in this runtime")
    return d


def ensure_skill_on_path():
    d = skill_dir_path()
    if d not in sys.path:
        sys.path.insert(0, d)


def fetch_and_clean_group(group_id: str, out_fasta_path: str):
    ensure_skill_on_path()
    import orthodb_fetch
    return orthodb_fetch.fetch_and_clean_group(group_id, out_fasta_path)


def run_soybean_expression_analysis(records, cache_dir: str):
    ensure_skill_on_path()
    import soybean_expression
    return soybean_expression.run_soybean_expression_analysis(records, cache_dir)


def resolve_cds_for_glyma_genes(resolved_records):
    ensure_skill_on_path()
    import soybean_selection
    return soybean_selection.resolve_cds_for_glyma_genes(resolved_records)


def tree_label(odb_id: str, organism: str):
    ensure_skill_on_path()
    import soybean_selection
    return soybean_selection.tree_label(odb_id, organism)


def prune_tree_to_labels(full_tree_newick: str, target_labels):
    ensure_skill_on_path()
    import soybean_selection
    return soybean_selection.prune_tree_to_labels(full_tree_newick, target_labels)


def prepare_tree_build_command(fasta_path: str, group_id: str, out_json_path: str, modal_python: str = None):
    """Returns the exact CLI command to run via the `bash` tool for the
    MAFFT+IQ-TREE Modal stage (this stage must NOT be called directly from
    this kernel -- see SKILL.md)."""
    ensure_skill_on_path()
    import os
    modal_python = modal_python or "python"
    script = os.path.join(skill_dir_path(), "build_tree_modal.py")
    return f"{modal_python} {script} {fasta_path} {group_id} {out_json_path}"


def prepare_selection_command(cds_fasta_path: str, tree_nwk_path: str, out_json_path: str, modal_python: str = None):
    """Returns the exact CLI command to run via the `bash` tool for the
    PRANK+HyPhy aBSREL Modal stage."""
    ensure_skill_on_path()
    import os
    modal_python = modal_python or "python"
    script = os.path.join(skill_dir_path(), "modal_absrel.py")
    return f"{modal_python} {script} {cds_fasta_path} {tree_nwk_path} {out_json_path}"


def load_json_result(path: str):
    import json
    with open(path) as f:
        return json.load(f)


def plot_absrel_branches(absrel_result_json, detected_labels, out_path: str):
    ensure_skill_on_path()
    import report
    return report.plot_absrel_branches(absrel_result_json, detected_labels, out_path)


def write_markdown_report(group_summary, expression_result, selection_result, out_path: str, figure_paths):
    ensure_skill_on_path()
    import report
    return report.write_markdown_report(group_summary, expression_result, selection_result, out_path, figure_paths)
