#!/usr/bin/env python
"""Superseded by the workflow.

This runner decided its own step order and had no way to tell whether a
downstream output was older than the input it came from. That order now lives
in ``workflow/Snakefile`` as file dependencies, and each step it used to call is
a rule invoking ``pipeline/steps.py``.

    snakemake -s workflow/Snakefile --cores 8 --use-conda results/integration_module/redundancy_scores.csv

See docs/PIPELINE.md. Kept as a pointer rather than deleted so the commands in
older notes resolve to an explanation instead of a traceback.
"""

import sys

TARGET = "results/integration_module/redundancy_scores.csv"

print(__doc__.replace("{target}", TARGET), file=sys.stderr)
raise SystemExit(
    "This entry point no longer runs anything. Use:\n"
    f"    snakemake -s workflow/Snakefile --cores 8 --use-conda results/integration_module/redundancy_scores.csv\n"
)
