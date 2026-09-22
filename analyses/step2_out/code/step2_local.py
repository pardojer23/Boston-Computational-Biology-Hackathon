"""Step 2 of Li, Zhang & Schmitz 2025 (Plant Cell, koaf279), run locally.

Identical pipeline to dupgen_step2_modal.py — same perl expressions ported from
duplicated_gene_sets_age_server.sh, same family filtering from
07042025_tree_infor_glyma.Rmd. Moved off Modal because the whole step is 15 MB
of input and a tree walk: it is a light step under the agreed split.
"""

import csv
import glob
import io
import os
import re
import subprocess
import sys
import time

FAM = "https://data.soybase.org/LEGUMES/Fabaceae/genefamilies/legume.fam3.VLMQ"
HSH = "legume.fam3.VLMQ.sup1A_hsh.tsv.gz"
TREES = "legume.fam3.VLMQ.sup1B_trees.tar.gz"
WORK = "step2_work"
OUT = "step2_out"

log = []


def note(s):
    log.append(s)
    print(s, flush=True)


def sh(cmd, cwd=None):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        raise RuntimeError(f"failed: {cmd}\n{r.stderr[-1500:]}")
    return r.stdout


os.makedirs(WORK, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

if not os.path.exists(f"{WORK}/hsh.tsv"):
    sh(f"curl -sSL --max-time 600 {FAM}/{HSH} -o hsh.tsv.gz && gunzip -f hsh.tsv.gz", cwd=WORK)
if not glob.glob(f"{WORK}/**/Legume.fam3*", recursive=True):
    sh(f"curl -sSL --max-time 900 {FAM}/{TREES} -o trees.tar.gz && tar -xzf trees.tar.gz", cwd=WORK)

tree_files = sorted(f for f in glob.glob(f"{WORK}/**/Legume.fam3*", recursive=True)
                    if os.path.isfile(f))
note(f"tree files: {len(tree_files)}")
note(f"hsh.tsv lines: {sum(1 for _ in open(f'{WORK}/hsh.tsv'))}")

# ---- tree walk: verbatim port of duplicated_gene_sets_age_server.sh ----------
SIMPLIFY = r"""s/\(/</g; s/\)/>/g; s/([<,])(\w+)\.[^:]+:\d+\.\d+/$1$2/g; s/>\d+\.\d+/>/g; s/:\d+\.\d+//g"""
COUNT = r"""while (/<(glyma),(glyma)>/g) { print "<$1,$2>" }"""
EXTRACT = r"""
  while (/\((glyma[^:]+):[^,]+,(glyma[^:]+):[^\)]+\)/gi) {
    if (!$seen1++) { ($a1,$a2) = ($1,$2); }
    elsif (!$seen2++) { ($b1,$b2) = ($1,$2); }
  }
  END { if ($a1 && $a2 && $b1 && $b2) { print "$a1\t$a2\t$b1\t$b2\n"; } }
"""

sets_rows, n_two, n_other = [], 0, 0
t0 = time.time()
for i, fam in enumerate(tree_files):
    with open(fam) as fh:
        simp = subprocess.run(["perl", "-pe", SIMPLIFY], capture_output=True,
                              text=True, stdin=fh).stdout
    pairs = subprocess.run(["perl", "-nle", COUNT], input=simp,
                           capture_output=True, text=True).stdout.split()
    if len(pairs) != 2:
        n_other += 1
        continue
    n_two += 1
    with open(fam) as fh:
        out = subprocess.run(["perl", "-ne", EXTRACT], capture_output=True,
                             text=True, stdin=fh).stdout.strip()
    if out:
        sets_rows.append([os.path.basename(fam)] + out.split("\t"))
    if (i + 1) % 2000 == 0:
        el = time.time() - t0
        note(f"  {i+1}/{len(tree_files)} families | {el:.0f}s | eta {el/(i+1)*(len(tree_files)-i-1):.0f}s")

note(f"families with exactly two glyma clades: {n_two} (other: {n_other})")
note(f"glyma_output_sets rows: {len(sets_rows)}")

with open(f"{OUT}/glyma_output_sets.txt", "w") as fh:
    fh.write("fam\tgeneA1\tgeneA2\tgeneB1\tgeneB2\n")
    for r in sets_rows:
        fh.write("\t".join(r) + "\n")

# ---- family filtering: port of 07042025_tree_infor_glyma.Rmd ----------------
GID = re.compile(r"Glyma\.\d{2}G\d{6}")
fam_members = {}
with open(f"{WORK}/hsh.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 2:
            continue
        if "glyma.Wm82.gnm4.ann1." in p[1]:
            fam_members.setdefault(p[0], []).append(p[1])
note(f"families containing soybean genes: {len(fam_members)}")

four = {k: v for k, v in fam_members.items() if len(v) == 4}
note(f"families with exactly 4 soybean genes: {len(four)}")

wgd = set()
with open("wgd_unique_genes_lexicographic.csv") as fh:
    for r in csv.DictReader(fh):
        wgd.add(r["Dup"])
note(f"WGD-derived genes supplied from step 1: {len(wgd)}")

four_wgd = {}
for k, members in four.items():
    ids = [GID.search(m).group() for m in members if GID.search(m)]
    if len(ids) == 4 and all(g in wgd for g in ids):
        four_wgd[k] = ids
note(f"families whose 4 soybean genes are ALL WGD-derived: {len(four_wgd)}")

keep, missing_fam = [], 0
for r in sets_rows:
    fam_key = r[0]
    cand = [fam_key, fam_key.replace("Legume.", "legume."), os.path.splitext(fam_key)[0]]
    hit = next((c for c in cand if c in four_wgd), None)
    if hit is None:
        missing_fam += 1
        continue
    ext = [GID.search(g).group() if GID.search(g) else "" for g in r[1:5]]
    keep.append(r + ext)
note(f"final duplicated gene SETS: {len(keep)} (dropped {missing_fam} not in 4x-WGD families)")

with open(f"{OUT}/duplicated_gene_sets_updated.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["fam", "geneA1", "geneA2", "geneB1", "geneB2",
                "geneA1_extracted", "geneA2_extracted",
                "geneB1_extracted", "geneB2_extracted"])
    w.writerows(keep)

open(f"{OUT}/step2_run.log", "w").write("\n".join(log) + "\n")
print("WROTE", OUT)
