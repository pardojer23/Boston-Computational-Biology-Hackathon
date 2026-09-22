# Handoff — plan for the next run

Rewritten after module 4 landed. Everything the previous version of this file
planned has been built; what follows is what is left, and what you need to know
before touching it. The README is the reference document — this file only says
what to do next and why.

## 0. State in one paragraph

All four modules run end to end: sequence (`run_local.py`), structure
(`run_structure.py`), expression (`soy_globin_expression.py`) and integration
(`run_integration.py`, figure in `plot_integration.py`). Total runtime under a
minute with references and the count cache warm. The four pair tables join on a
canonical key that is now defined once in `soy_globin_core`, and module 4
reduces them to `results/integration_module/redundancy_scores.csv` — 21 pairs,
every raw and normalised component, symmetric and directional scores under two
normalisation scopes. The scientific limitation is single and specific: **there
is no cell-type resolution anywhere**, so the expression term is built from two
tissue means, and since every molecular measure is saturated within the Lb
clade, the whole ranking rests on that one quantity.

## 1. The one thing worth doing next

**Cell-type resolution on GSE226149.** README §5.1 has the detail; the short
version is that it converts `dose_ratio` — one number per pair from two tissue
means — into a per-cell distribution, and it is the only route to
`coexpression_overlap` and to the infected-vs-uninfected contrast the whole
design was built around. Budget 2–4 h.

Steps, in order:

1. Cluster the GSE226149 barcodes. The matrices are already in
   `data/expression/`, and `work/expression/pseudobulk_counts.csv.gz` is the
   summed table, not the per-cell one — you will need to re-read the matrices
   for this.
2. Annotate clusters against `results/expression_module/marker_profiles.csv`.
   **Known gap:** the uninfected-interstitial marker `Glyma.10G121524`
   (uricase-2 / nodulin-35) is absent from the GSE226149 feature list. Pick a
   substitute for that population before you start, not after.
3. Emit per-cell-type CPM per gene, then recompute `E` from cell-level
   co-expression overlap and cell-level dose, keeping the existing
   tissue-resolution columns beside the new ones so the two can be compared.
4. Re-run module 4. It will pick the new columns up through the same join; do
   not hand-edit `redundancy_scores.csv`.

This is also the only remaining step where **Modal is genuinely justified** —
holding the full cell × gene matrix for clustering is memory-bound, not
FLOP-bound. Everything else in this pipeline runs locally in seconds; see
README §7 for the stage-by-stage revision of the original Modal plan.

## 2. Decisions that are yours, not the code's

1. **Is `Glyma.10G198900` an ingroup member?** Module 2 settled the identity
   question: it is a named leghemoglobin, **GmLb5** (`LGB5_SOYBN`, A0A0R0HW51).
   It sits inside the chr10 array between Lbc3 and Lbc1 and is 77–80% identical
   to the focal four. Currently it is treated as an outgroup, so it appears in
   the 12 `Lb-other` pairs and not in the six scored `Lb-Lb` pairs. Promoting it
   to a fifth focal gene is a one-line change to `core.FOCAL_GENES` plus a
   re-run of all four modules — but note it is the weakest structural row in the
   family (pLDDT 81.6, and the only AFDB model whose sequence differs in length
   from the a4 primary transcript, +17 aa), so its pocket and TM values carry
   more uncertainty than the rest.
2. **Synteny confirmation for Lbc2.** Lbc2 is on Gm20 yet is the *most* similar
   member to Lbc3 on every molecular measure, which is why duplication mode is
   carried as an annotation and deliberately excluded from the score (README
   §4.5). MCScanX would establish Gm10/Gm20 homoeology properly, ~45 min, and is
   not budgeted. Asserting WGD origin from the literature may be enough for the
   hackathon — your call, but do not label anything `wgd` in an output file
   until one of the two has happened.
3. **α.** The committed score uses α = 0.40 on the molecular term. The top two
   focal pairs are 0.011 apart and swap at α ≈ 0.52, so the identity of the
   single most redundant pair is unresolved. If you want it resolved, that is an
   argument for better expression data, not for a better weight.

## 3. What not to re-litigate

These were decided with evidence and the evidence is in the repo. Reopen them
only with new data, not new preference.

- **GSE270392 is unusable for pairwise expression statistics on this family.**
  Not a judgement — a measurement. Lba ranks 45,386 of 52,594 there; the
  intended negative control Hb2 is the top globin, including inside the
  annotated infected-cell cluster; the atlas's own per-cell-state table filters
  out Lba, Lbc1 and GmLb5 entirely; and every focal-pair co-detection overlap
  sits within 1.2–2.3× of the independence expectation while the Hb2 control
  scores 2.6× higher than the best real pair. README §3.1 has the numbers. The
  dataset is still fine for qualitatively confirming that Lb transcript
  localizes to infected cells via Lbc2/Lbc3. One qualifier on that rejection:
  all seven family IDs matched the GSE270392 matrix rownames as **strings**,
  which is not the same as confirming they denote the same loci. Wm82 a2/a3/a4
  IDs share the `Glyma.10G199100` format without always referring to the same
  gene, and the proper check is to map through the GFF3 `ancestorIdentifier`
  attribute. It was not done, because the rejection rests on abundance rather
  than on identity — but if you ever use GSE270392 quantitatively, do that
  mapping first.
- **The score is multiplicative, not additive.** Two genes that never share a
  cell cannot buffer each other's loss.
- **Spearman over the five libraries is not a usable co-expression statistic
  here.** It is exactly 1.000 for all six focal pairs. It is emitted and unused.
- **The family is declared in one place.** `core.GENE_SYMBOLS` and
  `core.FOCAL_GENES`; `core.label_of` builds every label;
  `core.canonical_pair_order` enumerates every pair. It used to be declared
  twice and the two copies silently cost 6 of 21 rows on any join. Do not
  reintroduce a local family dict, and do not build pairs with a bare
  `itertools.combinations`.

## 4. Things that will waste your time if you don't know them

- Phytozome needs a JGI login; use the SoyBase URLs in
  `soy_globin_core.PROTEOME_URL` / `GFF3_URL`.
- GSE226149's `features.tsv.gz` uses `GLYMA_10G199100`, not `Glyma.10G199100`.
  Go through `soy_globin_expression.to_geo_id()` or your lookup silently matches
  nothing.
- Hugging Face weights do not come from `huggingface.co` — see README §7. The
  failure mode looks like a corrupt cache, not a blocked domain.
- Bioconda's `iqtree` is v3 and installs as `iqtree`, not `iqtree2`.
- `environment.yml` has a pip section now (`tmtools`), so the solve is
  two-stage. `matplotlib-base` is in there too; both were missing before and
  both are needed.
- The Claude Science → Modal dispatch path was broken during these runs
  (README §9.1), which is why the committed ESM2 matrix is fp32 on CPU. A
  separate user-managed `modal` env works for `modal run`.
- `soy_globin_core.py` imports nothing from Modal, on purpose. Keep it that
  way: it is what let these runs finish while Modal dispatch was down.
- Put new modules in their own `results/<name>_module/` directory and join on
  `label_a`/`label_b`. No module writes into another module's directory.
- `run_local.py --reuse-esm2` remaps a reused cosine matrix by gene ID, so
  adding a symbol to `core.GENE_SYMBOLS` does not invalidate it. If you change
  a *sequence*, that path is wrong — recompute.
