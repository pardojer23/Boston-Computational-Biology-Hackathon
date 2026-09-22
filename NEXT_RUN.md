# Handoff — plan for the next run

Read `README.md` first; it records what exists and why. This file is only
what to do next.

---

## 0. State in one paragraph

Sequence module is done and its outputs are committed under
`results/sequence_module/`. The family is **7 proteins** (4 focal Lbs, a 5th
chr10 cluster member `Glyma.10G198900`, and a chr11 phytoglobin pair that
roots the tree). Two results constrain everything downstream: within the Lb
clade both sequence identity (91.7–95.2%) and ESM2 cosine distance
(0.0006–0.0024) are **saturated** and discriminate nothing; and Lbc2, the
lone chr20 copy, is the *most* identical member to Lbc3 (95.2%), so genomic
adjacency does not predict sequence similarity here. Structure, expression
and integration modules do not exist yet.

---

## 1. Do these three checks before writing code

1. **Is Modal dispatch fixed?** Call `list_compute`; if `byoc:modal` is
   listed, try one trivial `submit_job` before planning around it. During the
   last run the provider was enabled but its helper env
   (`~/.claude-science/conda/envs/compute-provider-modal/`) had never been
   created, so every submit failed. If it is still broken, everything below
   still works — run the GPU stages via `modal run` from the user's own
   `modal` conda env instead, or on CPU where tolerable.
2. **Re-run the sequence module and confirm it reproduces.**
   `python pipeline/run_local.py --reuse-esm2` should finish in ~8 s and
   leave `results/sequence_module/` unchanged apart from the timestamp. If
   it doesn't, resolve that before layering anything on top.
3. **Ask the human the two open questions in §6** — both change the score,
   and neither is answerable from the data.

---

## 2. Module 2 — structure (~1.5 h)

**Goal:** a pairwise structural-distance matrix over the same 7 labels, so it
drops into `paralog_pairs.csv` beside `pid_aligned` and
`esm2_cosine_distance`.

1. Try **AlphaFold DB** first — soybean is in the proteome set, and a
   download beats a prediction. Look up by UniProt accession; map
   `Glyma.*` → UniProt via UniProt's ID-mapping API. Expect misses for the
   less-characterised members.
2. Predict whatever is missing with **ESMFold** on Modal (A10G is enough for
   161-residue proteins; this is genuinely Modal-worthy). Keep the weights in
   the existing `soy-globin-hf-cache` Volume.
3. Superpose all-vs-all and record **TM-score** (US-align or TM-align;
   `tmtools` if you want it in-process) plus backbone RMSD. TM-score is the
   one to carry forward — RMSD is length-sensitive and these proteins differ
   by up to 17 residues.
4. Add the **distal histidine and proximal F8 histidine** positions and the
   heme-pocket residues as an annotation column. For leghemoglobins the
   functionally interesting variation is in the pocket, not the global fold,
   and global TM-score will almost certainly be saturated (>0.95) across all
   four Lbs — expect this and plan to report pocket-residue identity as the
   structural discriminator instead.

**Prediction, stated in advance so it can be falsified:** global structural
similarity will be as uninformative within the Lb clade as sequence identity
was. If TM-score is >0.95 for all six Lb pairs, say so and down-weight it
rather than pretending it contributes.

**Output:** `results/structure_module/structure_pairs.csv` with
`label_a,label_b,tm_score,rmsd,pocket_identity` and the structures under
`results/structure_module/pdb/`.

---

## 3. Module 3 — expression (~2 h, the one that matters)

**Goal:** per-gene expression profiles across nodule cell types, and a
pairwise co-expression / overlap measure. This module carries the
discriminating signal; budget accordingly.

Source: **Zhang et al. 2024 Cell**, GEO **GSE270392** — nodule snRNA-seq with
annotated infected cells. Cross-check against soybean-atlas.com and SoyBase
where convenient.

1. Pull the processed matrix and cell-type annotations from GEO. Check what
   format they deposited before planning memory: if it is a full raw count
   matrix, subset to the 7 genes immediately — you do not need the other
   52,865, and a 7-gene subset makes the rest of the module trivial locally.
   **Only reach for Modal if the object must be loaded whole**; the
   constraint here is memory, not FLOPs.
2. Confirm the atlas uses Wm82.a4.v1 gene IDs. If it is on a different
   annotation version (a2/a3 IDs look like `Glyma.10G199100` too but do not
   always refer to the same locus), map through the `ancestorIdentifier`
   attribute already present in the GFF3 — `parse_gff_genes()` currently
   discards it, so extend that function rather than writing a new parser.
3. Compute per gene: mean expression and detection rate **per annotated cell
   type**, with infected vs uninfected nodule cells kept separate; and a
   tissue-specificity index (τ) across cell types.
4. Compute per pair: Spearman correlation of the cell-type profiles, plus a
   **co-expression overlap** — the fraction of cells expressing both relative
   to the fraction expressing either. These measure different things and the
   score below uses the second.

**Output:** `results/expression_module/gene_celltype_profiles.csv` and
`expression_pairs.csv` keyed on the same `label_a,label_b`.

---

## 4. Module 4 — integration and the redundancy score (~1 h)

### Use a gated, not an additive, score

The obvious form is a weighted sum of similarities. **It is wrong for this
question.** Two genes that are 95% identical but never expressed in the same
cell are not functionally redundant — they cannot buffer each other's loss.
An additive score gives that pair a high value on the strength of sequence
alone. Redundancy is a conjunction, so the score should be multiplicative:

```
R(i,j) = M(i,j) ** alpha  *  E(i,j) ** beta,     alpha + beta = 1
```

where

- `M` = **molecular interchangeability** ∈ [0,1] — can protein i do protein
  j's job? Combine the saturated features here, since they agree:
  `M = w1*norm(pid_aligned) + w2*(1 - norm(esm2_cosine)) + w3*tm_score`
  with `w1=w2=w3=1/3`. Normalise across the observed family range, not
  0–100, or every Lb pair collapses to ~1.0.
- `E` = **expression overlap** ∈ [0,1] — are they ever in the same place at
  the same time? Use the co-expression overlap from §3, not the correlation.

Recommended `alpha = 0.35`, `beta = 0.65`. The asymmetry is not arbitrary: it
follows from the measured fact that `M` has almost no variance within the Lb
clade (§1.5, README §2b) while `E` is the only term that can separate them.
State the weights in the output, and emit the unweighted components too so a
reader can reweight without rerunning.

### Duplication mode is a modifier, not a term

Do not put `duplication_mode` in the score. Use it to *interpret* the result:
tandem pairs that score high are the candidates for genuine buffering;
dispersed (likely homoeologous) pairs that score high suggest sub-
functionalisation has not happened since the WGD. Report it as a column.

### Validate against something

The score is worth nothing without a sanity check. The strongest cheap one:
**Lba is the dominant nodule leghemoglobin** and the Lbc genes are
co-expressed with it — a score that does not put the Lba/Lbc pairs at the top
of the ranking is telling you something is wrong. Published RNAi/CRISPR work
knocking out the soybean Lb cluster is the other external check; if you cite
it, fetch and verify the paper rather than recalling the phenotype.

**Output:** `results/redundancy_scores.csv` (all pairs, all components, the
final `R`, the weights used) and a single summary figure — a 7×7 heatmap of
`R` with the tree as a side dendrogram, or a scatter of `M` vs `E` with pairs
labelled and duplication mode as the marker. One figure, not a gallery.

---

## 5. Time budget

Assumes the sequence module is already done (it is).

| block | time |
|---|---|
| §1 checks + Modal sanity | 0.25 h |
| Module 2 — structure | 1.5 h |
| Module 3 — expression | 2.0 h |
| Module 4 — integration + figure | 1.0 h |
| Writing up: extend README, record what broke | 0.5 h |
| Slack | 0.75 h |
| **total** | **6.0 h** |

If time runs short, **cut the structure module, not the expression module.**
§2's prediction is that structure adds little within the Lb clade; the atlas
is where the answer is. A defensible `R` can be computed with `M` from
sequence alone.

---

## 6. Open questions for the human

1. **Is `Glyma.10G198900` a leghemoglobin?** It sits inside the chr10 tandem
   array between Lbc3 and Lbc1, is 77–80% identical to the four focal genes,
   and carries a full-length globin domain — but the user's focal list omits
   it. Is it a known pseudogene, a fifth Lb under another name, or simply not
   of interest? This changes whether it belongs in the ingroup or is a fifth
   test case.
2. **Synteny confirmation for Lbc2.** Should the next run do the MCScanX
   analysis to establish Gm10/Gm20 homoeology properly, or is asserting WGD
   origin from the literature acceptable for the hackathon? The first costs
   ~45 min and is not in the budget above.

## 7. Things that will waste your time if you don't know them

- Phytozome needs a JGI login; use the SoyBase URLs in
  `soy_globin_core.PROTEOME_URL` / `GFF3_URL`.
- Hugging Face weights do not come from `huggingface.co` — see README §4. The
  failure mode looks like a corrupt cache, not a blocked domain.
- Bioconda's `iqtree` is v3 and installs as `iqtree`, not `iqtree2`.
- `soy_globin_core.py` imports nothing from Modal, on purpose. Keep it that
  way: it is what let the last run finish while Modal dispatch was down.
- Do not add a stage that writes into `results/sequence_module/` from outside
  `run_local.py`; put new modules in their own directory and join on
  `label_a`/`label_b`.
