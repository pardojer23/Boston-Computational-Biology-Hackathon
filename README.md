# Soybean leghemoglobin redundancy — Boston Computational Biology Hackathon

Predicting ortholog/paralog **functional redundancy** in soybean (*Glycine max*),
using the leghemoglobin family as the test case. The design is four modules —
sequence, structure, expression, integration — feeding one weighted redundancy
score.

## Status

| module | state | outputs |
|---|---|---|
| 1 — sequence | **done** | `results/sequence_module/` |
| 2 — structure | **done** | `results/structure_module/` |
| 3 — expression | **done, at pseudobulk resolution only** | `results/expression_module/` |
| 4 — integration / redundancy score | **not started** | — |

Three pair tables now exist on a common key (`label_a`, `label_b`) and are
ready to join: `paralog_pairs.csv`, `structure_pairs.csv`,
`expression_pairs.csv`. §5 lists exactly what is still missing, including the
two expression statistics that were requested and are not yet computed.

`NEXT_RUN.md` is the earlier handoff. Its §2 (structure) and §3 (expression)
have since been executed and are superseded by §2/§3 below; its §4 (the score
design) and §6 (open questions for the human) still stand.

Focal genes (Wm82.a4.v1 IDs):

| symbol | gene ID | chromosome |
|---|---|---|
| Lba | `Glyma.10G199100` | Gm10 |
| Lbc1 | `Glyma.10G199000` | Gm10 |
| Lbc2 | `Glyma.20G191200` | Gm20 |
| Lbc3 | `Glyma.10G198800` | Gm10 |

---

## 1. Module 1 — sequence

Runner: `pipeline/run_local.py`. Wall time 11.3 s with references cached.

### 1.1 Reference data

Phytozome requires a JGI login, so the Wm82.a4.v1 files were taken from
**SoyBase**, which serves the identical annotation under its own release name
`Wm82.gnm4.ann1.T8TQ`:

- `glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz` — 52,872 primary-transcript proteins
- `glyma.Wm82.gnm4.ann1.T8TQ.gene_models_exons.gff3.gz` — gene models

The **PF00042 (Globin)** HMM came from the InterPro API
(`/entry/pfam/PF00042?annotation=hmm`), resolving to **PF00042.29**, model
length 117, gathering threshold **GA = 22**.

sha256 of all three inputs is recorded in
`results/sequence_module/manifest.json`. The downloads are *not* committed —
`pipeline/run_local.py` re-fetches them into `data/`.

> Identifier convention: proteins are `glyma.Wm82.gnm4.ann1.Glyma.10G199100.1`,
> GFF3 seqids are `glyma.Wm82.gnm4.Gm10`. `soy_globin_core.gene_of()` strips
> both down to the bare `Glyma.10G199100` used everywhere else.

### 1.2 Family selection — 7 proteins

`hmmsearch --cut_ga` against the full primary proteome, unioned with the four
focal gene IDs (belt and braces — none of them actually needed rescuing).

| gene | symbol | length | bit score | E-value | PF00042 coverage |
|---|---|---|---|---|---|
| `Glyma.11G121800` | — | 161 | 79.6 | 7.1e-22 | 0.99 |
| `Glyma.11G121700` | — | 157 | 76.7 | 5.7e-21 | 0.99 |
| `Glyma.20G191200` | **Lbc2** | 145 | 69.2 | 1.2e-18 | 0.98 |
| `Glyma.10G198900` | — | 151 | 67.4 | 4.3e-18 | 0.98 |
| `Glyma.10G199000` | **Lbc1** | 144 | 66.7 | 7.1e-18 | 0.98 |
| `Glyma.10G199100` | **Lba** | 144 | 66.0 | 1.2e-17 | 0.98 |
| `Glyma.10G198800` | **Lbc3** | 145 | 63.2 | 8.9e-17 | 0.98 |

Seven felt low for a paleopolyploid, so a relaxed search (`-E 10`) was run to
check the cutoff wasn't truncating the family. It isn't: the 8th-best protein
in the genome scores **13.5** (E = 0.23), a ~50-bit gap below the weakest real
hit. The three non-focal members are a fifth chr10 cluster gene
(`Glyma.10G198900`) and a chr11 tandem pair (`Glyma.11G121700/121800`) that
are non-symbiotic phytoglobins — these are what root the Lb clade.

Module 2 later resolved all three against UniProt: `Glyma.10G198900` is
**GmLb5 / leghemoglobin 5** (`LGB5_SOYBN`, A0A0R0HW51), and the chr11 pair are
the anaerobic nitrite reductases **Hb1** (`NSHB1_SOYBN`) and **Hb2**
(`NSHB2_SOYBN`). That answers half of open question 1 in `NEXT_RUN.md` §6 —
`Glyma.10G198900` is a named leghemoglobin, not an unannotated locus — but not
whether it belongs in the ingroup, which is still the human's call.

**This set does not include class-3 truncated hemoglobins**, which match
**PF01152** (Bacterial-like globin), not PF00042. If you want a deeper
outgroup, that is the one-line change.

### 1.3 Alignment and phylogeny

MAFFT **L-INS-i** (`--localpair --maxiterate 1000`) → 162 columns. IQ-TREE
with ModelFinder, 1000 ultrafast bootstrap replicates, 1000 SH-aLRT
replicates, fixed seed 20240601. Best-fit model **Q.PLANT+I**.

```
(Glyma.20G191200_Lbc2:0.0153934781,((Glyma.10G199000_Lbc1:0.0257420508,
((Glyma.11G121800:0.1010816599,Glyma.11G121700:0.1164176392)100/100:1.6993936020,
Glyma.10G198900:0.0888325425)79.7/70:0.1463466693)45.4/47:0.0154733173,
Glyma.10G199100_Lba:0.0478835571)31.1/39:0.0126731124,Glyma.10G198800_Lbc3:0.0419615774);
```

Node labels are `SH-aLRT/UFBoot`. The phytoglobin split is 100/100 on a very
long branch (1.70), as expected for an outgroup. **Within the Lb clade the
support is poor** — 45.4/47 and 31.1/39 — which is not a bug: four sequences
at >91% identity over 162 columns do not contain enough signal to resolve
their branching order. Do not build an argument on the internal Lb topology.

### 1.4 Pairwise identity

Computed from the MSA (not from independent pairwise alignments), so it is
consistent with the tree. Percent identity over columns where **both**
sequences have a residue:

|  | Lbc2 | Lbc1 | Lba | Lbc3 | 11G121800 | 11G121700 | 10G198900 |
|---|---|---|---|---|---|---|---|
| **Lbc2** (Gm20) | — | 93.8 | 93.1 | **95.2** | 43.4 | 45.5 | 80.0 |
| **Lbc1** (Gm10) | 93.8 | — | 92.4 | 91.7 | 43.8 | 45.8 | 79.9 |
| **Lba** (Gm10) | 93.1 | 92.4 | — | 92.4 | 43.1 | 45.1 | 77.8 |
| **Lbc3** (Gm10) | **95.2** | 91.7 | 92.4 | — | 42.8 | 44.1 | 77.2 |

A second convention (`pid_shorter`, identities ÷ shorter ungapped length) is
in `pairwise_identity_pairs.csv`; it matters only for partial sequences, which
this family doesn't have.

### 1.5 ESM2-650M embeddings

`facebook/esm2_t33_650M_UR50D`, final hidden layer, **mean-pooled over residue
tokens only** — BOS, EOS and PAD are masked out before averaging, which
matters because sequences of different length are batched together. Cosine
distance = 1 − cosine similarity of the 1280-d pooled vectors.

| pair class | cosine distance |
|---|---|
| within the four Lbs | 0.0006 – 0.0024 |
| any Lb vs `Glyma.10G198900` | 0.0050 – 0.0058 |
| any Lb vs `Glyma.11G121700` | 0.0731 – 0.0758 |
| any Lb vs `Glyma.11G121800` | 0.1426 – 0.1447 |

**This stage ran on CPU, not on a GPU** (see §9.1). The Modal path uses fp16 on
an A10G and will differ in the fourth decimal; the committed matrix is fp32.

### 1.6 Duplication mode

From GFF3 gene order alone — every gene on a sequence is ranked by start
coordinate, and a pair is called by how many protein-coding genes sit between
its members:

- `tandem` — same sequence, ≤ 10 intervening genes
- `proximal` — same sequence, more intervening genes but ≤ 1 Mb apart
- `dispersed` — different sequence, or beyond the proximal window

Result: **7 tandem, 14 dispersed, 0 proximal** across 21 pairs. The chr10
cluster is four contiguous genes, all on the minus strand:

| gene | symbol | start | end | rank on Gm10 |
|---|---|---|---|---|
| `Glyma.10G198800` | Lbc3 | 43,082,011 | 43,083,137 | 1836 |
| `Glyma.10G198900` | GmLb5 | 43,085,930 | 43,089,182 | 1837 |
| `Glyma.10G199000` | Lbc1 | 43,090,608 | 43,091,995 | 1838 |
| `Glyma.10G199100` | Lba | 43,095,323 | 43,097,067 | 1839 |

`Glyma.20G191200` (Lbc2) is the single Gm20 copy at 42,955,083–42,956,301.

---

## 2. Module 2 — structure

Runner: `pipeline/run_structure.py`. Wall time 12.8 s. **No GPU was used and
no structure was predicted** — `NEXT_RUN.md` budgeted ESMFold on Modal for
this stage and it turned out to be unnecessary.

### 2.1 Structures come from AlphaFold DB

Every member of the family has a **reviewed** UniProt entry, so every member
has an AFDB model. Genes were resolved to accessions through the UniProt
search API (taxon 3847), then models pulled through
`alphafold.ebi.ac.uk/api/prediction/<acc>`.

| label | UniProt | entry | protein name | mean pLDDT |
|---|---|---|---|---|
| Lba | `P02238` | LGB3_SOYBN | Leghemoglobin 3 / GmLba | 96.5 |
| Lbc1 | `P02235` | LGB2_SOYBN | Leghemoglobin 2 / GmLbc1 | 95.8 |
| Lbc2 | `P02236` | LGB4_SOYBN | Leghemoglobin 4 / GmLbc2 | 95.7 |
| Lbc3 | `P02237` | LGB1_SOYBN | Leghemoglobin 1 / GmLbc3 | 96.1 |
| `Glyma.10G198900` | `A0A0R0HW51` | LGB5_SOYBN | Leghemoglobin 5 / GmLb5 | 81.6 |
| `Glyma.11G121700` | `I1LJI1` | NSHB1_SOYBN | Anaerobic nitrite reductase Hb1 | 94.4 |
| `Glyma.11G121800` | `Q42785` | NSHB2_SOYBN | Anaerobic nitrite reductase Hb2 | 93.8 |

All seven are AFDB model_v6 (created 2025-08-01). Note the older AFDB *file*
URL pattern (`AF-<acc>-F1-model_v4.pdb`) no longer resolves — go through the
API and use the URL it hands back.

**AFDB models are keyed on the UniProt sequence, which is not always the
Wm82.a4 primary transcript.** Each model is therefore aligned back to the a4
protein and the identity and length delta recorded in `structures.csv`. Six of
seven are byte-identical to the a4 sequence. The exception is
**`Glyma.10G198900` / GmLb5**: the UniProt sequence is 168 aa against 151 aa
for the a4 primary transcript — 100% identical over the 151 aligned columns,
i.e. a 17-residue extension rather than a different gene. Its pLDDT is also
the lowest in the set (81.6, with 7.7% of residues very low). Treat GmLb5
structural numbers as the least reliable row in the table.

### 2.2 The heme pocket is transferred from a crystal structure

AFDB models are apo, so there is no ligand to measure contacts against.
Residues with any heavy atom within **5.0 Å** of any heme heavy atom in
**1BIN** (soybean leghemoglobin-a, 2.20 Å) define the pocket: **23 residues**.
Those positions are transferred onto every family member through the MAFFT
alignment module 1 already produced — all 23 mapped, none dropped.

As an independent check, UniProt's own heme-binding annotations for P02238
were fetched separately: **4 of 5 annotated sites fall inside the
crystal-derived pocket**. The annotations are a check on the transfer, not its
source.

### 2.3 Global fold is saturated; the pocket is not

`NEXT_RUN.md` §2 predicted, in advance, that global TM-score would be
uninformative within the Lb clade (>0.95 for all six pairs). **That prediction
held.** TM-score is normalised by the shorter chain — the conservative choice,
since it cannot be inflated by one protein being a fragment of the other.

| pair | TM-score | RMSD (Å) | pocket identity (23 residues) |
|---|---|---|---|
| Lbc2 – Lbc3 | **0.9976** | 0.221 | **100.0** |
| Lbc2 – Lbc1 | 0.9915 | 0.428 | 95.7 |
| Lbc1 – Lbc3 | 0.9904 | 0.460 | 95.7 |
| Lbc2 – Lba | 0.9863 | 0.919 | 91.3 |
| Lba – Lbc3 | 0.9844 | 0.915 | 91.3 |
| Lbc1 – Lba | 0.9816 | 0.992 | 95.7 |
| any Lb – GmLb5 | 0.939 – 0.952 | 1.05 – 1.61 | 82.6 – 87.0 |
| any Lb – chr11 Hb1/Hb2 | 0.820 – 0.918 | 1.42 – 2.41 | 52.2 – 56.5 |

Two things to carry forward:

1. **Global TM-score adds nothing within the clade** (0.9816–0.9976, a range of
   0.016). Fold it into the molecular-interchangeability term if you like, but
   do not expect it to move the ranking.
2. **Pocket identity does vary** — 91.3% to 100% — and it varies *in the same
   direction as sequence identity*: Lbc2–Lbc3 is the top pair on all three
   measures (95.2% sequence identity, 0.9976 TM, identical pocket), and both
   pairs involving Lba are the bottom (91.3%). Across 23 positions that is 2
   substitutions versus 0, so it is a weak discriminator — but it is the only
   structural quantity here that discriminates at all, and
   `pocket_residues.csv` shows exactly which positions differ.

---

## 3. Module 3 — expression (pseudobulk)

Runner: `pipeline/soy_globin_expression.py` (module and CLI in one file).

### 3.1 The atlas named in the plan was rejected, with a diagnostic

`NEXT_RUN.md` §3 specified **GSE270392** (Zhang et al. 2024 *Cell*) — snRNA-seq
of *early* nodule with annotated infected cells. It is not usable for this
question, and the reason is measurable rather than a matter of taste: nuclear
RNA is depleted of abundant stable cytoplasmic transcripts, and leghemoglobin
is the extreme case of such a transcript.

In GSE270392's 29,436 nuclei:

| gene | total counts | rank (of 52,594) |
|---|---|---|
| Lba | 668 | 45,386 |
| Lbc1 | 777 | 44,218 |
| Lbc2 | 1,096 | 40,807 |
| Lbc3 | 3,038 | 21,481 |
| GmLb5 | 865 | 43,245 |
| Hb1 | 1,733 | 34,031 |
| **Hb2** (intended negative control) | **7,663** | **4,343** |

Median gene total in that matrix is 2,496, so three of the four Lbs sit below
the median gene, and the most abundant globin in the sample is the chr11
non-symbiotic haemoglobin that was supposed to be the outgroup. Co-expression
there would rest on 9–33 co-detected nuclei per Lb pair.

**GSE226149 is used instead** — mature nodule + root, protoplast scRNA-seq —
as **pseudobulk**: counts summed over every barcode in a library. Five
libraries, 14.7–48.6 M counts each: 2 nodule (`GSM7065810/11`), 3 root
(`GSM7065807/08/09`). In it, Lba is **rank 1 of 57,147 genes** and alone
carries 1.86% of all counts in the sample.

![Dataset diagnostic](results/expression_module/dataset_diagnostic.png)

Panel **b** is the reason to trust the switch: transcript share of total
leghemoglobin computed here (Lba 51.5%, Lbc1 21.7%, Lbc3 14.4%, Lbc2 12.4%)
recovers the published **protein** rank order of the four isoforms. The
published protein values plotted alongside come from the literature and are
**not sourced in-repo** — see §9.7.

### 3.2 What pseudobulk can and cannot give

Cell calling, clustering and cell-type annotation deliberately do **not**
happen in this module: the published cell-type labels are not in the GEO
submission, and re-deriving them would substitute a different piece of work
for the one the design asks for. The consequence is stated rather than worked
around — `detection_rate` and `coexpression_overlap` are per-*cell* quantities
and pseudobulk has no cells, so both are emitted as **NA with a reason
string**, not replaced by a look-alike statistic.

Per-gene profiles (CPM over summed library counts):

| gene | symbol | nodule mean CPM | root mean CPM | τ |
|---|---|---|---|---|
| `Glyma.10G199100` | Lba | 19,452 | 0.06 | >0.9999 |
| `Glyma.10G199000` | Lbc1 | 8,577 | 0.02 | >0.9999 |
| `Glyma.10G198800` | Lbc3 | 6,774 | 0.01 | >0.9999 |
| `Glyma.20G191200` | Lbc2 | 5,055 | 0.01 | >0.9999 |
| `Glyma.10G198900` | GmLb5 | 1.95 | 0.00 | 1.000 |
| `Glyma.11G121800` | Hb2 | 28.1 | 10.4 | 0.630 |
| `Glyma.11G121700` | Hb1 | 0.45 | 0.66 | 0.319 |

The four Lbs are nodule-exclusive to five orders of magnitude; the chr11
phytoglobins are not tissue-specific at all. **τ is computed over two tissue
means, so it collapses to a nodule-vs-root contrast** — do not read it as a
tissue-specificity index in the usual sense.

`marker_profiles.csv` carries four nodule markers as controls. NOD26
(`Glyma.08G120100`, symbiosome membrane, infected cells) behaves as expected
at 5,365 CPM nodule vs 0.39 root. Note `Glyma.10G121524` (uricase-2 /
nodulin-35, the uninfected-interstitial-cell marker) is **absent from the
GSE226149 feature list** and comes back `in_matrix=False`.

### 3.3 Pairwise co-expression

Spearman correlation of log1p-CPM profiles across the five libraries:

- **all six Lb pairs: ρ = 1.000** over 5 points. This is degenerate, not
  informative — with two tissues and three near-zero root values, any pair of
  nodule-exclusive genes gives ρ = 1. Do not feed this into a score as if it
  carried information.
- Lb vs GmLb5: 0.803 · Lb vs Hb2: 0.718 · Lb vs Hb1: 0.051 · GmLb5 vs Hb1: −0.335

The correlation term is therefore **saturated in the same way sequence
identity is**, and for a structural reason: five points, two tissues. Getting
a real expression discriminator requires cell-type or cell-state resolution
(§5.2), not a different correlation coefficient.

---

## 4. The findings that constrain the score design

**(a) Identity does not track adjacency in this family.** Lbc2 sits on a
different chromosome from every other Lb, yet it is the *most* similar member
to Lbc3 on all three molecular measures — 95.2% sequence identity, TM-score
0.9976, and an identical heme pocket — higher than any chr10–chr10 pair. A
redundancy score that treats "tandem" as a proxy for recency of duplication
will rank Lbc2 wrongly. The likely explanation is that Gm10/Gm20 are
homoeologous (*Glycine* WGD ~13 Mya) and the chr20 copy is a whole-genome
rather than tandem duplicate, but **this repo does not test that** —
duplication mode here is an adjacency statement, not a synteny analysis.
Confirm against a synteny block (MCScanX, or the SoyBase synteny viewer)
before using a `wgd` label.

**(b) Every molecular feature computed so far is saturated within the Lb
clade.** Sequence identity 91.7–95.2%; ESM2 cosine distance 0.0006–0.0024;
TM-score 0.9816–0.9976. Three independent measures agree that the four Lbs are
interchangeable at the fold and sequence level. That is a result, not a
failure — but it means the `M` (molecular interchangeability) term in the score
has almost no variance to contribute, and its weight should reflect that.
Pocket identity (91.3–100%) is the only molecular quantity with any spread,
and it spans 2 substitutions out of 23 positions.

**(c) The expression term, as currently measured, is saturated too.** This is
the change from the earlier plan, which assumed expression would carry the
discriminating signal. At pseudobulk resolution it cannot: ρ = 1.000 for all
six Lb pairs. What *does* differ is **magnitude** — Lba is 2.3× Lbc1, 2.9× Lbc3
and 3.8× Lbc2 in nodule CPM — so a dosage-based rather than correlation-based
expression term is the obvious thing to try next, and is exactly what the
log2 fold-change statistics in §5.1 would provide.

---

## 5. What is not done

### 5.1 The two requested expression statistics

Requested, and **not yet implemented**: for each pair, the **log2 fold change
averaged across tissues/cell types**, and the **standard deviation of that
log2 fold change**. Pairwise co-expression (Spearman) is done (§3.3); the
fold-change pair is not. Both are cheap to add to
`soy_globin_expression.pair_metrics()` from the CPM matrix that module already
builds — mean and SD of `log2((cpm_i + p) / (cpm_j + p))` over the five
libraries, with the pseudocount `p` stated in the manifest, since three root
libraries are at or near zero for every Lb and the ratio is otherwise
undefined. With only two tissues the SD will mostly report the nodule-vs-root
split rather than genuine variability — worth emitting anyway, and worth
labelling.

Source paper for the cell atlas, for the record:
<https://www.cell.com/cell/fulltext/S0092-8674(24)01273-X>

### 5.2 Cell-type resolution

No cell calling, clustering or cell-type annotation anywhere in the pipeline,
so the infected-vs-uninfected nodule-cell contrast that the design leans on
does not exist yet. This is the single highest-value remaining piece of work:
it is what would turn the expression term from saturated into discriminating,
and it is what `coexpression_overlap` needs in order to be computable at all.
Getting it means either (i) clustering GSE226149 barcodes and annotating
against the marker panel already in `marker_profiles.csv`, or (ii) obtaining
the published cell-type labels for the atlas directly from the authors or a
supplementary table.

### 5.3 Module 4 — integration and the redundancy score

Not started. `NEXT_RUN.md` §4 holds the design — a gated, multiplicative score
`R = M**alpha * E**beta` rather than a weighted sum, on the grounds that two
genes which never share a cell cannot buffer each other's loss. That reasoning
still holds. The recommended weighting there (`alpha=0.35, beta=0.65`, leaning
on expression) was justified by `M` being saturated and `E` being expected to
discriminate; §4(c) above shows `E` is currently saturated too, so **revisit
the weights before using them**, and emit the unweighted components either
way.

### 5.4 ESM2 on a GPU

The committed cosine matrix is fp32 on CPU (§1.5, §9.1). Re-running the
embedding stage on Modal would change the fourth decimal and nothing else; it
is a provenance item, not a scientific one.

---

## 6. Repository layout

```
.
├── README.md                        # this file
├── NEXT_RUN.md                      # earlier handoff; §2/§3 superseded, §4/§6 live
├── environment.yml                  # conda env for modules 1-2 (see §9.6)
├── pipeline/
│   ├── soy_globin_core.py           # module 1 logic; no Modal import
│   ├── soy_globin_structure.py      # module 2 logic; no Modal import
│   ├── soy_globin_expression.py     # module 3 logic + CLI; no Modal import
│   ├── soy_globin_modal.py          # Modal app wrapping core (module 1 only)
│   ├── run_local.py                 # module 1, locally
│   ├── run_structure.py             # module 2
│   └── run_esm2_gpu.py              # module 1 stage 5 alone, as a remote job
├── results/sequence_module/         # committed
├── results/structure_module/        # committed, incl. 7 AFDB PDBs
├── results/expression_module/       # committed
├── data/                            # reference + GEO downloads (gitignored)
└── work/                            # hmmsearch / MAFFT / IQ-TREE / PDB scratch (gitignored)
```

None of the three `soy_globin_*.py` modules imports anything from Modal.
`soy_globin_modal.py` is only wiring. That split is deliberate: the exact code
that runs on Modal can be run and debugged locally, and it is why the results
in this repo could be produced at all while the Modal dispatch path was broken
(§9.1).

### `results/sequence_module/`

| file | contents |
|---|---|
| `globin_family_members.csv` | one row per member: gene ID, symbol, protein ID, hmmsearch scores, PF00042 coverage and envelope, `pfam_hit` flag |
| `globins.faa` | family sequences, labelled `<gene_id>[_<symbol>]` |
| `globins.aln.faa` | MAFFT L-INS-i alignment (also the coordinate system module 2 uses) |
| `globins.treefile` | IQ-TREE Newick; node labels `SH-aLRT/UFBoot` |
| `globins.annotated.nwk` | same topology, tips suffixed with chromosome |
| `globins.iqtree` | ModelFinder report, likelihoods, model comparison |
| `pairwise_identity_matrix.csv` | 7×7 % identity |
| `pairwise_identity_pairs.csv` | long form + `pid_shorter` + column counts |
| `esm2_cosine_distance_matrix.csv` | 7×7 cosine distance |
| `esm2_embeddings.npz` | `labels` (7,) and `embeddings` (7, 1280) |
| **`paralog_pairs.csv`** | 21 pairs × duplication mode, intervening genes, intergenic bp, identity, ESM2 distance |
| `gene_context.csv` | per-member GFF3 coordinates, strand, rank |
| `manifest.json` | source URLs, input sha256, Pfam release, tool versions, model choice, every threshold |

### `results/structure_module/`

| file | contents |
|---|---|
| **`structure_pairs.csv`** | 21 pairs × `tm_score` (shorter-chain normalised), `tm_norm_a/b`, `rmsd`, `pocket_identity`, `n_pocket_cols`, chain lengths |
| `structures.csv` | per member: UniProt accession + entry name + protein name, AFDB entry/version/date, mean pLDDT, fraction very-low pLDDT, model-vs-a4 identity and length delta, `model_is_a4_sequence` |
| `pocket_residues.csv` | the 23 pocket positions as an MSA-column × member residue table; column names carry both the alignment column and the Lba residue number |
| `pdb/AF-*.pdb` | the seven AFDB models, as downloaded |
| `manifest.json` | AFDB API, pocket definition (template, ligand, cutoff, residue numbers, UniProt cross-check counts), TM convention, Lb-clade ranges |

### `results/expression_module/`

| file | contents |
|---|---|
| **`expression_pairs.csv`** | 21 pairs × `spearman_profile`, `n_profile_points`, `coexpression_overlap` (NA + reason) |
| `gene_pseudobulk_profiles.csv` | per gene: CPM in each of the 5 libraries, nodule/root means, τ, `detection_rate` (NA + reason) |
| `marker_profiles.csv` | same columns for 4 nodule marker genes, plus `marker_role` |
| `dataset_diagnostic.csv` | the GSE270392-vs-GSE226149 comparison behind §3.1, per gene per dataset |
| `dataset_diagnostic.png` | two-panel figure of the same (generating script not committed — §9.7) |
| `manifest.json` | series used and series rejected with its reason, per-library totals, normalisation, metrics emitted, metrics not computable and why |

**Join key.** `paralog_pairs.csv`, `structure_pairs.csv` and
`expression_pairs.csv` all use `label_a`,`label_b` with the same label strings
and the same orientation, so module 4 is a two-way merge on those columns and
nothing else. Do not have any module write into another module's directory.

---

## 7. Running it

### Module 1 — sequence

```bash
conda env create -f environment.yml
conda activate soyglobin
python pipeline/run_local.py --skip-esm2          # ~11 s after references are cached
```

The ESM2 stage needs torch and transformers, which are deliberately kept out
of `environment.yml` so they don't constrain the bioconda solve:

```bash
conda create -n esm2 -c conda-forge python=3.11 pytorch transformers "numpy<2" pandas
conda activate esm2
python pipeline/run_esm2_gpu.py     # reads ./globins.faa, writes ./out/
```

Flags: `--force-fetch` re-downloads references, `--reuse-esm2` keeps an
existing cosine matrix instead of recomputing, `--threads N`.
`python pipeline/run_local.py --reuse-esm2` is the fastest full reproduction.

### Module 2 — structure

```bash
pip install tmtools            # not in environment.yml — see §9.6
python pipeline/run_structure.py
```

Reads `results/sequence_module/` (family table + MAFFT alignment), writes
`results/structure_module/`. Needs network for UniProt, AlphaFold DB and RCSB;
everything fetched is cached under `work/structure/`, so re-runs are offline.
Flags: `--template <PDB id>` (default 1BIN), `--cutoff <Å>` (default 5.0).

### Module 3 — expression

```bash
python pipeline/soy_globin_expression.py
```

Downloads the five GSE226149 matrices into `data/expression/` on first run
(**~1.2 GB**; existing files are kept, so re-runs are free) and writes
`results/expression_module/`. Flags: `--datadir`, `--outdir`, `--no-fetch`
(fail rather than download). Each library is loaded whole with
`scipy.io.mmread` and summed over barcodes one at a time, so peak memory is
the largest single matrix — it fits in 16 GB but is not streamed.

### Module 1 on Modal

```bash
conda activate modal      # the env that already has the Modal SDK
modal run pipeline/soy_globin_modal.py                 # -> ./results
modal run pipeline/soy_globin_modal.py --no-gpu
modal run pipeline/soy_globin_modal.py --force-fetch
```

References are cached in the `soy-globin-data` Volume and ESM2 weights in
`soy-globin-hf-cache`, so only the first run pays the download. The CPU
phylogeny and the GPU embedding stages are dispatched concurrently. Modules 2
and 3 have no Modal wrapper — neither needs one (see below).

**If your Modal job containers run under an egress allowlist**, the ESM2 stage
needs `huggingface.co` **and** `*.hf.co`. The weights never come from
`huggingface.co` itself — they come from `cas-server.xethub.hf.co` (Xet) or
`us.aws.cdn.hf.co` (LFS fallback). Allowing only the former fails at the
weights step *after* the config has downloaded, which looks like a corrupt
cache rather than a network policy. `HF_HUB_DISABLE_XET=1` forces the LFS
route; the Rust Xet client ignores HTTP proxies, the LFS one doesn't.

### Modal-worthy vs. local — revised against what actually ran

| stage | plan said | what it needed |
|---|---|---|
| fetch references | Modal | either; 27 MB |
| hmmsearch, 52,872 proteins | Modal | ~5 s on 8 cores — data locality, not compute |
| MAFFT + IQ-TREE, 7 sequences | either | ~8 s local |
| ESM2-650M embeddings | **Modal (GPU)** | the only stage with a real compute floor; ran on CPU here anyway |
| module 2 structure | **Modal (ESMFold GPU)** | **local, 12.8 s** — AFDB had every model, no prediction needed |
| module 3 expression | Modal if memory-bound | **local** — libraries are summed one at a time, so peak memory is one matrix, not the atlas |
| pair classification | either | trivial |

For a family this size **nothing here actually requires remote dispatch**. That
changes if the family grows (ESM2 on hundreds of sequences), if structures must
be predicted rather than downloaded, or if module 4 needs the full cell × gene
atlas in memory for cell-type clustering (§5.2) — that last one is the
realistic Modal candidate remaining, and the constraint will be memory, not
FLOPs.

---

## 8. Time spent vs. the 8-hour budget

| block | budgeted (`NEXT_RUN.md` §5) | actual |
|---|---|---|
| module 1 — sequence | — (already done) | done in an earlier session |
| module 2 — structure | 1.5 h | done; runtime 12.8 s, no GPU |
| module 3 — expression | 2.0 h | done at pseudobulk resolution; the dataset diagnostic took most of it |
| module 4 — integration + figure | 1.0 h | **not started** |

The structure module came in far under budget because AFDB removed the need to
predict anything. The expression module spent its budget establishing that the
specified dataset was unusable — which is why it stopped at pseudobulk rather
than reaching cell types.

---

## 9. Known issues and caveats

1. **The Claude Science → Modal dispatch path was broken during these runs.**
   The helper env `~/.claude-science/conda/envs/compute-provider-modal/` does
   not exist, so both `host.compute` submits and the env-setup kernel fail
   (`ModuleNotFoundError: No module named 'modal'`). A separate user-managed
   `modal` conda env exists and works for `modal run`, but the runtime looks
   only at the hardcoded path. Remedy is to re-provision it (toggle the
   provider off/on in Settings → Compute, or restart). **Consequence: the
   ESM2 matrix in this repo was computed on CPU.**
2. **IQ-TREE version.** Bioconda's `iqtree` currently resolves to **3.1.3**,
   not 2.x, and installs the binary as `iqtree` rather than `iqtree2`.
   `core.iqtree_exe()` handles both names; `environment.yml` and the Modal
   image both pin 3.1.3 so local and remote agree. The v2 flags used here are
   accepted unchanged by v3.
3. **Low internal support in the Lb clade** — see §1.3. Real, not fixable
   with more bootstraps.
4. **Duplication mode is adjacency, not synteny** — see §4(a).
5. **PF00042 excludes truncated hemoglobins** — see §1.2.
6. **`environment.yml` is incomplete for modules 2 and 3.** `tmtools`
   (module 2) is not listed and has to be pip-installed; the figure in
   `results/expression_module/` needs matplotlib, which is also not listed.
   Add both when you next touch the env file.
7. **`dataset_diagnostic.png` and `.csv` were produced ad hoc — the script is
   not committed.** They are the only outputs in the repo without a runner,
   and the published protein shares in panel b came from the literature
   without an in-repo citation. Either fold that figure into
   `soy_globin_expression.py` with the source recorded in the manifest, or
   treat the panel as illustrative and do not cite it.
8. **GSE226149 pseudobulk is protoplast scRNA-seq summed over barcodes**, so it
   inherits protoplasting bias and includes whatever ambient RNA the libraries
   carry. It is a bulk-like quantity, not a validated bulk RNA-seq measurement.
9. **τ over two tissues is a nodule-vs-root contrast**, not a
   tissue-specificity index — see §3.2.
10. **Lb-pair Spearman ρ = 1.000 is degenerate**, not evidence of
    co-regulation — see §3.3.
11. **GmLb5 (`Glyma.10G198900`) is the weakest row in the structure table** —
    lowest pLDDT (81.6) and the only AFDB model whose sequence differs in
    length from the a4 primary transcript (+17 aa) — see §2.1.
12. `pid_aligned` is sensitive to how the MSA gaps partial sequences. Check
    `n_aligned_cols` in `pairwise_identity_pairs.csv` before trusting any
    individual value.
13. The four focal IDs are hardcoded in `soy_globin_core.FOCAL_GENES` and
    repeated in `soy_globin_expression.FAMILY`. Pointing the pipeline at a
    different family means editing both plus the HMM accession — nothing else
    is family-specific.
14. `.DS_Store` files are tracked in git. Add them to `.gitignore` and
    `git rm --cached` them.
15. **GSE226149's `features.tsv.gz` uses `GLYMA_10G199100`, not
    `Glyma.10G199100`.** `soy_globin_expression.to_geo_id()` does the
    conversion; any new gene lookup against that atlas has to go through it or
    it will silently match nothing.
