# Soybean leghemoglobin redundancy — Boston Computational Biology Hackathon

Predicting ortholog/paralog **functional redundancy** in soybean (*Glycine max*),
using the leghemoglobin family as the test case. The full design is four
modules — sequence, structure, expression, integration — feeding one weighted
redundancy score.

**Status: the sequence module is built, run, and its outputs are in
`results/sequence_module/`. The other three modules are not started.**
`NEXT_RUN.md` is the handoff for whoever (or whatever) picks this up next.

Focal genes (Wm82.a4.v1 IDs):

| symbol | gene ID | chromosome |
|---|---|---|
| Lba | `Glyma.10G199100` | Gm10 |
| Lbc1 | `Glyma.10G199000` | Gm10 |
| Lbc2 | `Glyma.20G191200` | Gm20 |
| Lbc3 | `Glyma.10G198800` | Gm10 |

---

## 1. What was actually done

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

**This stage ran on CPU, not on a GPU** (see §5). The Modal path uses fp16 on
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
| `Glyma.10G198900` | — | 43,085,930 | 43,089,182 | 1837 |
| `Glyma.10G199000` | Lbc1 | 43,090,608 | 43,091,995 | 1838 |
| `Glyma.10G199100` | Lba | 43,095,323 | 43,097,067 | 1839 |

`Glyma.20G191200` (Lbc2) is the single Gm20 copy at 42,955,083–42,956,301.

---

## 2. The two findings that should change the score design

**(a) Identity does not track adjacency in this family.** Lbc2 sits on a
different chromosome from every other Lb, yet it is the *most* identical
sequence to Lbc3 (95.2%) — higher than any chr10–chr10 pair. A redundancy
score that treats "tandem" as a proxy for recency of duplication will rank
Lbc2 wrongly. The likely explanation is that Gm10/Gm20 are homoeologous
(*Glycine* WGD ~13 Mya) and the chr20 copy is a whole-genome rather than
tandem duplicate, but **this repo does not test that** — duplication mode
here is an adjacency statement, not a synteny analysis. Confirm against a
synteny block (MCScanX, or the SoyBase synteny viewer) before using a `wgd`
label.

**(b) Both sequence features are saturated within the Lb clade.** 91.7–95.2%
identity and cosine distances in the third-to-fourth decimal. Neither
discriminates among Lba/Lbc1/Lbc2/Lbc3 — they only separate the Lbs from the
phytoglobins, which we already knew. The sequence module is therefore good at
*defining the family* and weak as a *redundancy discriminator*, and the weight
it gets in the integrated score should reflect that. Expression (especially
the infected-cell fraction of the nodule atlas) and structure have to carry
the discriminating signal.

---

## 3. Repository layout

```
.
├── README.md                      # this file
├── NEXT_RUN.md                    # handoff / plan for the next session
├── environment.yml                # conda env for everything except ESM2
├── pipeline/
│   ├── soy_globin_core.py         # all pipeline logic; no Modal import
│   ├── soy_globin_modal.py        # Modal app wrapping core
│   ├── run_local.py               # run the whole module locally
│   └── run_esm2_gpu.py            # stage 5 alone, as a remote job
├── results/sequence_module/       # committed outputs
├── data/                          # reference downloads (gitignored)
└── work/                          # hmmsearch / MAFFT / IQ-TREE scratch (gitignored)
```

All pipeline logic lives in `soy_globin_core.py` and it imports nothing from
Modal. `soy_globin_modal.py` is only wiring. That split is deliberate: it
means the exact code that runs on Modal can be run and debugged locally, and
it is why the results in this repo could be produced at all while the Modal
dispatch path was broken.

### Outputs in `results/sequence_module/`

| file | contents |
|---|---|
| `globin_family_members.csv` | one row per member: gene ID, symbol, protein ID, hmmsearch scores, PF00042 coverage and envelope, `pfam_hit` flag |
| `globins.faa` | family sequences, labelled `<gene_id>[_<symbol>]` |
| `globins.aln.faa` | MAFFT L-INS-i alignment |
| `globins.treefile` | IQ-TREE Newick; node labels `SH-aLRT/UFBoot` |
| `globins.annotated.nwk` | same topology, tips suffixed with chromosome |
| `globins.iqtree` | ModelFinder report, likelihoods, model comparison |
| `pairwise_identity_matrix.csv` | 7×7 % identity |
| `pairwise_identity_pairs.csv` | long form + `pid_shorter` + column counts |
| `esm2_cosine_distance_matrix.csv` | 7×7 cosine distance |
| `esm2_embeddings.npz` | `labels` (7,) and `embeddings` (7, 1280) |
| **`paralog_pairs.csv`** | **the integration table** — 21 pairs × duplication mode, intervening genes, intergenic bp, identity, ESM2 distance |
| `gene_context.csv` | per-member GFF3 coordinates, strand, rank |
| `manifest.json` | source URLs, input sha256, Pfam release, tool versions, model choice, every threshold |

`paralog_pairs.csv` is the file the integration module should consume.

---

## 4. Running it

### Locally

```bash
conda env create -f environment.yml
conda activate soyglobin
python pipeline/run_local.py --skip-esm2          # ~8 s after references are cached
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

### On Modal

```bash
conda activate modal      # the env that already has the Modal SDK
modal run pipeline/soy_globin_modal.py                 # -> ./results
modal run pipeline/soy_globin_modal.py --no-gpu
modal run pipeline/soy_globin_modal.py --force-fetch
```

References are cached in the `soy-globin-data` Volume and ESM2 weights in
`soy-globin-hf-cache`, so only the first run pays the download. The CPU
phylogeny and the GPU embedding stages are dispatched concurrently.

**If your Modal job containers run under an egress allowlist**, the ESM2 stage
needs `huggingface.co` **and** `*.hf.co`. The weights never come from
`huggingface.co` itself — they come from `cas-server.xethub.hf.co` (Xet) or
`us.aws.cdn.hf.co` (LFS fallback). Allowing only the former fails at the
weights step *after* the config has downloaded, which looks like a corrupt
cache rather than a network policy. `HF_HUB_DISABLE_XET=1` forces the LFS
route; the Rust Xet client ignores HTTP proxies, the LFS one doesn't.

### Modal-worthy vs. local

| stage | where | why |
|---|---|---|
| fetch references | Modal | 27 MB lands next to the Volume instead of crossing the wire twice |
| hmmsearch, 52,872 proteins | Modal | ~5 s on 8 cores — dispatched for data locality, not compute |
| MAFFT + IQ-TREE | either | ~8 s for 7 sequences |
| **ESM2-650M embeddings** | **Modal (GPU)** | the only stage with a real compute floor |
| pair classification | either | trivial |

For a family this size the GPU stage is the sole justification for remote
dispatch. That changes for the structure module (ESMFold) and the atlas
module, where memory rather than FLOPs is the binding constraint.

---

## 5. Known issues and caveats

1. **The Claude Science → Modal dispatch path was broken during this run.**
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
4. **Duplication mode is adjacency, not synteny** — see §2(a).
5. **PF00042 excludes truncated hemoglobins** — see §1.2.
6. `pid_aligned` is sensitive to how the MSA gaps partial sequences. Check
   `n_aligned_cols` in `pairwise_identity_pairs.csv` before trusting any
   individual value.
7. The four focal IDs are hardcoded in `soy_globin_core.FOCAL_GENES`. Point
   the module at a different family by editing that dict and the HMM
   accession — nothing else is family-specific.
