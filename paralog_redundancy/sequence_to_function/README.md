# Sequence-to-function divergence of the leghemoglobin family

Follow-on to the identity-based redundancy analysis in the parent directory. Two foundation
models are applied to the same five genes, and both are compared against the
De Kegel & Ryan 2019 percent-identity features.

| | ESM-2 (protein) | Evo 2 (DNA) |
|---|---|---|
| Input | 145-168 aa protein | 435-507 nt CDS |
| What its space encodes | fold and active-site similarity | which genomic sequence is plausible |
| Sees synonymous change | **no, by construction** | yes |
| Sees codon usage | no | yes |
| Sees genomic position | no | yes (1 Mb context) |
| Per-residue attribution | yes (masked marginals) | no (sequence likelihood) |
| Status here | **run, results included** | **job ready, not yet run** |

## ESM-2 — complete

`esm_divergence.py` (CPU, ~8 min for this family with `esm2_t33_650M_UR50D`) computes two things:
mean-pooled embedding cosine distance per pair, and the masked-marginal log-odds
(Meier et al. 2021) at the positions where two aligned paralogs differ, evaluated in each
protein's own context.

Headline result: **the embedding distance is percent identity restated** (gene-level Spearman
rho = -1.000), so a protein language model adds nothing to the redundancy call in a family this
similar. The constraint score, however, is **orthogonal to identity** (gene-level rho = +0.05) and
exposes disagreements identity cannot — notably LB2-LB4, which is 93.1% identical yet carries
per-substitution severity comparable to the 67-69% identical LB5 pairs, and LB4, which is every
other gene's nearest neighbour by identity but has the most negative own-context log-odds in the
family. Full account with caveats in `METHODS_esm_vs_identity.md`.

Reproduce:

```bash
pip install fair-esm torch
# weights: dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt (+ -contact-regression.pt)
#          into $TORCH_HOME/hub/checkpoints/
python esm_divergence.py                      # ESM_MODEL=esm2_t6_8M_UR50D for a fast smoke test
```

## Evo 2 — job ready, results NOT in this directory

Evo 2 needs >=24 GB CUDA VRAM for the 7B model and has no CPU path, so it was not run in the
environment that produced everything else here. **`results/` therefore contains no Evo 2
scores — they are absent because the job has not been run, not because they were lost.**

`evo2_sequences.json` holds the 65 sequences to score: the 5 reference CDS plus, for each ordered
pair in both directions, three chimeras built by swapping aligned codons from the donor —
`all` (every differing codon), `syn` (synonymous only) and `nonsyn` (amino-acid-changing only).
Burden is then `ll(chimera) - ll(reference)` from `Evo2.score_sequences`.

The `syn` chimeras are the control that makes the two models directly comparable: they translate
to a byte-identical protein, so ESM-2 must score them as unchanged, while Evo 2 still assigns a
likelihood shift.

```bash
conda activate modal
modal run evo2_modal_app.py --sequences evo2_sequences.json     # writes evo2_scores.json
modal run evo2_modal_app.py --sequences evo2_sequences.json --model evo2_1b_base   # cheap smoke test
```

First image build compiles flash-attn (20-40 min); if it fails, install flash-attn separately with
`--no-build-isolation` before `evo2`. First run pulls ~15 GB of weights into a Modal Volume that
later runs reuse. Defaults to A100-40GB. The local entrypoint's reporting was dry-run tested
against synthetic scores; the GPU path itself is untested.

## Why the DNA model can see more here

Codon-level partition of every pairwise difference (`results/leghemoglobin_codon_divergence.csv`),
computed from the canonical-transcript CDS in `cds/`:

| pair    |   protein % id |   codon % id |   syn |   nonsyn |   nt subs |   frac syn |
|:--------|---------------:|-------------:|------:|---------:|----------:|-----------:|
| LB1-LB4 |          95.17 |        84.83 |    15 |        7 |        24 |     0.6818 |
| LB2-LB4 |          93.75 |        85.42 |    12 |        9 |        22 |     0.5714 |
| LB3-LB4 |          93.06 |        87.5  |     8 |       10 |        21 |     0.4444 |
| LB1-LB3 |          92.36 |        86.11 |     9 |       11 |        23 |     0.45   |
| LB2-LB3 |          92.36 |        87.5  |     7 |       11 |        19 |     0.3889 |
| LB1-LB2 |          91.67 |        81.94 |    14 |       12 |        28 |     0.5385 |
| LB5-LB4 |          80    |        64.83 |    22 |       29 |        57 |     0.4314 |
| LB5-LB2 |          79.86 |        63.19 |    24 |       29 |        57 |     0.4528 |
| LB5-LB3 |          77.78 |        63.19 |    21 |       32 |        60 |     0.3962 |
| LB1-LB5 |          77.24 |        61.38 |    23 |       33 |        65 |     0.4107 |

Across the 10 pairs, **155 of 338 changed codons (45.9%) are synonymous** — nearly
half the divergence in this family sits outside a protein model's field of view. The effect is
largest exactly where the identity analysis saw least divergence: LB1 and LB4, the 95.2%-identical
WGD ohnolog pair that the De Kegel features call the family's most redundant, differ at 22 codons
and 24 nucleotide positions, two-thirds of those changes silent.

## Layout

```
esm_divergence.py                     ESM-2 analysis (CPU)
evo2_modal_app.py                     Evo 2 scoring job (Modal GPU)
evo2_sequences.json                   65 sequences: 5 references + 60 chimeras
METHODS_esm_vs_identity.md            full ESM-2 write-up, results and caveats
cds/                                  canonical-transcript CDS, verified ATG..stop
results/                              ESM-2 outputs, the joined comparison tables, the figure
```
