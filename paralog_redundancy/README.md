# Paralog redundancy from Ensembl Compara

A single-file pipeline that quantifies **paralog redundancy** for a gene family, implementing the
method of:

> De Kegel B, Ryan CJ (2019) *Paralog buffering contributes to the variable essentiality of genes in
> cancer cell lines.* PLoS Genetics 15(10):e1008466. [doi:10.1371/journal.pgen.1008466](https://doi.org/10.1371/journal.pgen.1008466)

Given a few seed gene IDs, it resolves the paralog family from Ensembl Compara, scores sequence
identity in both directions for every pair, applies the paper's inclusion filter, summarises each
gene by the redundancy available to it, and labels each pair as a whole-genome (WGD) or small-scale
(SSD) duplicate.

Works for any species served by the unified Ensembl REST endpoint — vertebrates and the Ensembl
Genomes divisions (Plants, Fungi, Metazoa, Protists).

## Method

| Step | What the pipeline does |
|---|---|
| 1 | Fetch within-species paralogues for each seed from Compara; optionally expand to the seeds' partners, giving the connected component of the family in the paralogy graph. |
| 2 | Record protein sequence identity **in both directions** for every pair — the percent of A's own sequence matched in B, and of B's in A. These differ whenever the proteins differ in length, and the asymmetry is the point: one symmetric number would assign a pair a single redundancy level instead of two. |
| 3 | Keep pairs sharing `--min-identity` (default 20%) **in both directions** where both genes are protein-coding. |
| 4 | Summarise each gene by the number of pairs it belongs to and the **maximum percent of its own sequence matched in any paralogue**. Genes absent from the summary are singletons. |
| 5 | Label pairs WGD or SSD; a gene is WGD if it belongs to *any* WGD pair. |

**On step 5.** The original consumes two curated human ohnolog lists, which do not exist for most
species. Those lists are themselves built by synteny comparison, so this pipeline reproduces the
criterion directly: for a cross-chromosome pair it collects genes in a `--window` interval around
each locus, fetches paralogues for one window, and calls a homeologous block when enough reciprocal
anchors are found and their gene order is collinear. Passing `--wgd-node` additionally requires
Compara's inferred duplication node for the pair to be one of the labels you list, which prevents
out-paralogs that merely *sit* in a homeologous block from being labelled ohnologs. Same-chromosome
pairs are SSD by definition.

The node test is **exact string matching**, not a taxonomic-depth comparison: Compara node labels
carry no ordering in the REST response, so the script cannot work out on its own which of two
labels is older. List every node label you consider WGD-era — they are printed in the
`duplication_node` column of `paralog_pairs.csv`, so a first pass with `--no-synteny` tells you
which labels your family actually uses.

Both conditions are written to `paralog_pairs.csv` as their own columns —
`in_homeologous_block` (positional) and `duplication_node_concordant` (temporal) — with `dup_mode`
being WGD only where both hold. Read `dup_mode` for the classification; a row can legitimately show
`in_homeologous_block = True` with `dup_mode = SSD`.

## Install

```bash
pip install -r requirements.txt   # pandas, numpy, matplotlib
```

`curl` must be on `PATH` — it is used instead of a Python HTTP client because some environments
route outbound traffic through a proxy that `requests`/`urllib` cannot reach.

## Usage

```bash
python paralog_redundancy.py \
    --species glycine_max \
    --genes GLYMA_10G198800 GLYMA_10G198900 GLYMA_10G199000 GLYMA_10G199100 GLYMA_20G191200 \
    --wgd-node Glycine_subgen._Soja \
    --out results/
```

Useful flags:

| Flag | Effect |
|---|---|
| `--min-identity` | inclusion floor in percent, required in both directions (paper: 20) |
| `--window` | half-width in bp of the synteny window (default 400 kb) |
| `--wgd-node` | one or more Compara duplication-node labels considered WGD-era, matched exactly; a pair at any other node stays SSD |
| `--no-synteny` | skip the synteny test entirely — every pair becomes SSD. Fast |
| `--no-expand` | query Compara for the seed genes only, not their partners |
| `--cache` | HTTP cache directory (default `.ensembl_cache`) |

**Runtime.** The Compara homology endpoint takes roughly 10–15 s per gene, and the synteny test
queries every gene in one window, so a first run with synteny on an ~80-gene window takes
20–30 minutes. Every response is cached on disk, so re-runs and parameter sweeps are near-instant,
and `--no-synteny` finishes in under a minute. The endpoint also drops connections intermittently
(`SSL_read: unexpected eof`); `fetch()` retries and validates each body, so a run recovers without
manual intervention.

## Outputs

Written to `--out`:

| File | Contents |
|---|---|
| `gene_redundancy.csv` | per-gene features: pair count, max/mean/min percent of own sequence matched, closest paralogue, duplication mode |
| `paralog_pairs.csv` | every pair with identity in both directions, duplication node, arrangement, WGD/SSD call, and whether it passed the filter |
| `synteny_anchors_chrX_chrY.csv` | the reciprocal anchor relations supporting each block call |
| `run_metadata.json` | parameters, assembly, Ensembl endpoint, block statistics, timestamp |
| `paralog_redundancy.png` | three panels: the synteny block, the asymmetric identity matrix, per-gene redundancy |

## Worked example: soybean leghemoglobins

`examples/soybean_leghemoglobin/` holds the full output for the *Glycine max* leghemoglobin family,
generated by the command above. Summary:

* Compara places the five leghemoglobins in a single closed clique — every gene in 4 pairs, all 10
  pairs clearing the filter, no singletons.
* *LB1*, *LB4*, *LB2*, *LB3* each retain a backup copy matching 93–95% of their own sequence; *LB5*
  is the outlier at 69.0%.
* *LB5* encodes a 168-aa protein against 144–145 aa for the rest, so 80.0% of *LB4*'s sequence is
  matched in *LB5* while only 69.0% of *LB5*'s is matched in *LB4* — the bidirectional step earning
  its keep.
* 67 of 83 genes flanking *LB4* on chromosome 20 have a paralogue in the matching chromosome-10
  window, collinear at Spearman ρ = −0.95 (inverted orientation), establishing the homeologous
  block. Final split: 3 WGD pairs (*LB4* against *LB1*, *LB2*, *LB3*) and 7 SSD pairs, with *LB5*
  the only SSD-only gene.

`docs/METHODS_soybean_leghemoglobin.md` is the full write-up, including how each step of the paper
maps onto soybean and what the analysis does *not* cover.

## Scope

This implements the **redundancy-quantification half** of the source study. The other half — the
association between these features and variable gene essentiality, and the synthetic-lethality
inference — needs loss-of-function fitness screens across many genetic backgrounds. Where no such
resource exists, these features describe redundancy *potential*, not measured functional
compensation.

Family membership is Compara's, exactly as the source study takes it. A family defined by a
domain-wide HMM or a custom gene tree may be broader.
