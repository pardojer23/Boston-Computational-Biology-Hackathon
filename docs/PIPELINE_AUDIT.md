# Step-by-step audit — what each step does, what it actually depends on

Companion to `pipeline_dag_current.png`. One block per executable step, as the
code stands at commit `a2dfe7e`. Each block reads:

- **does** — the computation
- **in** / **out** — files and in-memory objects crossing the step boundary
- **params** — everything that changes the output
- **cost / env / net** — wall time, which conda env, whether it touches the network
- **refine** — what has to change for this to be a pipeline step rather than a
  line inside a runner. `—` means nothing; the step is already clean.

"Undeclared dependency" throughout means: the step genuinely needs an upstream
output, but resolves it through a module-level constant built from
`Path(__file__).parent.parent`, so it cannot be redirected, overridden in a test,
or checked for staleness.

---

## Module 1 — sequence (`pipeline/run_local.py`, logic in `soy_globin_core.py`)

### 1.1 fetch references — `core.fetch_inputs`
- **does** download the Wm82.gnm4.ann1.T8TQ primary proteome and gene-model GFF3
  from SoyBase, and PF00042 from the InterPro API; gunzip the HMM.
- **in** three URLs (`core.PROTEOME_URL`, `GFF3_URL`, `PFAM_HMM_URL`).
- **out** `data/proteins_primary.faa.gz`, `data/gene_models_exons.gff3.gz`,
  `data/PF00042.hmm` (+ `.gz`).
- **params** `force` (re-download).
- **cost** 27 MB, seconds. **net** yes. **env** soyglobin.
- **refine** cache invalidation is file *presence*, not checksum, and the sha256s
  are computed later in `run_local` for the manifest rather than compared to a
  pinned expectation. Pin the three sha256s in config and fail on mismatch —
  SoyBase re-releases under the same URL. Also: this is the only step whose
  inputs are recorded anywhere, and only for module 1.

### 1.2 family selection — `core.select_family`
- **does** `hmmsearch --cut_ga` PF00042 over 52,872 proteins, union with
  `core.FOCAL_GENES`, build the member table and the labelled FASTA.
- **in** proteome, HMM (1.1).
- **out** `work/hmmer/PF00042.{tblout,domtblout}`; in-memory `seqs`, `members`
  → `results/sequence_module/{globins.faa, globin_family_members.csv}`.
- **params** `--cut_ga` (fixed, not exposed), `cpu`, `FOCAL_GENES`, HMM accession.
- **cost** ~5 s on 8 cores. **net** no.
- **refine** the relaxed `-E 10` sanity search that established the ~50-bit gap
  below the weakest real hit (README §1.2) was run by hand and is not in the
  code. It is the evidence that the family is not truncated, so it belongs in
  the step as a recorded check, not in prose. The threshold should be a param
  (`--cut_ga` vs an E-value) since changing the family means changing it.

### 1.3a alignment — `core.run_mafft`
- **does** MAFFT L-INS-i, 162 columns.
- **in** `globins.faa` (1.2). **out** `work/phylo/globins.aln.faa`.
- **params** `--localpair --maxiterate 1000` (hardcoded), threads.
- **cost** ~1 s. **refine** this alignment is the coordinate system for 1.4, 2.4b
  and 2.5. It is the highest-fan-out artefact in the repo and its mode is
  hardcoded inside the function. Make it a declared parameter so a change is
  visible in the manifest.

### 1.3b phylogeny — `core.run_iqtree`
- **does** ModelFinder + 1000 UFBoot + 1000 SH-aLRT, seed 20240601.
- **in** the MSA. **out** `work/phylo/globins.{treefile,iqtree,contree,...}`.
- **params** `-m MFP -B 1000 --alrt 1000 --seed 20240601`, threads.
- **cost** ~7 s. **refine** nothing downstream consumes the tree — 1.7 relabels
  it and it is never read again. It is a *terminal* branch of the DAG, not a
  dependency of the score. Worth stating explicitly so nobody assumes the score
  is phylogeny-aware. Internal Lb support is 31–47 (README §1.3); keep it as a
  reported figure, not an input.

### 1.4 pairwise identity — `core.pairwise_identity`
- **does** % identity over columns where both sequences have a residue, from the
  MSA (so consistent with the tree), plus `pid_shorter`.
- **in** the MSA. **out** `pairwise_identity_{matrix,pairs}.csv`, long form to 1.6.
- **cost** ms. **refine** —

### 1.5 ESM2 embeddings — `core.esm2_embeddings`
- **does** `esm2_t33_650M_UR50D`, final layer, mean-pooled over residue tokens
  (BOS/EOS/PAD masked).
- **in** `seqs` from 1.2 — **not** the MSA.
- **out** `esm2_cosine_distance_matrix.csv`, `esm2_embeddings.npz`.
- **params** model id, batch size, dtype (fp32 CPU vs fp16 GPU).
- **cost** ~10 min CPU incl. download; seconds on a GPU. **env** `esm2`
  (torch + transformers), **not** soyglobin. **net** yes (HF weights).
- **refine** this is the one genuine environment seam in the pipeline, and the
  runner handles it with three mutually exclusive flags
  (`--skip-esm2` / `--reuse-esm2` / compute). `--reuse-esm2` reads a matrix from
  the *output* directory, relabels it by gene ID, and records
  `"provenance": "reused from a previous run (see that run's manifest)"` — i.e.
  the committed values have no recoverable provenance beyond that sentence.
  Make it a first-class step with its own declared env and a real output hash;
  drop `--skip` (a missing input should be the orchestrator's problem, and
  skipping silently propagates NaN — see 4.3).

### 1.6 genomic context + pair classification — `core.parse_gff_genes`, `classify_pairs`
- **does** rank genes per seqid by start; call each pair tandem (≤10 intervening),
  proximal (same seqid, ≤1 Mb), or dispersed; attach identity and ESM2 distance.
- **in** GFF3 (1.1), members (1.2), identity long (1.4), cosine matrix (1.5,
  optional).
- **out** `paralog_pairs.csv` (21 rows — the module's contract with module 4),
  `gene_context.csv`.
- **params** `TANDEM_MAX_INTERVENING=10`, `PROXIMAL_MAX_BP=1e6`.
- **refine** if 1.5 was skipped, `esm2_cosine_distance` is absent/NaN and this
  table still validates as 21 rows. The failure surfaces four steps later as a
  NaN `M_seq`. Make the ESM2 column a declared requirement of the table schema.

### 1.7 annotated tree — `core.annotate_newick`
- **does** append chromosome to tip labels.
- **in** treefile (1.3b), `gene_context.csv` (1.6). **out** `globins.annotated.nwk`.
- **refine** terminal, presentational. Fine as-is; belongs in a `report` stage
  rather than the sequence stage.

---

## Module 2 — structure (`pipeline/run_structure.py`, logic in `soy_globin_structure.py`)

### 2.1 gene → UniProt — `struc.resolve_uniprot`
- **does** UniProt search API, taxon 3847, one query per gene.
- **in** family gene list, read from `results/sequence_module/globin_family_members.csv`
  via a module-level `SEQ_DIR` constant.
- **out** in-memory accession table → columns of `structures.csv`.
- **cost** seconds. **net** yes, **uncached** — every run re-queries UniProt.
- **refine** the accession mapping is the least stable thing in the pipeline (a
  UniProt release can change the best hit) and it is the only network step with
  no cache. Persist `gene_id → accession` as a versioned output with the UniProt
  release stamped on it, and let the resolution be overridable by a pinned table
  in config. Raise-on-missing is already correct.

### 2.2 AlphaFold DB models — `struc.fetch_afdb`
- **does** `alphafold.ebi.ac.uk/api/prediction/<acc>`, download the PDB the API
  hands back, compute mean pLDDT.
- **in** accessions (2.1). **out** `results/structure_module/pdb/AF-*.pdb` (7).
- **net** yes; cached **by presence in the output directory**.
- **refine** output dir doubling as download cache is the one place where "delete
  the results and re-run" also deletes the cache, and where a stale model
  survives a config change. Move the cache to `work/structure/afdb/` and copy
  into `results/`. Record `afdb_model_version` per file as a dependency — these
  are model_v6, created 2025-08-01; a v7 release would silently change TM-scores.

### 2.3 AFDB-vs-a4 reconciliation
- **does** align each model's chain sequence to the a4 primary transcript, record
  identity, length delta, `model_is_a4_sequence`.
- **in** PDBs (2.2), `globins.faa` (1.2). **out** `structures.csv`.
- **refine** this is a *check* with no threshold: GmLb5's model is +17 aa and the
  run proceeds. Give it a declared tolerance and a warn/fail policy, so the
  caveat (README §2.1, §9.11) is enforced rather than narrated.

### 2.4a heme pocket definition — `struc.heme_pocket_residues`
- **does** residues with any heavy atom ≤5.0 Å of any heme heavy atom in 1BIN.
- **in** RCSB `1BIN.pdb`, cached at `work/structure/1BIN.pdb`.
- **out** 23 residue numbers (in-memory).
- **params** `--template 1BIN`, `--cutoff 5.0`, ligand `HEM`.
- **refine** independent of everything upstream — it depends only on the template
  and the cutoff. In a DAG it is a **root**, not a child of module 1, and should
  be cached on `(template, ligand, cutoff)` so a cutoff sweep is cheap.

### 2.4b pocket transfer onto the MSA
- **does** map the 23 crystal positions to alignment columns through the Lba row
  of the MSA; emit the column × member residue table.
- **in** pocket resnums (2.4a), MSA (1.3a). **out** `pocket_residues.csv`, `cols`.
- **refine** **family hardcode**: the reference row is found with
  `next(l for l in labels if l.endswith("_Lba"))`. Rename the symbol or point the
  pipeline at another family and this raises `StopIteration` with no message.
  The template's reference member belongs in config next to the template PDB id.

### 2.4c UniProt ligand-site cross-check
- **does** fetch P02238's annotated binding sites, map them through the same
  alignment, count how many fall inside the crystal pocket (4 of 5).
- **in** accessions (2.1), cols (2.4b). **out** two integers in the manifest.
- **refine** a validation step whose result is recorded but not asserted. It is a
  genuine check — give it a pass condition.

### 2.5 all-vs-all superposition
- **does** `tmtools.tm_align` over 21 pairs; TM-score normalised by the shorter
  chain; pocket identity over the mapped columns.
- **in** PDB chains (2.2), pocket cols (2.4b), canonical pair order (`core`).
- **out** `structure_pairs.csv` (21 rows — contract with module 4).
- **cost** ~1 s. **refine** —

---

## Module 3 — expression (`pipeline/soy_globin_expression.py`, module and CLI in one file)

### 3.1 fetch GSE226149 — `expr.fetch_inputs`
- **does** download 5 libraries × 3 CellRanger files from the GEO FTP under
  their GEO filenames.
- **out** `data/expression/*` (~1.2 GB).
- **params** `SAMPLES` (accession → stem, tissue).
- **net** yes; cached by presence. **refine** no checksums on 1.2 GB of input;
  a truncated download is detected only by `mmread` failing later. Record sizes
  or md5 from the GEO SOFT record.

### 3.2 pseudobulk — `expr.pseudobulk_library`, `build_pseudobulk`
- **does** `scipy.io.mmread` each matrix whole, sum over barcodes, assemble
  genes × 5 libraries, cache, then CPM per library.
- **in** the five matrices + feature lists (3.1).
- **out** `work/expression/pseudobulk_counts.csv.gz` (0.6 MB), CPM in memory.
- **cost** minutes on first run; peak memory = largest single matrix (160 MB gz,
  fits 16 GB, not streamed). **refine** the cache is validated against
  `list(SAMPLES)` only — change a library's *file* and the stale cache is reused
  silently. Key the cache on the input file checksums. This is also the step that
  would be replaced wholesale by cell-level quantification (README §5.1), so it
  should be the DAG's designated swap point: `counts` and `cpm` are the only
  things downstream consumes.

### 3.3 gene profiles — `expr.gene_profiles`
- **does** per-gene CPM in each library, nodule/root means, τ; `detection_rate`
  emitted as NA with a reason string.
- **in** CPM (3.2) and **the family definition, read from module 1's
  `globin_family_members.csv` through the module-level `SEQ_MEMBERS` constant** —
  an undeclared cross-module dependency (red edge in the DAG).
- **out** `gene_pseudobulk_profiles.csv` — consumed by 4.2, not by 4.1.
- **refine** pass `members_csv` down from the CLI; `load_family()` already takes
  the argument, the caller in `run()` just doesn't use it. τ over two tissue
  means is a nodule-vs-root contrast (README §9.9) — rename the column so the
  caveat travels with the number instead of living in prose.

### 3.4 pair metrics — `expr.pair_metrics`
- **does** Spearman over 5 libraries, log2FC mean/SD (all and nodule-only),
  `dose_ratio`, directional `cover_a_by_b` / `cover_b_by_a`;
  `coexpression_overlap` as NA + reason.
- **in** CPM (3.2), family (same undeclared path as 3.3).
- **out** `expression_pairs.csv` (21 rows — contract with module 4).
- **params** `LOG2FC_PSEUDOCOUNT_CPM=1.0`, `NODULE_LIBS`/`ROOT_LIBS` derived from
  `SAMPLES`.
- **refine** emits `spearman_profile` = 1.000 for all six focal pairs, which is
  degenerate and deliberately unused downstream (README §3.3, §9.10). A column
  that is known-uninformative should be marked as such in the schema, not left
  for a reader to notice in §9. The NA-with-reason convention here is the right
  pattern and should be the repo-wide standard.

### 3.5 marker profiles
- **does** the same profile columns for four nodule markers.
- **in** CPM (3.2), `MARKERS` constant. **out** `marker_profiles.csv`.
- **refine** `Glyma.10G121524` (uricase-2, the uninfected-interstitial marker) is
  absent from the GSE226149 feature list and returns `in_matrix=False`. The
  marker panel is the annotation basis for the cell-type work in README §5.1, so
  the gap needs a substitute chosen *before* that step, not after.

### 3.6 dataset diagnostic — **no runner committed**
- **does** the GSE270392-vs-GSE226149 comparison behind the dataset rejection.
- **out** `dataset_diagnostic.{csv,png}` — committed, with no generator
  (README §9.7). Panel b also plots published protein shares with no in-repo
  citation.
- **refine** the two orphan outputs in the repo. Either fold into this module as
  a real step with the literature source recorded in the manifest, or drop the
  panel. As it stands the *evidence for the module's central design decision* is
  the one thing that cannot be reproduced.

---

## Module 4 — integration (`run_integration.py`, `soy_globin_integration.py`, `plot_integration.py`)

### 4.1 join — `integ.load_pairs`
- **does** outer-merge the three pair tables on `(label_a, label_b)`; verify
  canonical orientation; **raise** on any non-`both` row.
- **in** `paralog_pairs.csv` (1.6), `structure_pairs.csv` (2.5),
  `expression_pairs.csv` (3.4). **out** 21-row joined frame.
- **refine** this is the strongest step in the repo and the model for everything
  else: an explicit contract that fails loudly. Generalise it into a shared
  `contracts.py` — declared columns, dtypes, row count, orientation — applied at
  every module boundary, not just this one.

### 4.2 components — `integ.add_components`, `tissue_overlap`
- **does** `esm2_similarity = -cosine`; histogram-intersection tissue overlap
  from the gene profiles; `pair_class` from the sequence module's `both_focal`.
- **in** joined pairs (4.1), `gene_pseudobulk_profiles.csv` (3.3).
- **refine** **family hardcode**: `pair_class` literals `"Lb-Lb"` / `"Lb-other"`
  / `"other-other"`, and `TISSUE_COLS = {"nodule": ..., "root": ...}` names the
  two tissues. Both should be derived (`focal-focal`, and tissue columns from the
  expression manifest).

### 4.3 score — `integ.score`
- **does** min-max the molecular features over all 21 pairs; three equally
  weighted sub-axes (sequence / fold / pocket) → `M`; `E = tissue_overlap ×
  dose_ratio`; `R = M^α · E^β`; both directional variants; scopes `family` and
  `clade` (floored at 0.05).
- **in** 4.2. **out** `redundancy_scores.csv`, 21 rows × 32 columns.
- **params** `ALPHA_M=0.40`, `M_AXIS_WEIGHTS`, `CLADE_FLOOR=0.05`.
- **cost** ms. **refine** `--alpha` mutates the module global
  `integ.ALPHA_M` from the runner before calling `run()`. It works, and α does
  land in the output CSV, but it is the only parameter in the repo that is set by
  reassigning another module's global. Thread it through as an argument. The
  min-max scope being *all 21 pairs* means adding or removing one family member
  rescales every score — that is correct behaviour, and it is exactly why
  `FOCAL_GENES` needs to be a config input with a recorded hash rather than a
  source constant.

### 4.4 weight sensitivity
- **does** sweep α over 21 steps, record top pair and spread.
- **in** 4.3. **out** `weight_sensitivity.csv`. **refine** —

### 4.5 validation checks
- **does** three checks with evidence into the manifest: focal/outgroup
  separation, directional coverage vs measured abundance, which pair ranks top.
- **in** 4.3, gene profiles. **out** `manifest.json.validation_checks`.
- **refine** these are recorded but non-blocking — `run_integration.py` prints
  `passed=` and exits 0 either way. Decide per check whether it gates the run.
  This is the only place in the pipeline where a scientific assertion is
  evaluated at all.

### 4.6 figure — `plot_integration.py`
- **does** two panels; prints its own text-bbox overlap report.
- **in** `redundancy_scores.csv`, `weight_sensitivity.csv`.
- **out** `redundancy_summary.png`. **refine** —

---

## Cross-cutting findings

1. **No orchestrator, and the run order exists only in README §7.** Four runners
   invoked by hand. Nothing knows that `redundancy_scores.csv` is older than the
   `paralog_pairs.csv` it was computed from. Changing `core.FOCAL_GENES` and
   re-running only module 4 produces a clean-looking, wrong table.
2. **Provenance is per-module and does not chain.** Four independent
   `manifest.json`s; only module 1 records input checksums; no manifest records
   the *identity of the upstream outputs it consumed*. There is no run id, so two
   partial re-runs cannot be told apart in `results/`.
3. **Library modules resolve paths from `__file__`.** `expr.SEQ_MEMBERS`,
   `expr.COUNTS_CACHE`, and `run_structure`'s `SEQ_DIR`/`OUT_DIR`/`WORK` are
   module-level constants. Paths belong to the CLI; libraries should take them.
4. **Caches and outputs are mixed.** `results/structure_module/pdb/` is both.
5. **Cache validity is presence or a weak key.** References by presence; the
   pseudobulk cache by column names only; AFDB by filename. No checksums.
6. **The family is not declared in one place, despite the README's claim.**
   `core.FOCAL_GENES` and `core.GENE_SYMBOLS` are single-source, but pointing at
   another family additionally requires editing the `_Lba` reference lookup
   (2.4b), the `"Lb-Lb"` class literals (4.2), `TISSUE_COLS` (4.2),
   `SAMPLES` and `MARKERS` (module 3), the HMM accession and the pocket template.
7. **Two environments, no lockfile.** `environment.yml` is unpinned for
   biopython/pandas/numpy/scipy, has a two-stage pip solve for `tmtools`, and
   the ESM2 env is described only in a comment. Neither is captured as a lock.
8. **No tests.** The one defect class the repo has already been bitten by — a
   silent join loss — is guarded at runtime in one function and by nothing else.
   The canonical-pair contract, the label round-trip
   (`label_of`/`gene_from_label`), and the 21-row invariant are all testable in
   milliseconds.
9. **Two orphan outputs** (3.6) and one terminal branch no downstream step reads
   (1.3b), neither marked as such.
10. **`.DS_Store` is tracked** (3 files, 2 modified in the working tree).

## What stays as it is

- The gated (multiplicative) score, the two normalisation scopes, the
  NA-with-a-reason convention, `load_pairs` raising on row loss, the canonical
  pair ordering in `core`, and the deliberate exclusion of `spearman_profile`
  and `duplication_mode` from the score. All decided with evidence in-repo.
- The module/runner split that keeps Modal out of the library code.
