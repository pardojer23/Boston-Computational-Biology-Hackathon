#!/usr/bin/env python3
"""Paralog redundancy features for a gene family, after De Kegel & Ryan (2019).

Implements the paralog-redundancy quantification of:

    De Kegel B, Ryan CJ (2019) Paralog buffering contributes to the variable
    essentiality of genes in cancer cell lines. PLoS Genet 15(10):e1008466.
    doi:10.1371/journal.pgen.1008466

Method, as given in that paper's Methods ("Paralog data"; "Whole genome vs.
small-scale duplicates"):

  1. Paralog relationships come from Ensembl Compara.
  2. Protein sequence identity is taken in BOTH directions for every pair
     (percent of A matched in B, and of B matched in A; these differ because
     the proteins differ in length).
  3. Pairs are kept only if they share >= MIN_IDENTITY percent in both
     directions and both genes are protein-coding.
  4. Each gene is summarised by (i) the number of paralog pairs it is in and
     (ii) the maximum percent of its OWN sequence matched in any paralog.
  5. Pairs are labelled WGD or SSD; a gene is WGD if it is in any WGD pair.

Step 5 in the original consumes curated human ohnolog lists. Those lists are
built by synteny comparison, so for species without one this script reproduces
the criterion directly: a cross-chromosome pair is WGD only if the two loci sit
in a collinear block of reciprocal paralogues AND Compara's inferred duplication
node for that pair is one of the labels passed to --wgd-node.

The node test is exact string matching against the labels you supply, not a
taxonomic-depth comparison -- Compara node labels carry no ordering in the REST
response, so the script cannot tell by itself which of two labels is older. Pass
every node label you consider WGD-era (they are printed in the duplication_node
column of paralog_pairs.csv); any pair at a label you did not list stays SSD.

Works for any species in Ensembl (vertebrates or Ensembl Genomes divisions)
reachable through the unified REST endpoint.

Usage
-----
    python paralog_redundancy.py --species glycine_max \
        --genes GLYMA_10G198800 GLYMA_10G198900 GLYMA_10G199000 \
                GLYMA_10G199100 GLYMA_20G191200 \
        --wgd-node Glycine_subgen._Soja --out results/

All HTTP responses are cached under --cache (default .ensembl_cache), so
re-runs and parameter sweeps cost no network traffic.
"""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REST = "https://rest.ensembl.org"


# --------------------------------------------------------------------------- #
# HTTP                                                                        #
# --------------------------------------------------------------------------- #
def fetch(endpoint: str, cache: Path, key: str, fmt: str = "json", tries: int = 4) -> object:
    """GET an Ensembl REST endpoint, caching the body on disk.

    curl is used rather than requests/urllib because some sandboxes route
    outbound traffic through a proxy the Python HTTP stack cannot reach.
    --http1.1 avoids an HTTP/2 header error seen on EBI-hosted endpoints, and
    the retry loop covers the intermittent 'SSL_read: unexpected eof' this
    endpoint returns under sustained querying.
    """
    ext = "json" if fmt == "json" else "fa"
    path = cache / f"{key}.{ext}"
    if path.exists() and path.stat().st_size > 0:
        try:
            return json.loads(path.read_text()) if fmt == "json" else path.read_text()
        except json.JSONDecodeError:
            path.unlink()

    header = "Content-type:application/json" if fmt == "json" else "Content-type:text/x-fasta"
    cache.mkdir(parents=True, exist_ok=True)
    last = None
    for attempt in range(tries):
        proc = subprocess.run(
            ["curl", "-sSL", "--http1.1", "--retry", "2", "--max-time", "120",
             "-H", header, f"{REST}{endpoint}"],
            capture_output=True, text=True,
        )
        body = proc.stdout
        if proc.returncode == 0 and body and not body.lstrip().startswith("<"):
            if fmt != "json":
                path.write_text(body)
                return body
            try:
                parsed = json.loads(body)
            except json.JSONDecodeError as exc:
                last = exc
            else:
                path.write_text(body)
                return parsed
        else:
            last = proc.stderr.strip() or "HTML or empty body"
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"failed to fetch {endpoint}: {last}")


def lookup_gene(gene: str, cache: Path) -> dict:
    return fetch(f"/lookup/id/{gene}", cache / "lookup", gene)


def paralogues(species: str, gene: str, cache: Path) -> list[dict]:
    d = fetch(f"/homology/id/{species}/{gene}?type=paralogues;sequence=none;cigar_line=0",
              cache / "homology", gene)
    return d["data"][0]["homologies"] if d.get("data") else []


def protein_length(gene: str, cache: Path) -> int:
    fa = fetch(f"/sequence/id/{gene}?type=protein;multiple_sequences=1", cache / "protein",
               gene, fmt="fasta")
    return len("".join(l.strip() for l in fa.splitlines() if not l.startswith(">")))


def region_genes(species: str, chrom: str, start: int, end: int, cache: Path) -> list[dict]:
    start = max(1, start)
    return fetch(
        f"/overlap/region/{species}/{chrom}:{start}-{end}?feature=gene;biotype=protein_coding",
        cache / "region", f"{chrom}_{start}_{end}")


# --------------------------------------------------------------------------- #
# Steps 1-3: family, bidirectional identity, filter                           #
# --------------------------------------------------------------------------- #
def build_family(species: str, seeds: list[str], cache: Path, expand: bool = True,
                 max_family: int = 200) -> tuple[dict, dict]:
    """Return (directed identity records, gene metadata) for the paralog family.

    Compara is queried for each seed; if expand, every newly seen partner is
    queried too, so pairs among non-seed members are not missed. The result is
    the connected component of the seeds in the within-species paralogy graph.
    """
    directed: dict[tuple[str, str], dict] = {}
    seen: set[str] = set()
    queue = list(dict.fromkeys(seeds))
    while queue:
        gene = queue.pop(0)
        if gene in seen:
            continue
        seen.add(gene)
        for h in paralogues(species, gene, cache):
            src, tgt = h["source"], h["target"]
            if tgt.get("species") != src.get("species"):
                continue  # paralogues only, never orthologues
            directed[(src["id"], tgt["id"])] = {
                "perc_id_source": src["perc_id"],
                "perc_id_target": tgt["perc_id"],
                "duplication_node": h.get("taxonomy_level"),
                "homology_type": h.get("type"),
            }
            if expand and tgt["id"] not in seen and len(seen) < max_family:
                queue.append(tgt["id"])
    meta = {}
    for gene in sorted({g for pair in directed for g in pair} | seen):
        d = lookup_gene(gene, cache)
        meta[gene] = {"chrom": str(d["seq_region_name"]), "start": d["start"], "end": d["end"],
                      "strand": d["strand"], "biotype": d["biotype"],
                      "display_name": d.get("display_name") or gene}
    return directed, meta


def pair_table(directed: dict, meta: dict, min_identity: float) -> pd.DataFrame:
    """Undirected pairs carrying identity in both directions, plus the filter."""
    genes = sorted({g for pair in directed for g in pair})
    rows = []
    for a, b in itertools.combinations(genes, 2):
        fwd, rev = directed.get((a, b)), directed.get((b, a))
        if fwd is None and rev is None:
            continue
        rec = fwd if fwd is not None else rev
        # perc_id belongs to whichever gene was the 'source' of the record used
        id_a = rec["perc_id_source"] if fwd is not None else rec["perc_id_target"]
        id_b = rec["perc_id_target"] if fwd is not None else rec["perc_id_source"]
        same_chrom = meta[a]["chrom"] == meta[b]["chrom"]
        rows.append({
            "gene_A": a, "gene_B": b,
            "symbol_A": meta[a]["display_name"], "symbol_B": meta[b]["display_name"],
            "perc_id_A_in_B": id_a, "perc_id_B_in_A": id_b,
            "min_perc_id": min(id_a, id_b),
            "duplication_node": rec["duplication_node"],
            "chrom_A": meta[a]["chrom"], "chrom_B": meta[b]["chrom"],
            "intergenic_kb": (abs(meta[a]["start"] - meta[b]["start"]) / 1e3
                              if same_chrom else np.nan),
            "passes_filter": (min(id_a, id_b) >= min_identity
                              and meta[a]["biotype"] == "protein_coding"
                              and meta[b]["biotype"] == "protein_coding"),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Step 5: duplication mode from synteny                                       #
# --------------------------------------------------------------------------- #
def synteny_block(species: str, locus_a: tuple[str, int], locus_b: tuple[str, int],
                  window: int, cache: Path, min_anchors: int = 5,
                  min_rho: float = 0.5) -> dict:
    """Test whether two loci sit in a collinear block of reciprocal paralogues.

    Genes are collected in a +/-window interval around each locus; paralogues of
    every gene in the second window are fetched and those landing in the first
    window are 'anchors'. A block is called when there are >= min_anchors and
    |Spearman rho| of anchor order >= min_rho (the sign only reports whether the
    blocks are in the same or inverted orientation).
    """
    (chrom_a, pos_a), (chrom_b, pos_b) = locus_a, locus_b
    win_a = {g["id"]: g for g in region_genes(species, chrom_a, pos_a - window,
                                              pos_a + window, cache)}
    win_b = {g["id"]: g for g in region_genes(species, chrom_b, pos_b - window,
                                              pos_b + window, cache)}
    rows = []
    for gene in win_b:
        for h in paralogues(species, gene, cache):
            if h["target"]["id"] in win_a:
                rows.append({"gene_window_B": gene, "gene_window_A": h["target"]["id"],
                             "perc_id_B": h["source"]["perc_id"],
                             "perc_id_A": h["target"]["perc_id"],
                             "duplication_node": h.get("taxonomy_level")})
    anchors = pd.DataFrame(rows)
    out = {"chrom_A": chrom_a, "chrom_B": chrom_b, "window_bp": window,
           "genes_window_A": len(win_a), "genes_window_B": len(win_b),
           "n_anchors": len(anchors), "anchored_genes_B": 0,
           "spearman_rho": np.nan, "is_block": False, "anchors": anchors}
    if anchors.empty:
        return out
    ord_a = {g: i for i, g in enumerate(sorted(win_a, key=lambda x: win_a[x]["start"]))}
    ord_b = {g: i for i, g in enumerate(sorted(win_b, key=lambda x: win_b[x]["start"]))}
    anchors["rank_A"] = anchors.gene_window_A.map(ord_a)
    anchors["rank_B"] = anchors.gene_window_B.map(ord_b)
    anchors["pos_A_mb"] = anchors.gene_window_A.map(lambda g: win_a[g]["start"] / 1e6)
    anchors["pos_B_mb"] = anchors.gene_window_B.map(lambda g: win_b[g]["start"] / 1e6)
    rho = anchors[["rank_A", "rank_B"]].corr(method="spearman").iloc[0, 1]
    out.update(anchored_genes_B=int(anchors.gene_window_B.nunique()),
               spearman_rho=float(rho),
               is_block=bool(len(anchors) >= min_anchors and abs(rho) >= min_rho),
               anchors=anchors)
    return out


def assign_duplication_mode(pairs: pd.DataFrame, meta: dict, species: str, cache: Path,
                            window: int, wgd_nodes: list[str] | None, run_synteny: bool) -> dict:
    """Label each pair WGD or SSD. Same-chromosome pairs are SSD by definition.

    Two independent conditions are materialised as their own columns so a reader
    of the CSV can see why any given pair was called:

      in_homeologous_block       - positional: the pair spans a collinear block
      duplication_node_concordant - temporal: Compara's duplication node for the
                                    pair is one of wgd_nodes

    dup_mode is WGD only where both are true. A pair can therefore sit in a
    homeologous block and still be SSD, which is the correct call for an
    out-paralog whose duplication predates the WGD.
    """
    pairs["in_homeologous_block"] = False
    pairs["arrangement"] = ""
    blocks: dict[tuple[str, str], dict] = {}
    for i, r in pairs.iterrows():
        if r.chrom_A == r.chrom_B:
            pairs.at[i, "arrangement"] = "same chromosome (tandem / local duplicate)"
            continue
        key = tuple(sorted((r.chrom_A, r.chrom_B)))
        if run_synteny and key not in blocks:
            blocks[key] = synteny_block(
                species, (r.chrom_A, meta[r.gene_A]["start"]),
                (r.chrom_B, meta[r.gene_B]["start"]), window, cache)
        blk = blocks.get(key)
        in_block = bool(blk and blk["is_block"])
        pairs.at[i, "in_homeologous_block"] = in_block
        pairs.at[i, "arrangement"] = (f"homeologous block (chr{key[0]}/chr{key[1]})" if in_block
                                      else "different chromosomes, no block support")
    pairs["duplication_node_concordant"] = (pairs.duplication_node.isin(wgd_nodes)
                                            if wgd_nodes else True)
    pairs["dup_mode"] = np.where(pairs.in_homeologous_block
                                 & pairs.duplication_node_concordant, "WGD", "SSD")
    return blocks


# --------------------------------------------------------------------------- #
# Step 4: per-gene summary                                                    #
# --------------------------------------------------------------------------- #
def gene_table(pairs: pd.DataFrame, meta: dict, cache: Path) -> pd.DataFrame:
    kept = pairs[pairs.passes_filter]
    rows = []
    for gene in sorted({*kept.gene_A, *kept.gene_B}):
        sub = kept[(kept.gene_A == gene) | (kept.gene_B == gene)]
        own = [(r.perc_id_A_in_B if r.gene_A == gene else r.perc_id_B_in_A,
                r.gene_B if r.gene_A == gene else r.gene_A) for r in sub.itertuples()]
        best = max(own)
        rows.append({
            "gene": gene, "symbol": meta[gene]["display_name"],
            "chrom": meta[gene]["chrom"], "start": meta[gene]["start"],
            "end": meta[gene]["end"], "strand": meta[gene]["strand"],
            "protein_len_aa": protein_length(gene, cache),
            "n_paralog_pairs": len(sub),
            "max_perc_id_own_seq": round(best[0], 4),
            "closest_paralog": meta[best[1]]["display_name"],
            "mean_perc_id_own_seq": round(float(np.mean([o[0] for o in own])), 4),
            "min_perc_id_own_seq": round(float(min(o[0] for o in own)), 4),
            "dup_mode": "WGD" if (sub.dup_mode == "WGD").any() else "SSD",
            "n_wgd_pairs": int((sub.dup_mode == "WGD").sum()),
            "n_ssd_pairs": int((sub.dup_mode == "SSD").sum()),
        })
    return (pd.DataFrame(rows)
            .sort_values("max_perc_id_own_seq", ascending=False)
            .reset_index(drop=True))


# --------------------------------------------------------------------------- #
# Figure                                                                      #
# --------------------------------------------------------------------------- #
def make_figure(pairs: pd.DataFrame, genes: pd.DataFrame, blocks: dict, out: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 6, "ytick.labelsize": 6,
                         "axes.spines.top": False, "axes.spines.right": False})
    WGD_C, SSD_C, GREY = "#1f6fb4", "#d1832a", "#8a8a8a"
    block = next((b for b in blocks.values() if b["is_block"]), None)
    ncols = 3 if block is not None else 2
    fig = plt.figure(figsize=(3.75 * ncols, 3.5))
    gs = fig.add_gridspec(1, ncols, width_ratios=[1.05, 1.0, 1.15][-ncols:], wspace=0.42)
    axes = [fig.add_subplot(gs[0, i]) for i in range(ncols)]
    letters = iter("abc")

    if block is not None:
        ax, a = axes.pop(0), block["anchors"]
        fam = set(genes.gene)
        focal = a.gene_window_A.isin(fam) | a.gene_window_B.isin(fam)
        ax.scatter(a.pos_A_mb[~focal], a.pos_B_mb[~focal], s=14, c=GREY, lw=0, alpha=.75,
                   label="other anchor pairs")
        ax.scatter(a.pos_A_mb[focal], a.pos_B_mb[focal], s=34, c=WGD_C, lw=0, zorder=3,
                   label="family anchors")
        ax.set_xlabel(f"position on chromosome {block['chrom_A']} (Mb)")
        ax.set_ylabel(f"position on chromosome {block['chrom_B']} (Mb)")
        ax.set_title("Collinear paralogue anchors mark a\nhomeologous block", loc="left")
        ax.legend(frameon=False, loc="upper right", handletextpad=.3, borderpad=.1)
        ax.margins(0.05)
        ax.annotate(f"$\\rho$ = {block['spearman_rho']:.2f}\n"
                    f"{block['anchored_genes_B']}/{block['genes_window_B']} genes anchored",
                    xy=(.03, .06), xycoords="axes fraction", fontsize=6, color=GREY)
        ax.text(-0.20, 1.02, next(letters), transform=ax.transAxes, fontweight="bold",
                fontsize=10, va="bottom", ha="right")

    # --- identity matrix, ordered by chromosome then coordinate
    ax = axes.pop(0)
    order = list(genes.sort_values(["chrom", "start"]).symbol)
    idx = {s: i for i, s in enumerate(order)}
    n = len(order)
    M = np.full((n, n), np.nan)
    wgd_cells = []
    for r in pairs[pairs.passes_filter].itertuples():
        if r.symbol_A not in idx or r.symbol_B not in idx:
            continue
        ia, ib = idx[r.symbol_A], idx[r.symbol_B]
        M[ia, ib], M[ib, ia] = r.perc_id_A_in_B, r.perc_id_B_in_A
        if r.dup_mode == "WGD":
            wgd_cells += [(ia, ib), (ib, ia)]
    lo = max(0.0, float(np.floor(np.nanmin(M) / 5) * 5 - 5))
    im = ax.imshow(M, cmap="viridis", vmin=lo, vmax=100)
    mid = lo + 0.7 * (100 - lo)
    for i in range(n):
        ax.add_patch(Rectangle((i - .5, i - .5), 1, 1, facecolor="#eeeeee", edgecolor="none"))
        for j in range(n):
            if not np.isnan(M[i, j]) and n <= 14:
                ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center", fontsize=6,
                        color="white" if M[i, j] < mid else "black")
    for (i, j) in wgd_cells:
        ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False, ec=WGD_C, lw=1.6))
    ax.set_xticks(range(n), order, rotation=90 if n > 8 else 0)
    ax.set_yticks(range(n), order)
    ax.set_xlabel("matched in this paralogue\nblue outline = WGD pair")
    ax.set_ylabel("% of this gene's sequence")
    ax.set_title("Identity is asymmetric: each gene is\nscored in both directions", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.04)
    cb.set_label("% identity", fontsize=6)
    cb.ax.tick_params(labelsize=6)
    ax.text(-0.20, 1.02, next(letters), transform=ax.transAxes, fontweight="bold",
            fontsize=10, va="bottom", ha="right")

    # --- per-gene redundancy
    ax = axes.pop(0)
    gp = genes.sort_values("max_perc_id_own_seq")
    y = np.arange(len(gp))
    ax.hlines(y, 0, gp.max_perc_id_own_seq, color=GREY, lw=.9, zorder=1)
    ax.scatter(gp.max_perc_id_own_seq, y, s=60, zorder=3,
               c=[WGD_C if m == "WGD" else SSD_C for m in gp.dup_mode])
    for yi, r in zip(y, gp.itertuples()):
        ax.annotate(f"{r.max_perc_id_own_seq:.1f}%  (vs {r.closest_paralog})",
                    xy=(r.max_perc_id_own_seq + 2.5, yi), va="center", fontsize=6)
    ax.set_yticks(y, gp.symbol)
    ax.set_xlim(0, 132)
    ax.set_ylim(-.75, len(gp) - .3)
    ax.set_xlabel("max % of gene's own sequence matched in any paralogue")
    ax.set_title("Redundancy available to each family\nmember", loc="left")
    for c, lab in ((WGD_C, "gene in a WGD pair"), (SSD_C, "gene in SSD pairs only")):
        ax.scatter([], [], c=c, s=40, label=lab)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.62, 0.42),
              handletextpad=.3, borderpad=.1)
    ax.text(-0.30, 1.02, next(letters), transform=ax.transAxes, fontweight="bold",
            fontsize=10, va="bottom", ha="right")

    path = out / "paralog_redundancy.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    return path


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", required=True, help="Ensembl species name, e.g. glycine_max")
    p.add_argument("--genes", required=True, nargs="+", help="seed gene IDs (Ensembl)")
    p.add_argument("--out", default="results", type=Path)
    p.add_argument("--cache", default=".ensembl_cache", type=Path)
    p.add_argument("--min-identity", type=float, default=20.0,
                   help="inclusion floor: percent identity required in BOTH directions "
                        "(paper: 20)")
    p.add_argument("--window", type=int, default=400_000,
                   help="half-width of the synteny window, bp")
    p.add_argument("--wgd-node", default=None, nargs="+", metavar="NODE",
                   help="one or more Compara duplication-node labels considered WGD-era, e.g. "
                        "Glycine_subgen._Soja. Matched exactly; a pair at any other node stays "
                        "SSD even if it spans a homeologous block. Omitted = synteny alone "
                        "decides. Node labels appear in the duplication_node output column.")
    p.add_argument("--no-synteny", action="store_true",
                   help="skip the synteny test (every pair becomes SSD); much faster")
    p.add_argument("--no-expand", action="store_true",
                   help="query Compara for the seed genes only, not their partners")
    p.add_argument("--no-figure", action="store_true")
    args = p.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    print(f"[1/5] family from {len(args.genes)} seed gene(s) in {args.species}", flush=True)
    directed, meta = build_family(args.species, args.genes, args.cache, expand=not args.no_expand)
    print(f"      {len(meta)} genes, {len(directed)} directed homology records", flush=True)

    print("[2/5] bidirectional identity + inclusion filter", flush=True)
    pairs = pair_table(directed, meta, args.min_identity)
    if pairs.empty:
        sys.exit("no paralog pairs found for these seeds")
    print(f"      {len(pairs)} pairs, {int(pairs.passes_filter.sum())} pass "
          f">= {args.min_identity}% in both directions", flush=True)

    print("[3/5] duplication mode" + ("" if not args.no_synteny else " (synteny skipped)"),
          flush=True)
    blocks = assign_duplication_mode(pairs, meta, args.species, args.cache, args.window,
                                     args.wgd_node, run_synteny=not args.no_synteny)
    print(f"      {int((pairs.dup_mode == 'WGD').sum())} WGD / "
          f"{int((pairs.dup_mode == 'SSD').sum())} SSD pairs", flush=True)

    print("[4/5] per-gene summary", flush=True)
    genes = gene_table(pairs, meta, args.cache)

    cols = ["gene_A", "gene_B", "symbol_A", "symbol_B", "perc_id_A_in_B", "perc_id_B_in_A",
            "min_perc_id", "passes_filter", "chrom_A", "chrom_B", "intergenic_kb",
            "in_homeologous_block", "duplication_node", "duplication_node_concordant",
            "dup_mode", "arrangement"]
    pairs[cols].to_csv(args.out / "paralog_pairs.csv", index=False)
    genes.to_csv(args.out / "gene_redundancy.csv", index=False)
    for key, blk in blocks.items():
        if blk["is_block"]:
            blk["anchors"].to_csv(args.out / f"synteny_anchors_chr{key[0]}_chr{key[1]}.csv",
                                  index=False)
    (args.out / "run_metadata.json").write_text(json.dumps({
        "method_source": "De Kegel B, Ryan CJ (2019) PLoS Genet 15:e1008466",
        "species": args.species, "seed_genes": args.genes,
        "min_identity": args.min_identity, "window_bp": args.window,
        "wgd_node": args.wgd_node, "synteny": not args.no_synteny,
        "ensembl_rest": REST,
        "assembly": lookup_gene(args.genes[0], args.cache).get("assembly_name"),
        "blocks": {f"chr{k[0]}_chr{k[1]}": {m: v for m, v in b.items() if m != "anchors"}
                   for k, b in blocks.items()},
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }, indent=2, default=str))

    print("[5/5] figure" if not args.no_figure else "[5/5] figure skipped", flush=True)
    if not args.no_figure:
        make_figure(pairs, genes, blocks, args.out)

    print(f"\nwrote {args.out}/")
    print(genes[["symbol", "chrom", "protein_len_aa", "n_paralog_pairs",
                 "max_perc_id_own_seq", "closest_paralog", "dup_mode"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
