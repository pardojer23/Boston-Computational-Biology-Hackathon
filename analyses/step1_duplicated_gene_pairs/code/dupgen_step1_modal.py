"""Step 1 of Li, Zhang & Schmitz 2025 (Plant Cell, 10.1093/plcell/koaf279):
identify soybean duplicated gene pairs with DupGen_finder, common bean as outgroup.

Runs on Modal. DIAMOND is pinned to v2.1.8, the version named in the authors'
step1 notes, so BLAST hit sets are not subject to version drift.

  modal run dupgen_step1_modal.py
"""

import modal

app = modal.App("soybean-dupgen-step1")

DIAMOND_URL = (
    "https://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz"
)
SOYBASE = "https://data.soybase.org"
GLYMA_DIR = f"{SOYBASE}/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"
PHAVU_DIR = f"{SOYBASE}/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d"

GLYMA_PEP = "glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz"
GLYMA_BED = "glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.bed.gz"
PHAVU_PEP = "phavu.G19833.gnm2.ann1.PB8d.protein_primary.faa.gz"
PHAVU_BED = "phavu.G19833.gnm2.ann1.PB8d.gene_models_main.bed.gz"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("build-essential", "perl", "git", "wget", "ca-certificates")
    .run_commands(
        f"wget -q {DIAMOND_URL} -O /tmp/diamond.tar.gz",
        "tar -xzf /tmp/diamond.tar.gz -C /tmp && install -m 755 /tmp/diamond /usr/local/bin/diamond",
        "git clone --depth 1 https://github.com/qiao-xin/DupGen_finder.git /opt/DupGen_finder",
        "cd /opt/DupGen_finder && make && chmod 775 *.pl",
    )
    .env({"PATH": "/opt/DupGen_finder:/usr/local/bin:/usr/bin:/bin"})
)


@app.function(image=image, cpu=16.0, memory=32768, timeout=7200)
def run_step1():
    import gzip
    import os
    import re
    import subprocess
    import shutil

    os.makedirs("/work/raw", exist_ok=True)
    os.chdir("/work")
    log = []

    def sh(cmd, **kw):
        log.append(f"$ {cmd}")
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
        if r.stdout:
            log.append(r.stdout[-3000:])
        if r.returncode != 0:
            log.append(f"[stderr] {r.stderr[-3000:]}")
            raise RuntimeError(f"failed ({r.returncode}): {cmd}\n{r.stderr[-2000:]}")
        return r.stdout

    log.append(sh("diamond --version"))
    log.append("DupGen_finder: " + str(os.path.isfile("/opt/DupGen_finder/DupGen_finder.pl")))

    # ---- inputs -------------------------------------------------------------
    for d, f in ((GLYMA_DIR, GLYMA_PEP), (GLYMA_DIR, GLYMA_BED),
                 (PHAVU_DIR, PHAVU_PEP), (PHAVU_DIR, PHAVU_BED)):
        sh(f"wget -q {d}/{f} -O raw/{f}")

    def read_pep(path):
        """Return {id: sequence}; ids are the first whitespace-delimited token."""
        seqs, cur = {}, None
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    cur = line[1:].split()[0]
                    seqs[cur] = []
                else:
                    seqs[cur].append(line.strip())
        return {k: "".join(v) for k, v in seqs.items()}

    glyma_pep = read_pep(f"raw/{GLYMA_PEP}")
    phavu_pep = read_pep(f"raw/{PHAVU_PEP}")
    log.append(f"glyma primary proteins: {len(glyma_pep)}")
    log.append(f"phavu primary proteins: {len(phavu_pep)}")

    def write_fasta(seqs, path):
        with open(path, "w") as fh:
            for k, v in seqs.items():
                fh.write(f">{k}\n")
                for i in range(0, len(v), 60):
                    fh.write(v[i:i + 60] + "\n")

    write_fasta(glyma_pep, "glyma.pep")
    write_fasta(phavu_pep, "phavu.pep")

    def bed_to_gff(bed_path, keep_ids):
        """DupGen_finder gene-position format: chr <tab> gene <tab> start <tab> end.
        One row per primary-transcript gene, keeping the BED's own chromosome labels
        (which is what reproduces the published Location strings)."""
        rows, seen = [], set()
        with gzip.open(bed_path, "rt") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                chrom, start, end, tx = p[0], p[1], p[2], p[3]
                if tx in keep_ids and tx not in seen:
                    seen.add(tx)
                    rows.append((chrom, tx, start, end))
        return rows

    glyma_rows = bed_to_gff(f"raw/{GLYMA_BED}", set(glyma_pep))
    phavu_rows = bed_to_gff(f"raw/{PHAVU_BED}", set(phavu_pep))
    log.append(f"glyma gff rows: {len(glyma_rows)} | phavu gff rows: {len(phavu_rows)}")

    def write_gff(rows, path):
        with open(path, "w") as fh:
            for r in rows:
                fh.write("\t".join(r) + "\n")

    write_gff(glyma_rows, "glyma.gff")
    write_gff(phavu_rows, "phavu.gff")
    write_gff(glyma_rows + phavu_rows, "glyma_phavu.gff")

    # ---- DIAMOND ------------------------------------------------------------
    PARAMS = "-p 16 --sensitive --max-target-seqs 5 --evalue 1e-10 --quiet"
    sh("diamond makedb --in glyma.pep -d glyma --quiet")
    sh("diamond makedb --in phavu.pep -d phavu --quiet")
    sh(f"diamond blastp -d glyma -q glyma.pep -o glyma.blast {PARAMS}")
    sh(f"diamond blastp -d phavu -q glyma.pep -o glyma_phavu.blast {PARAMS}")
    for f in ("glyma.blast", "glyma_phavu.blast"):
        log.append(f"{f}: {sum(1 for _ in open(f))} hits")

    # ---- filter: identity > 40, coverage > 0.7 ------------------------------
    # coverage = alignment length (col 4) / query protein length, per the
    # authors' 07032025_blast_filter.Rmd.
    lengths = {k: len(v) for k, v in glyma_pep.items()}

    def filter_blast(src, dst, numeric=True):
        kept = 0
        with open(src) as fin, open(dst, "w") as fout:
            for line in fin:
                p = line.rstrip("\n").split("\t")
                if len(p) < 12:
                    continue
                qlen = lengths.get(p[0])
                if not qlen:
                    continue
                cov = int(p[3]) / qlen
                if cov <= 0.7:
                    continue
                # numeric=True is the filter the Methods text describes.
                # numeric=False reproduces R's read.delim2(dec=",") behaviour,
                # where pident stays character and ">" is a string comparison.
                ok = (float(p[2]) > 40) if numeric else (p[2] > "40")
                if ok:
                    fout.write("\t".join(p[:12]) + "\n")
                    kept += 1
        return kept

    variants = {}
    for label, numeric in (("numeric", True), ("lexicographic", False)):
        d = f"/work/data_{label}"
        os.makedirs(d, exist_ok=True)
        n1 = filter_blast("glyma.blast", f"{d}/glyma.blast", numeric)
        n2 = filter_blast("glyma_phavu.blast", f"{d}/glyma_phavu.blast", numeric)
        for g in ("glyma.gff", "glyma_phavu.gff"):
            shutil.copy(g, f"{d}/{g}")
        log.append(f"filter[{label}]: glyma {n1} hits, glyma_phavu {n2} hits")
        variants[label] = d

    # ---- DupGen_finder ------------------------------------------------------
    results = {}
    for label, d in variants.items():
        out = f"/work/results_{label}"
        sh(f"DupGen_finder.pl -i {d} -t glyma -c phavu -o {out}")
        files = sorted(os.listdir(out))
        log.append(f"DupGen_finder[{label}] produced: {files}")
        for fn in files:
            path = os.path.join(out, fn)
            if os.path.isfile(path):
                results[f"{label}/{fn}"] = open(path, "rb").read()

    results["run.log"] = "\n".join(str(x) for x in log).encode()
    return results


@app.local_entrypoint()
def main():
    import os

    out = run_step1.remote()
    base = "modal_out"
    for name, data in out.items():
        path = os.path.join(base, name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "wb") as fh:
            fh.write(data)
        print(f"{len(data):>12,}  {path}")
