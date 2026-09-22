"""
Map the soybean (Glycine max) leghemoglobin paralogs pulled from OrthoDB group
706508at2759 onto the Glyma.XXGXXXXXX locus IDs used by GSE270392's expression
matrices. These two ID systems do not agree directly: OrthoDB's own
"GLYMA_##G######v4" tag is only populated for some records, and GSE270392 uses
plain "Glyma.##G######" (no version suffix, dot instead of underscore).

Mapping strategy, in the order actually needed to resolve all 5 records:
  1. For records that already carry an OrthoDB "GLYMA_##G######v4" tag,
     reformat directly to "Glyma.##G######".
  2. For records with only a placeholder gene_id (e.g. "N-50", "N-2"), cross-
     reference NCBI's RefSeq annotation of the same genome assembly
     (GCF_000004515.6 = Wm82.a4.v1, matching OrthoDB's "v4" suffix) by
     gene_synonym -- NCBI's synonym list for a locus includes the exact same
     placeholder tag when historically assigned by the classic nodulin
     literature (e.g. LB2's synonyms include "N-50,Nodulin-50").
  3. Confirm/obtain the exact Glyma locus tag for those NCBI-identified genes
     by BLASTP against SoyBase's Wm82.a4.v1 protein database
     (https://sequenceserver.soybase.org), since NCBI's own annotation does
     not carry a Glyma cross-reference.

Usage:
    python map_leghemoglobin_ids.py
"""
import gzip
import json
import re
import urllib.request
import uuid

NCBI_ASSEMBLY_FTP = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/"
    "GCF_000004515.6_Glycine_max_v4.0/"
)
SOYBASE_BLAST_URL = "https://sequenceserver.soybase.org/"
WM82_A4_PROTEIN_DB_ID = "982cbf1cfaa76dd4e7589d7d1b2d60a0"  # Wm82.a4.v1 protein db


def fetch(url, data=None, headers=None):
    req = urllib.request.Request(url, data=data, headers=headers or {"User-Agent": "research-agent"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read()


def convert_orthodb_glyma_tag(gene_id):
    """"GLYMA_10G198800v4" -> "Glyma.10G198800" """
    if not gene_id or not gene_id.startswith("GLYMA_"):
        return None
    core = re.sub(r"v\d+$", "", gene_id[len("GLYMA_"):])
    return f"Glyma.{core}"


def find_ncbi_gene_by_synonym(placeholder_id):
    """Search the NCBI Glycine max v4.0 GFF3 gene_synonym fields for an exact
    match to a classic nodulin placeholder tag (e.g. "N-50")."""
    gff_url = NCBI_ASSEMBLY_FTP + "GCF_000004515.6_Glycine_max_v4.0_genomic.gff.gz"
    gff_text = gzip.decompress(fetch(gff_url)).decode()
    for line in gff_text.splitlines():
        if "\tgene\t" not in line:
            continue
        attrs = line.split("\t")[8]
        m = re.search(r"gene_synonym=([^;]+)", attrs)
        if m and placeholder_id in m.group(1).split(","):
            gene_id_m = re.search(r"GeneID:(\d+)", attrs)
            name_m = re.search(r"Name=([^;]+)", attrs)
            desc_m = re.search(r"description=([^;]+)", attrs)
            return {
                "ncbi_gene_id": gene_id_m.group(1) if gene_id_m else None,
                "symbol": name_m.group(1) if name_m else None,
                "description": desc_m.group(1) if desc_m else None,
            }
    return None


def blastp_soybase(query_fasta, db_id=WM82_A4_PROTEIN_DB_ID):
    """Submit a BLASTP search to SoyBase's public SequenceServer instance and
    return the parsed JSON result (blocks until the job finishes since these
    short protein queries typically complete in well under a minute)."""
    boundary = uuid.uuid4().hex

    def field(name, value):
        return f"--{boundary}\r\nContent-Disposition: form-data; name=\"{name}\"\r\n\r\n{value}\r\n".encode()

    body = b"".join([
        field("sequence", query_fasta),
        field("method", "blastp"),
        field("databases[]", db_id),
    ]) + f"--{boundary}--\r\n".encode()

    req = urllib.request.Request(
        SOYBASE_BLAST_URL, data=body, method="POST",
        headers={"User-Agent": "research-agent",
                 "Content-Type": f"multipart/form-data; boundary={boundary}"},
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        job_url = r.url

    import time
    for _ in range(20):
        raw = fetch(job_url + ".json")
        try:
            result = json.loads(raw)
        except json.JSONDecodeError:
            time.sleep(3)
            continue
        if result.get("queries") is not None:
            return result
        time.sleep(3)
    raise TimeoutError(f"BLAST job at {job_url} did not finish in time")


def best_hit_glyma_id(blast_result, query_id):
    for q in blast_result["queries"]:
        if q["id"].startswith(query_id):
            if not q["hits"]:
                return None, None
            top = q["hits"][0]
            hsp = top["hsps"][0]
            pct_identity = 100.0 * hsp["identity"] / q["length"]
            acc_parts = top["accession"].split(".")
            return acc_parts[0] + "." + acc_parts[1], pct_identity
    return None, None


def main():
    # The 5 Glycine max records returned by the OrthoDB group 706508at2759
    # fetch (see ../../orthodb_706508at2759_leghemoglobin/code/fetch_orthodb_group.py)
    gmax_records = [
        {"odb_id": "3847_0:006070", "gene_id": "GLYMA_10G198900v4", "description": "hypothetical protein",
         "seq": "MGAFTEKQEALVNSSFEAFKANLPHHSVVFFNSILEKAPAAKNMFSFLGDAVDPKNPKLAGHAEKLFGLVRDSAVQLQTKGLVVADATLGPIHTQKGVTDLQFAVVKEALLKTIKEAVGDKWSEELSNAWEVAYDEIAAAIKKAMAIGSLV"},
        {"odb_id": "3847_0:006217", "gene_id": "N-50", "description": "Leghemoglobin C1",
         "seq": "MGAFTEKQEALVSSSFEAFKANIPQYSVVFYNSILEKAPAAKDLFSFLANGVDPTNPKLTGHAEKLFALVRDSAGQLKTNGTVVADAALVSIHAQKAVTDPQFVVVKEALLKTIKEAVGGNWSDELSSAWEVAYDELAAAIKKA"},
        {"odb_id": "3847_0:0065e0", "gene_id": "GLYMA_10G198800v4", "description": "Leghemoglobin C3",
         "seq": "MGAFTDKQEALVSSSFEAFKTNIPQYSVVFYTSILEKAPVAKDLFSFLANGVDPTNPKLTGHAEKLFGLVRDSAGQLKASGTVVIDAALGSIHAQKAITDPQFVVVKEALLKTIKEAVGDKWSDELSSAWEVAYDELAAAIKKAF"},
        {"odb_id": "3847_0:0065fa", "gene_id": "N-2", "description": "Leghemoglobin A",
         "seq": "MVAFTEKQDALVSSSFEAFKANIPQYSVVFYTSILEKAPAAKDLFSFLANGVDPTNPKLTGHAEKLFALVRDSAGQLKASGTVVADAALGSVHAQKAVTDPQFVVVKEALLKTIKAAVGDKWSDELSRAWEVAYDELAAAIKKA"},
        {"odb_id": "3847_0:00ca2d", "gene_id": "GLYMA_20G191200v4", "description": "Leghemoglobin C2",
         "seq": "MGAFTEKQEALVSSSFEAFKANIPQYSVVFYTSILEKAPAAKDLFSFLSNGVDPSNPKLTGHAEKLFGLVRDSAGQLKANGTVVADAALGSIHAQKAITDPQFVVVKEALLKTIKEAVGDKWSDELSSAWEVAYDELAAAIKKAF"},
    ]

    query_fasta = "\n".join(f">{r['odb_id']}|{r['gene_id']}|{r['description']}\n{r['seq']}" for r in gmax_records)
    blast_result = blastp_soybase(query_fasta)

    final_map = []
    for r in gmax_records:
        direct = convert_orthodb_glyma_tag(r["gene_id"])
        ncbi_hit = None if direct else find_ncbi_gene_by_synonym(r["gene_id"])
        glyma_id, pct_identity = best_hit_glyma_id(blast_result, r["odb_id"])
        final_map.append({
            "odb_id": r["odb_id"],
            "orthodb_gene_id": r["gene_id"],
            "orthodb_description": r["description"],
            "ncbi_synonym_hit": ncbi_hit,
            "glyma_id_direct_from_orthodb_tag": direct,
            "glyma_id_from_blast": glyma_id,
            "blast_pct_identity": pct_identity,
        })
        print(json.dumps(final_map[-1], indent=2))

    with open("leghemoglobin_id_map.json", "w") as f:
        json.dump(final_map, f, indent=2)


if __name__ == "__main__":
    main()
