"""
Fetch coding-sequence (CDS) nucleotide sequences for the four functional
soybean (Glycine max) leghemoglobin paralogs used in the aBSREL branch-site
selection analysis, and validate each against the corresponding OrthoDB
protein record from group 706508at2759.

The gene->mRNA accession mapping was obtained from the NCBI RefSeq
Glycine_max_v4.0 (GCF_000004515.6) GFF3 annotation (see
../../gse270392_leghemoglobin_redundancy/code/map_leghemoglobin_ids.py for how
each Glyma ID / gene symbol was originally resolved from OrthoDB); it is
hardcoded here since it is a small, fixed lookup already established in that
prior step.

The fifth OrthoDB Glycine max record in this group (odb_id 3847_0:006070,
Glyma.10G198900) is a pseudogene in the current NCBI annotation (pseudo=true,
GeneID 100777546) with no valid CDS feature, and is excluded from this
codon-based selection analysis.

Usage:
    python fetch_soybean_cds.py
"""
import json
import urllib.request
from textwrap import wrap

MRNA_MAP = {
    "LB1": {"paralog": "Lbc3", "glyma_id": "Glyma.10G198800", "mrna": "NM_001248494.2",
            "gene_id": "100527391", "odb_id": "3847_0:0065e0"},
    "LB2": {"paralog": "Lbc1", "glyma_id": "Glyma.10G199000", "mrna": "NM_001358072.1",
            "gene_id": "100785236", "odb_id": "3847_0:006217"},
    "LB3": {"paralog": "Lba", "glyma_id": "Glyma.10G199100", "mrna": "NM_001248999.3",
            "gene_id": "100527427", "odb_id": "3847_0:0065fa"},
    "LB4": {"paralog": "Lbc2", "glyma_id": "Glyma.20G191200", "mrna": "NM_001248319.3",
            "gene_id": "100527379", "odb_id": "3847_0:00ca2d"},
}

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "research-agent"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def translate(seq):
    codons = wrap(seq, 3)
    aa = "".join(CODON_TABLE.get(c, "X") for c in codons if len(c) == 3)
    return aa.rstrip("*")


def parse_fasta_single(text):
    lines = text.strip().split("\n")
    return lines[0], "".join(lines[1:])


def tree_label(odb_id):
    return f"Glycine_max__{odb_id.replace(':', '_')}"


def main():
    # OrthoDB protein sequences for the 4 functional paralogs (group 706508at2759)
    orthodb_protein = {
        "GLYMA_10G198800v4": "MGAFTDKQEALVSSSFEAFKTNIPQYSVVFYTSILEKAPVAKDLFSFLANGVDPTNPKLTGHAEKLFGLVRDSAGQLKASGTVVIDAALGSIHAQKAITDPQFVVVKEALLKTIKEAVGDKWSDELSSAWEVAYDELAAAIKKAF",
        "N-50": "MGAFTEKQEALVSSSFEAFKANIPQYSVVFYNSILEKAPAAKDLFSFLANGVDPTNPKLTGHAEKLFALVRDSAGQLKTNGTVVADAALVSIHAQKAVTDPQFVVVKEALLKTIKEAVGGNWSDELSSAWEVAYDELAAAIKKA",
        "N-2": "MVAFTEKQDALVSSSFEAFKANIPQYSVVFYTSILEKAPAAKDLFSFLANGVDPTNPKLTGHAEKLFALVRDSAGQLKASGTVVADAALGSVHAQKAVTDPQFVVVKEALLKTIKAAVGDKWSDELSRAWEVAYDELAAAIKKA",
        "GLYMA_20G191200v4": "MGAFTEKQEALVSSSFEAFKANIPQYSVVFYTSILEKAPAAKDLFSFLSNGVDPSNPKLTGHAEKLFGLVRDSAGQLKANGTVVADAALGSIHAQKAITDPQFVVVKEALLKTIKEAVGDKWSDELSSAWEVAYDELAAAIKKAF",
    }
    odb_gene_id_lookup = {"LB1": "GLYMA_10G198800v4", "LB2": "N-50", "LB3": "N-2", "LB4": "GLYMA_20G191200v4"}

    records = []
    for key, info in MRNA_MAP.items():
        acc = info["mrna"]
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta_cds_na&retmode=text"
        header, seq = parse_fasta_single(fetch(url).decode())

        translated = translate(seq)
        odb_seq = orthodb_protein[odb_gene_id_lookup[key]]
        match = translated == odb_seq
        n_diff = sum(1 for a, b in zip(translated, odb_seq) if a != b) if len(translated) == len(odb_seq) else None

        print(f"{key} ({info['paralog']}, {info['glyma_id']}, {acc}): "
              f"{len(seq)} nt CDS, translation matches OrthoDB protein: {match}"
              + (f" ({n_diff} aa difference)" if not match and n_diff is not None else ""))

        records.append({**info, "cds_seq": seq, "translation_matches_orthodb": match})

    with open("leghemoglobin_gmax_cds.fasta", "w") as f:
        for r in records:
            f.write(f">{tree_label(r['odb_id'])}\n{r['cds_seq']}\n")

    with open("leghemoglobin_gmax_cds_metadata.json", "w") as f:
        json.dump(records, f, indent=2)

    print("\nExcluded: Glyma.10G198900 (odb_id 3847_0:006070) -- annotated as a "
          "pseudogene by NCBI (GeneID 100777546, pseudo=true), no valid CDS feature.")


if __name__ == "__main__":
    main()
