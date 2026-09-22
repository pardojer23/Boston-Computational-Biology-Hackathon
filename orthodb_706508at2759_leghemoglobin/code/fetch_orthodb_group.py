"""
Fetch protein sequences for an OrthoDB orthologous group and clean headers
into Newick-safe FASTA labels.

Usage:
    python fetch_orthodb_group.py 706508at2759 group_706508at2759_clean.fasta
"""
import sys
import re
import json
import urllib.request

ORTHODB_BASE = "https://data.orthodb.org/current"


def fetch(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode()


def fetch_group_info(group_id: str) -> dict:
    return json.loads(fetch(f"{ORTHODB_BASE}/group?id={group_id}"))


def fetch_group_fasta(group_id: str) -> str:
    return fetch(f"{ORTHODB_BASE}/fasta?id={group_id}")


def parse_fasta(fasta_text: str) -> list[dict]:
    records = []
    for block in fasta_text.strip().split(">")[1:]:
        lines = block.split("\n")
        header = lines[0]
        seq = "".join(lines[1:])
        m = re.match(r"(\S+)\s+(\{.*\})", header)
        odb_id, meta_str = m.group(1), m.group(2)
        meta = json.loads(meta_str)
        records.append({
            "odb_id": odb_id,
            "organism": meta.get("organism_name"),
            "gene_id": meta.get("pub_gene_id"),
            "description": meta.get("description"),
            "seq": seq,
        })
    return records


def make_label(rec: dict) -> str:
    org = rec["organism"].replace(" ", "_").replace(".", "")
    oid = rec["odb_id"].replace(":", "_")
    return f"{org}__{oid}"


def write_clean_fasta(records: list[dict], out_path: str) -> None:
    labels = [make_label(r) for r in records]
    assert len(set(labels)) == len(labels), "duplicate labels produced"
    with open(out_path, "w") as f:
        for r, lab in zip(records, labels):
            f.write(f">{lab}\n{r['seq']}\n")


def main() -> None:
    group_id = sys.argv[1]
    out_path = sys.argv[2]

    info = fetch_group_info(group_id)
    print(json.dumps(info.get("data", {}), indent=2)[:1000])

    fasta_text = fetch_group_fasta(group_id)
    records = parse_fasta(fasta_text)
    print(f"Fetched {len(records)} sequences for group {group_id} "
          f"across {len(set(r['organism'] for r in records))} organisms")

    write_clean_fasta(records, out_path)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
