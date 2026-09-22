"""
Fetch an OrthoDB orthologous-group's metadata and member protein sequences,
and clean headers into unique, Newick-safe FASTA labels.

Species-agnostic -- works for any OrthoDB group id.

Usage:
    python orthodb_fetch.py <group_id> <output_clean_fasta>
"""
import json
import re
import sys
import urllib.request

ORTHODB_BASE = "https://data.orthodb.org/current"


def fetch(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "research-agent"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode()


def fetch_group_info(group_id: str) -> dict:
    """Group metadata: name/description, level, gene count, species count."""
    return json.loads(fetch(f"{ORTHODB_BASE}/group?id={group_id}"))


def fetch_group_fasta(group_id: str) -> str:
    return fetch(f"{ORTHODB_BASE}/fasta?id={group_id}")


def parse_fasta(fasta_text: str) -> list:
    """Each OrthoDB FASTA header is '<odb_id> <json metadata>'."""
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
    org = (rec["organism"] or "unknown").replace(" ", "_").replace(".", "")
    oid = rec["odb_id"].replace(":", "_")
    return f"{org}__{oid}"


def tree_label_for_odb_id(odb_id: str, organism: str) -> str:
    """Reconstruct the same label from just organism + odb_id (for downstream
    lookups without re-parsing the whole FASTA)."""
    return make_label({"organism": organism, "odb_id": odb_id})


def write_clean_fasta(records: list, out_path: str) -> list:
    labels = [make_label(r) for r in records]
    if len(set(labels)) != len(labels):
        # disambiguate any residual collisions (e.g. duplicate-copy species)
        from collections import Counter
        seen = Counter()
        new_labels = []
        for lab in labels:
            seen[lab] += 1
            new_labels.append(lab if seen[lab] == 1 else f"{lab}_{seen[lab]}")
        labels = new_labels
    with open(out_path, "w") as f:
        for r, lab in zip(records, labels):
            f.write(f">{lab}\n{r['seq']}\n")
    return labels


def fetch_and_clean_group(group_id: str, out_fasta_path: str) -> dict:
    """One-shot: fetch group info + fasta, clean, write, and return a summary
    dict (also used by the pipeline orchestrator)."""
    info = fetch_group_info(group_id)
    fasta_text = fetch_group_fasta(group_id)
    records = parse_fasta(fasta_text)
    labels = write_clean_fasta(records, out_fasta_path)

    n_species = len(set(r["organism"] for r in records))
    return {
        "group_id": group_id,
        "group_name": info.get("data", {}).get("name"),
        "n_genes": len(records),
        "n_species": n_species,
        "records": records,
        "labels": labels,
        "clean_fasta_path": out_fasta_path,
    }


def main():
    group_id, out_path = sys.argv[1], sys.argv[2]
    summary = fetch_and_clean_group(group_id, out_path)
    print(f"Group {group_id}: {summary['group_name']!r}, "
          f"{summary['n_genes']} genes across {summary['n_species']} species")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
