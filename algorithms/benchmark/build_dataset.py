"""Build a merged TetR-family operator benchmarking dataset from two sources.

Inputs (paths configurable via CLI; default to ~/Downloads/):
  * groovDB_dataset.json — list of TF "sensors", each with uniprotID, accession
    (NCBI), alias, operators[]={sequence, doi}. groovDB labels family as 'TETR'.
  * Snowprint_benchmarking_dataset.xlsx — one row per protein, with NCBI
    accession, alias, Known operator sequence (mixed-case: uppercase core,
    lowercase flanking ~5-10 bp), Reference DOI. Family column is 'TetR'.

Output: a JSON file shaped like algorithms/benchmark/datasets/tetr_137.json
  [
    {
      "ncbi_accession": "WP_001224188.1",
      "uniprot_id": null | "A0A...",
      "alias": "GbaA",
      "operators": [
        {"sequence": "...", "doi": "10.1016/...", "source": "snowprint"},
        ...
      ],
      "sources": ["snowprint"] | ["groovDB"] | ["snowprint","groovDB"]
    },
    ...
  ]

Merge key is the NCBI accession (Snowprint's primary; groovDB has it too).
Where the same accession appears in both inputs, operators are concatenated
and both sources are listed. The UniProt id is taken from groovDB when
available — Snowprint's xlsx has no UniProt column.

Usage:
  # Regenerate the default tetr_137.json from ~/Downloads/ inputs:
  python -m algorithms.benchmark.build_dataset

  # Custom dataset:
  python -m algorithms.benchmark.build_dataset \\
      --groovdb path/to/groovDB.json \\
      --snowprint path/to/snowprint.xlsx \\
      --output algorithms/benchmark/datasets/my_set.json
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

DEFAULT_GROOVDB = Path.home() / "Downloads/groovDB_dataset.json"
DEFAULT_SNOWPRINT = Path.home() / "Downloads/Snowprint_benchmarking_dataset.xlsx"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "datasets" / "tetr_137.json"

_DOI_RE = re.compile(r"^(?:https?://(?:dx\.)?doi\.org/)?(10\.\d{4,9}/[^\s]+)$")


def _normalise_doi(doi: str | None) -> str | None:
    if not doi or not isinstance(doi, str):
        return None
    doi = doi.strip()
    m = _DOI_RE.match(doi)
    return m.group(1) if m else doi


def _clean_seq(seq: str | None) -> str | None:
    """Keep only ACGT (any case). Preserves mixed case so downstream code can
    distinguish Snowprint's lowercase-flanking convention if useful."""
    if not seq or not isinstance(seq, str):
        return None
    return re.sub(r"[^ACGTacgt]", "", seq) or None


def load_groovdb(path: Path) -> list[dict]:
    raw = json.loads(path.read_text())
    sensors = raw["sensors"]
    out: list[dict] = []
    for s in sensors:
        if (s.get("family") or "").upper() != "TETR":
            continue
        ops = []
        for o in s.get("operators") or []:
            seq = _clean_seq(o.get("sequence"))
            if not seq:
                continue
            ops.append({
                "sequence": seq,
                "doi": _normalise_doi(o.get("doi")),
                "source": "groovDB",
            })
        if not ops:
            continue
        out.append({
            "ncbi_accession": s.get("accession"),
            "uniprot_id": s.get("uniprotID"),
            "alias": s.get("alias"),
            "operators": ops,
        })
    return out


def load_snowprint(path: Path) -> list[dict]:
    df = pd.read_excel(path)
    df = df[df["Structural Family"] == "TetR"]
    out: list[dict] = []
    for _, row in df.iterrows():
        seq = _clean_seq(row.get("Known operator sequence"))
        if not seq:
            continue
        out.append({
            "ncbi_accession": row.get("Protein NCBI Accession ID"),
            "uniprot_id": None,
            "alias": row.get("Alias"),
            "operators": [{
                "sequence": seq,
                "doi": _normalise_doi(row.get("Reference DOI")),
                "source": "snowprint",
            }],
        })
    return out


def merge(groov: list[dict], snow: list[dict]) -> list[dict]:
    """Merge by NCBI accession; concatenate operators when a protein appears
    in both. When both sources contribute the same operator sequence (case-
    insensitive), keep one and tag its `source` as `["groovDB","snowprint"]`.
    Sources at the protein level is sorted ['groovDB','snowprint']."""
    by_acc: dict[str, dict] = {}
    for entry in groov + snow:
        acc = entry["ncbi_accession"]
        if not acc:
            continue
        if acc not in by_acc:
            by_acc[acc] = {
                "ncbi_accession": acc,
                "uniprot_id": entry.get("uniprot_id"),
                "alias": entry.get("alias"),
                "operators": list(entry["operators"]),
                "sources": set(),
            }
        else:
            existing = by_acc[acc]
            if not existing["uniprot_id"] and entry.get("uniprot_id"):
                existing["uniprot_id"] = entry["uniprot_id"]
            if not existing["alias"] and entry.get("alias"):
                existing["alias"] = entry["alias"]
            existing["operators"].extend(entry["operators"])
        for op in entry["operators"]:
            by_acc[acc]["sources"].add(op["source"])

    # Dedupe operators within each protein by uppercase sequence; merge sources.
    for e in by_acc.values():
        by_seq: dict[str, dict] = {}
        for op in e["operators"]:
            key = op["sequence"].upper()
            if key in by_seq:
                existing_src = by_seq[key]["source"]
                if isinstance(existing_src, str):
                    by_seq[key]["source"] = sorted({existing_src, op["source"]})
                else:
                    by_seq[key]["source"] = sorted(set(existing_src) | {op["source"]})
                if not by_seq[key].get("doi") and op.get("doi"):
                    by_seq[key]["doi"] = op["doi"]
            else:
                by_seq[key] = {**op}
        e["operators"] = list(by_seq.values())

    out = []
    for _, e in by_acc.items():
        e["sources"] = sorted(e["sources"])
        out.append(e)
    out.sort(key=lambda x: x["ncbi_accession"])
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--groovdb", type=Path, default=DEFAULT_GROOVDB,
                    help=f"Path to groovDB JSON export (default: {DEFAULT_GROOVDB})")
    p.add_argument("--snowprint", type=Path, default=DEFAULT_SNOWPRINT,
                    help=f"Path to Snowprint benchmarking xlsx (default: {DEFAULT_SNOWPRINT})")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                    help=f"Output JSON path (default: {DEFAULT_OUTPUT})")
    args = p.parse_args()

    if not args.groovdb.exists():
        raise SystemExit(f"groovDB not found: {args.groovdb}")
    if not args.snowprint.exists():
        raise SystemExit(f"Snowprint xlsx not found: {args.snowprint}")

    groov = load_groovdb(args.groovdb)
    snow = load_snowprint(args.snowprint)
    merged = merge(groov, snow)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(merged, indent=2))

    n = len(merged)
    in_both = sum(1 for e in merged if len(e["sources"]) > 1)
    only_groov = sum(1 for e in merged if e["sources"] == ["groovDB"])
    only_snow = sum(1 for e in merged if e["sources"] == ["snowprint"])
    total_ops = sum(len(e["operators"]) for e in merged)
    multi_op = sum(1 for e in merged if len(e["operators"]) > 1)

    print(f"groovDB TetR proteins:    {len(groov)}")
    print(f"Snowprint TetR proteins:  {len(snow)}")
    print(f"Merged unique proteins:   {n}")
    print(f"  in both sources:        {in_both}")
    print(f"  only groovDB:           {only_groov}")
    print(f"  only snowprint:         {only_snow}")
    print(f"Total operators:          {total_ops}")
    print(f"  proteins with >1 op:    {multi_op}")
    print(f"\nSaved to {args.output}")


if __name__ == "__main__":
    main()
