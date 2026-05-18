"""Stage: build a RefSeq → UniProt cross-reference index from the InterPro TSV.

The TSV's first column is a UniProt accession and the second is a semicolon-
separated list of linked RefSeq accessions. Inverting that lets us resolve a
user-supplied RefSeq accession back to a UniProt member of the family, which
is the key into our members.fasta sequence index. With both indexes loaded
at backend startup, RefSeq and UniProt queries skip the NCBI/UniProt round
trips that were dominating per-query latency (~3s of a ~3.4s RefSeq query).

Output: families/<key>/refseq_to_uniprot.json — a single JSON object.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def build_refseq_to_uniprot(tsv_path: Path, output_path: Path, force: bool = False) -> int:
    if output_path.exists() and not force:
        existing = json.loads(output_path.read_text())
        log.info("%s already exists (%d entries); skipping", output_path.name, len(existing))
        return len(existing)

    mapping: dict[str, str] = {}
    with tsv_path.open() as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            uid, refseq_col = parts[0], parts[1].strip(";")
            if not refseq_col:
                continue
            for ref in refseq_col.split(";"):
                ref = ref.strip()
                if ref:
                    # If a RefSeq appears in multiple TSV rows we keep the
                    # first UniProt mapping — RefSeq → UniProt is usually
                    # 1:1 within a single family.
                    mapping.setdefault(ref, uid)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(mapping))
    log.info("wrote %d RefSeq → UniProt mappings to %s", len(mapping), output_path)
    return len(mapping)
