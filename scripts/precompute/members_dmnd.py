"""Stage: build the family's `members_with_promoters.dmnd` Diamond DB.

Reads `members_predictions.jsonl` to pick the subset of family members that
have a non-null promoter, pulls their sequences out of `members.fasta`, writes
a filtered FASTA, and runs `diamond makedb` on it. The result is the new
smart-lookup index: BLASTing the query against it directly produces a set of
member hits, each of which has cached operon+promoter ready to use.
"""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)


def _members_with_promoters(predictions_path: Path) -> set[str]:
    keep: set[str] = set()
    with predictions_path.open() as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("promoter"):
                uid = r.get("uniprot_id")
                if uid:
                    keep.add(uid)
    return keep


def _normalise_fasta_id(raw: str) -> str:
    """Strip UniProt 'sp|P12345|NAME' style headers down to the accession."""
    if "|" in raw:
        parts = raw.split("|")
        if len(parts) >= 2 and parts[0] in ("sp", "tr"):
            return parts[1]
    return raw.split()[0]


def write_filtered_fasta(
    members_fasta: Path,
    keep_ids: set[str],
    output_fasta: Path,
) -> int:
    """Stream members.fasta and write only the entries whose accession is in
    `keep_ids` to output_fasta. Returns the count written."""
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    keep_current = False
    with members_fasta.open() as fin, output_fasta.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                acc = _normalise_fasta_id(line[1:].strip())
                keep_current = acc in keep_ids
                if keep_current:
                    fout.write(line)
                    written += 1
            elif keep_current:
                fout.write(line)
    log.info("wrote %d filtered sequences to %s", written, output_fasta)
    return written


def build_members_dmnd(
    members_fasta: Path,
    predictions_path: Path,
    output_fasta: Path,
    output_dmnd: Path,
    force: bool = False,
    diamond_bin: str = "diamond",
) -> int:
    """End-to-end builder: filter members.fasta to those with promoters, then
    run `diamond makedb`. Returns the number of sequences in the resulting DB."""
    if output_dmnd.exists() and not force:
        log.info("%s already exists; skipping", output_dmnd)
        return 0

    keep = _members_with_promoters(predictions_path)
    log.info("members with promoter: %d", len(keep))
    if not keep:
        raise RuntimeError("no members with promoters — did the predict stage run?")

    n = write_filtered_fasta(members_fasta, keep, output_fasta)
    if n == 0:
        raise RuntimeError(f"no sequences matched the keep set in {members_fasta}")

    cmd = [
        diamond_bin,
        "makedb",
        "--in",
        str(output_fasta),
        "-d",
        str(output_dmnd.with_suffix("")),
    ]
    log.info("running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.info("built %s with %d sequences", output_dmnd, n)
    return n
