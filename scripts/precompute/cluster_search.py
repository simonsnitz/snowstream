"""Stage: compute per-rep homologs from MMseqs2 cluster membership.

Replaces the BLAST step for precompute. For each representative, we run
`mmseqs easy-search` of the rep against the full members.fasta, then keep
only hits whose target is in the rep's own cluster (since cluster members
are by definition the relevant homologs). Output is shaped like the
BLAST homolog list so downstream stages can consume it unchanged.

Output: families/<key>/cluster_homologs.json
  { rep_uniprot_id: [{"uniprot_id", "identity", "coverage"}, ...] }
  Sorted by identity (descending), capped at max_homologs, self-hit excluded.
Resume: skip if output exists.
"""

from __future__ import annotations

import json
import logging
import subprocess
import tempfile
from pathlib import Path

log = logging.getLogger(__name__)


def search_cluster_homologs(
    representatives_fasta: Path,
    members_fasta: Path,
    representatives: list[dict],
    output_path: Path,
    max_homologs: int = 100,
    sensitivity: float = 7.5,
    force: bool = False,
    mmseqs_bin: str = "mmseqs",
) -> dict[str, list[dict]]:
    if output_path.exists() and not force:
        log.info("cluster_homologs.json already exists; skipping")
        return json.loads(output_path.read_text())

    allowed: dict[str, set[str]] = {}
    for rep in representatives:
        rep_id = rep["uniprot_id"]
        members = set(rep["cluster"]["members"])
        members.discard(rep_id)
        allowed[rep_id] = members

    results: dict[str, list[dict]] = {rep_id: [] for rep_id in allowed}

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_p = Path(tmpdir)
        out_m8 = tmpdir_p / "rep_vs_members.m8"
        scratch = tmpdir_p / "scratch"
        cmd = [
            mmseqs_bin,
            "easy-search",
            str(representatives_fasta),
            str(members_fasta),
            str(out_m8),
            str(scratch),
            "-s",
            str(sensitivity),
            "--format-output",
            "query,target,fident,qcov",
            "--max-seqs",
            str(max(max_homologs * 4, 400)),
        ]
        log.info("running: %s", " ".join(cmd))
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        with out_m8.open() as fh:
            for line in fh:
                parts = line.rstrip().split("\t")
                if len(parts) < 4:
                    continue
                q = _normalise_id(parts[0])
                t = _normalise_id(parts[1])
                if q not in allowed or t not in allowed[q]:
                    continue
                try:
                    fident = float(parts[2])
                    qcov = float(parts[3])
                except ValueError:
                    continue
                results[q].append(
                    {
                        "uniprot_id": t,
                        "identity": fident * 100.0,
                        "coverage": qcov * 100.0,
                    }
                )

    for rep_id, hits in results.items():
        best: dict[str, dict] = {}
        for h in hits:
            cur = best.get(h["uniprot_id"])
            if cur is None or h["identity"] > cur["identity"]:
                best[h["uniprot_id"]] = h
        ordered = sorted(best.values(), key=lambda h: h["identity"], reverse=True)
        results[rep_id] = ordered[:max_homologs]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(results, indent=2))
    total_hits = sum(len(v) for v in results.values())
    log.info(
        "wrote cluster homologs for %d reps (%d total hits) to %s",
        len(results),
        total_hits,
        output_path,
    )
    return results


def _normalise_id(raw: str) -> str:
    if "|" in raw:
        parts = raw.split("|")
        if len(parts) >= 2 and parts[0] in ("sp", "tr"):
            return parts[1]
    return raw
