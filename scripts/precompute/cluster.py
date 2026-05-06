"""Stage: cluster member sequences at a configurable identity threshold using MMseqs2.

Output: families/<key>/clusters.tsv
        Two columns: representative_id<TAB>member_id (MMseqs2's easy-cluster TSV format).
        For a singleton cluster the representative_id == member_id.
Resume: skip if output exists.

We invoke `mmseqs easy-cluster` and copy its `*_cluster.tsv` artifact to our
expected path, then clean up the temp work directory it produces.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

log = logging.getLogger(__name__)


def cluster(
    fasta_path: Path,
    output_path: Path,
    identity: float = 0.5,
    coverage: float = 0.8,
    force: bool = False,
    mmseqs_bin: str = "mmseqs",
) -> int:
    """Cluster `fasta_path` at the given identity threshold; write the cluster TSV.

    Returns the number of clusters (unique representatives).
    """
    if output_path.exists() and not force:
        reps = {line.split("\t", 1)[0] for line in output_path.read_text().splitlines() if line}
        log.info("clusters.tsv already exists (%d clusters); skipping", len(reps))
        return len(reps)
    if not fasta_path.exists():
        raise FileNotFoundError(f"missing prerequisite: {fasta_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_p = Path(tmpdir)
        prefix = tmpdir_p / "cluster"
        scratch = tmpdir_p / "scratch"
        cmd = [
            mmseqs_bin,
            "easy-cluster",
            str(fasta_path),
            str(prefix),
            str(scratch),
            "--min-seq-id",
            str(identity),
            "-c",
            str(coverage),
            "--cov-mode",
            "0",
        ]
        log.info("running: %s", " ".join(cmd))
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        cluster_tsv = prefix.with_name(prefix.name + "_cluster.tsv")
        if not cluster_tsv.exists():
            raise RuntimeError(f"mmseqs did not produce {cluster_tsv}")
        shutil.copyfile(cluster_tsv, output_path)

    reps = {line.split("\t", 1)[0] for line in output_path.read_text().splitlines() if line}
    log.info("wrote %d clusters to %s", len(reps), output_path)
    return len(reps)


def parse_clusters(clusters_tsv: Path) -> dict[str, list[str]]:
    """Parse the MMseqs2 TSV into {representative_id: [member_ids]}."""
    clusters: dict[str, list[str]] = {}
    for line in clusters_tsv.read_text().splitlines():
        if not line:
            continue
        rep, member = line.split("\t", 1)
        # Strip UniProt fasta-header prefixes if MMseqs included them
        # (e.g. "sp|P43506|...|" → "P43506").
        rep_id = _normalise_id(rep)
        member_id = _normalise_id(member)
        clusters.setdefault(rep_id, []).append(member_id)
    return clusters


def _normalise_id(raw: str) -> str:
    """Strip UniProt 'sp|P12345|NAME' style headers down to the accession."""
    if "|" in raw:
        parts = raw.split("|")
        if len(parts) >= 2 and parts[0] in ("sp", "tr"):
            return parts[1]
    return raw
