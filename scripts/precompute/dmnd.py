"""Stage: build representatives.fasta and representatives.dmnd from the
selected representatives. The Diamond DB is what PR 4's smart-lookup will
BLAST user input against.

Inputs:
- members.fasta (everyone)
- representatives.json (selected reps + evidence)

Outputs:
- representatives.fasta (subset of members.fasta containing only the reps)
- representatives.dmnd  (Diamond DB built from the above)
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)


def _parse_fasta(fasta_path: Path) -> dict[str, tuple[str, str]]:
    """Return {uniprot_id: (header_line, sequence)}."""
    out: dict[str, tuple[str, str]] = {}
    current_id: str | None = None
    current_header = ""
    current_chunks: list[str] = []
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            if current_id is not None:
                out[current_id] = (current_header, "".join(current_chunks))
            current_header = line
            current_id = _extract_uniprot(line)
            current_chunks = []
        else:
            current_chunks.append(line.strip())
    if current_id is not None:
        out[current_id] = (current_header, "".join(current_chunks))
    return out


def _extract_uniprot(header: str) -> str:
    """Extract a UniProt accession from a FASTA header line.

    UniProt's headers look like `>sp|P43506|FOO_BAR ...` or
    `>tr|A0A...|...`; we want the middle column. Falls back to the first
    whitespace-separated token (without leading `>`).
    """
    # `>sp|P43506|FOO_BAR …`
    body = header[1:] if header.startswith(">") else header
    if "|" in body:
        parts = body.split("|")
        if len(parts) >= 2 and parts[0] in ("sp", "tr"):
            return parts[1]
    return body.split()[0]


def write_representatives_fasta(
    members_fasta: Path,
    representatives_json: list[dict],
    output_path: Path,
) -> int:
    """Write a FASTA containing only the representatives' sequences."""
    if not members_fasta.exists():
        raise FileNotFoundError(f"missing prerequisite: {members_fasta}")
    parsed = _parse_fasta(members_fasta)
    rep_ids = [rec["uniprot_id"] for rec in representatives_json]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with output_path.open("w") as fh:
        for uid in rep_ids:
            entry = parsed.get(uid)
            if not entry:
                log.warning("no sequence found for representative %s", uid)
                continue
            header, seq = entry
            # Normalise the header to '>UNIPROT_ID' so Diamond emits clean
            # subject IDs in BLAST output.
            fh.write(f">{uid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
            written += 1
    log.info("wrote %d representative sequences to %s", written, output_path)
    return written


def build_dmnd(fasta_path: Path, dmnd_path: Path, force: bool = False, diamond_bin: str = "diamond") -> None:
    """Run `diamond makedb` to produce the representatives Diamond DB."""
    if dmnd_path.exists() and not force:
        log.info("%s already exists; skipping", dmnd_path)
        return
    if not fasta_path.exists():
        raise FileNotFoundError(f"missing prerequisite: {fasta_path}")

    cmd = [diamond_bin, "makedb", "--in", str(fasta_path), "-d", str(dmnd_path.with_suffix(""))]
    log.info("running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if not dmnd_path.exists():
        raise RuntimeError(f"diamond did not produce {dmnd_path}")
    log.info("built Diamond DB at %s", dmnd_path)
