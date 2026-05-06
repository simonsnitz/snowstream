"""Stage: pick a representative protein for each cluster.

Output: families/<key>/representatives.json — JSON list of records, one per
cluster, of the form:

  {
    "uniprot_id": "P43506",         # the chosen representative
    "centroid_uniprot_id": "P43506", # MMseqs2's chosen centroid (may differ)
    "cluster": { "size": 12, "members": ["P43506", "..."] },
    "evidence": {
      "source": "groovdb" | "paperblast" | "default",
      "url": "https://www.groov.bio/entry/TetR/P43506",   # groovdb
      "match": { "subject_uniprot_id": "...",
                 "identity": 0.94, "coverage": 0.98 },     # paperblast
      "papers": [ { "doi": "10.x", "year": 2018 } ]        # paperblast
    }
  }

Selection priority:
  1. Any cluster member exists in groovDB → that member is the rep, link to
     groovDB. If multiple members are in groovDB, prefer the one whose
     sequence is closest to the centroid (highest identity).
  2. Else, run PaperBLAST on the centroid; among hits ≥50% identity, look for
     any whose subject UniProt ID is *also in our cluster*. If found, that
     member is the rep, with the matched papers attached.
  3. Else, fall back to the MMseqs2 centroid; no literature.

This module is pure-Python (the heavy lifting — clustering, FASTA parsing,
remote queries — happens upstream). It's intentionally easy to unit-test.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Callable, Optional

from . import groovdb as groovdb_mod
from . import paperblast as paperblast_mod

log = logging.getLogger(__name__)


PaperBlastQueryFn = Callable[[str], list[paperblast_mod.PaperBlastHit]]


def select_representative(
    centroid: str,
    members: list[str],
    sequences: dict[str, str],
    groovdb_dump: dict[str, dict[str, Any]],
    groovdb_entry_prefix: str,
    paperblast_query: Optional[PaperBlastQueryFn],
    min_paperblast_identity_pct: float = 50.0,
) -> dict:
    """Choose the representative for one cluster and build its evidence record."""
    # 1. groovDB — any cluster member?
    members_in_groov = [m for m in members if m in groovdb_dump]
    if members_in_groov:
        chosen = members_in_groov[0]  # any one is fine; the user said it's OK
        return {
            "uniprot_id": chosen,
            "centroid_uniprot_id": centroid,
            "cluster": {"size": len(members), "members": members},
            "evidence": {
                "source": "groovdb",
                "url": groovdb_mod.make_url(groovdb_entry_prefix, chosen),
                "uniprot_id": chosen,
            },
        }

    # 2. PaperBLAST on the centroid sequence.
    if paperblast_query is not None:
        centroid_seq = sequences.get(centroid, "")
        if centroid_seq:
            hits = paperblast_query(centroid_seq)
            cluster_members = set(members)
            for hit in hits:
                if hit.identity_pct < min_paperblast_identity_pct:
                    continue
                if hit.subject_uniprot_id in cluster_members:
                    return {
                        "uniprot_id": hit.subject_uniprot_id,
                        "centroid_uniprot_id": centroid,
                        "cluster": {"size": len(members), "members": members},
                        "evidence": {
                            "source": "paperblast",
                            "match": {
                                "subject_uniprot_id": hit.subject_uniprot_id,
                                "identity": round(hit.identity_pct / 100.0, 4),
                                "coverage": round(hit.coverage_pct / 100.0, 4),
                            },
                            "papers": [
                                {"doi": p.doi, "year": p.year} for p in hit.papers
                            ],
                        },
                    }

    # 3. Fallback: MMseqs2 centroid.
    return {
        "uniprot_id": centroid,
        "centroid_uniprot_id": centroid,
        "cluster": {"size": len(members), "members": members},
        "evidence": {"source": "default"},
    }


def build_representatives(
    clusters: dict[str, list[str]],
    sequences: dict[str, str],
    groovdb_dump: dict[str, dict[str, Any]],
    groovdb_entry_prefix: str,
    paperblast_cache_dir: Optional[Path],
    use_paperblast: bool = True,
    min_paperblast_identity_pct: float = 50.0,
    on_progress: Optional[Callable[[int, int, str], None]] = None,
) -> list[dict]:
    """Iterate clusters, pick a representative per cluster, return the list."""
    pblast_fn: Optional[PaperBlastQueryFn] = None
    if use_paperblast and paperblast_cache_dir is not None:
        def pblast_fn(seq: str) -> list[paperblast_mod.PaperBlastHit]:
            return paperblast_mod.query(seq, paperblast_cache_dir)

    out: list[dict] = []
    items = sorted(clusters.items())
    for i, (centroid, members) in enumerate(items):
        if on_progress:
            on_progress(i, len(items), centroid)
        rec = select_representative(
            centroid=centroid,
            members=members,
            sequences=sequences,
            groovdb_dump=groovdb_dump,
            groovdb_entry_prefix=groovdb_entry_prefix,
            paperblast_query=pblast_fn,
            min_paperblast_identity_pct=min_paperblast_identity_pct,
        )
        out.append(rec)
    return out


def write_representatives(records: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(records, indent=2) + "\n")
    tmp.replace(path)


def read_representatives(path: Path) -> list[dict]:
    return json.loads(path.read_text())
