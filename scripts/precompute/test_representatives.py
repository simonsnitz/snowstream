"""Tests for the representative-selection stage.

These are the most important tests in the precompute suite — the rep
selection logic encodes the priority groovDB > PaperBLAST > centroid and
several edge cases that determine what ends up in predictions.jsonl.

Run: .venv/bin/python -m unittest scripts.precompute.test_representatives
"""

import unittest
from unittest.mock import MagicMock

from .paperblast import Paper, PaperBlastHit
from .representatives import select_representative


class SelectRepresentativeTests(unittest.TestCase):
    sequences = {
        "P43506": "MQRTAA",
        "Q9KS52": "MQRTBB",
        "P0A6F5": "MQRTCC",
        "X11111": "MQRTDD",
    }

    def test_groovdb_member_takes_precedence(self) -> None:
        """If any cluster member is in groovDB, it's the rep regardless of PaperBLAST."""
        groov_dump = {"Q9KS52": {"name": "TetR-X"}}
        # PaperBLAST returns a hit for a different cluster member; should be ignored.
        pblast = MagicMock(
            return_value=[
                PaperBlastHit("P43506", 99.0, 99.0, [Paper(doi="10.x", year=2020)])
            ]
        )

        rec = select_representative(
            centroid="P43506",
            members=["P43506", "Q9KS52", "P0A6F5"],
            sequences=self.sequences,
            groovdb_dump=groov_dump,
            groovdb_entry_prefix="https://www.groov.bio/entry/TetR/",
            paperblast_query=pblast,
        )

        self.assertEqual(rec["uniprot_id"], "Q9KS52")
        self.assertEqual(rec["evidence"]["source"], "groovdb")
        self.assertEqual(rec["evidence"]["url"], "https://www.groov.bio/entry/TetR/Q9KS52")
        # PaperBLAST should not have been queried.
        pblast.assert_not_called()

    def test_paperblast_match_in_cluster(self) -> None:
        """No groovDB; a PaperBLAST hit ≥50 % is in the cluster — that hit is the rep."""
        pblast = MagicMock(
            return_value=[
                PaperBlastHit("UNRELATED", 99.0, 99.0, [Paper("10.unrelated", 2018)]),
                PaperBlastHit("P0A6F5", 75.0, 98.0, [Paper("10.real", 2019)]),
            ]
        )

        rec = select_representative(
            centroid="P43506",
            members=["P43506", "Q9KS52", "P0A6F5"],
            sequences=self.sequences,
            groovdb_dump={},
            groovdb_entry_prefix="https://www.groov.bio/entry/TetR/",
            paperblast_query=pblast,
        )

        self.assertEqual(rec["uniprot_id"], "P0A6F5")
        self.assertEqual(rec["evidence"]["source"], "paperblast")
        self.assertEqual(rec["evidence"]["match"]["subject_uniprot_id"], "P0A6F5")
        self.assertEqual(rec["evidence"]["match"]["identity"], 0.75)
        self.assertEqual(rec["evidence"]["papers"], [{"doi": "10.real", "year": 2019}])

    def test_paperblast_below_threshold_ignored(self) -> None:
        """A PaperBLAST hit at <50 % identity is dropped; falls back to centroid."""
        pblast = MagicMock(
            return_value=[PaperBlastHit("P0A6F5", 40.0, 99.0, [])],  # 40% < 50%
        )

        rec = select_representative(
            centroid="P43506",
            members=["P43506", "P0A6F5"],
            sequences=self.sequences,
            groovdb_dump={},
            groovdb_entry_prefix="https://www.groov.bio/entry/TetR/",
            paperblast_query=pblast,
        )

        self.assertEqual(rec["uniprot_id"], "P43506")
        self.assertEqual(rec["evidence"]["source"], "default")

    def test_paperblast_hit_outside_cluster_ignored(self) -> None:
        """A PaperBLAST hit ≥50 % whose UniProt ID isn't in our cluster is ignored."""
        pblast = MagicMock(
            return_value=[PaperBlastHit("OUTSIDER", 99.0, 99.0, [Paper("10.x", 2020)])]
        )

        rec = select_representative(
            centroid="P43506",
            members=["P43506", "Q9KS52"],
            sequences=self.sequences,
            groovdb_dump={},
            groovdb_entry_prefix="https://www.groov.bio/entry/TetR/",
            paperblast_query=pblast,
        )

        self.assertEqual(rec["uniprot_id"], "P43506")
        self.assertEqual(rec["evidence"]["source"], "default")

    def test_no_paperblast_at_all(self) -> None:
        """When PaperBLAST is disabled (None), fall back to centroid silently."""
        rec = select_representative(
            centroid="P43506",
            members=["P43506", "Q9KS52"],
            sequences=self.sequences,
            groovdb_dump={},
            groovdb_entry_prefix="",
            paperblast_query=None,
        )
        self.assertEqual(rec["uniprot_id"], "P43506")
        self.assertEqual(rec["evidence"]["source"], "default")

    def test_record_includes_full_cluster_membership(self) -> None:
        rec = select_representative(
            centroid="P43506",
            members=["P43506", "Q9KS52", "P0A6F5"],
            sequences=self.sequences,
            groovdb_dump={},
            groovdb_entry_prefix="",
            paperblast_query=None,
        )
        self.assertEqual(rec["cluster"]["size"], 3)
        self.assertEqual(rec["cluster"]["members"], ["P43506", "Q9KS52", "P0A6F5"])
        self.assertEqual(rec["centroid_uniprot_id"], "P43506")


if __name__ == "__main__":
    unittest.main()
