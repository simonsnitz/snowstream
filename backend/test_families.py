"""Tests for backend/families.py — manifest loading + JSONL byte-offset lookup.

Run with: .venv/bin/python -m unittest backend.test_families
"""

import json
import tempfile
import time
import unittest
from pathlib import Path

from .families import FamilyRegistry, _build_index


def _write_manifest(fam_dir: Path, key: str = "tetr") -> None:
    fam_dir.mkdir(parents=True, exist_ok=True)
    (fam_dir / "family.json").write_text(
        json.dumps(
            {
                "key": key,
                "name": f"{key.upper()} family",
                "interpro_id": "IPR001647",
                "groovdb_entry_prefix": f"https://www.groov.bio/entry/{key}/",
                "cluster_identity_threshold": 0.5,
                "smart_lookup": {"min_identity_pct": 50, "min_coverage_pct": 95},
            }
        )
    )


def _write_jsonl(path: Path, records: list[dict]) -> None:
    path.write_text("".join(json.dumps(r) + "\n" for r in records))


class FamilyRegistryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name) / "families"
        self.fam_dir = self.root / "tetr"
        _write_manifest(self.fam_dir)

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_loads_manifest(self) -> None:
        reg = FamilyRegistry(self.root)
        fams = reg.list()
        self.assertEqual([f.key for f in fams], ["tetr"])
        self.assertEqual(fams[0].interpro_id, "IPR001647")

    def test_skips_underscore_prefixed_dirs(self) -> None:
        # _resources/ holds shared downloads (e.g. groovdb dump) and isn't a family.
        (self.root / "_resources").mkdir()
        (self.root / "_resources" / "groovdb.json").write_text("{}")
        reg = FamilyRegistry(self.root)
        self.assertEqual([f.key for f in reg.list()], ["tetr"])

    def test_lookup_missing_family(self) -> None:
        reg = FamilyRegistry(self.root)
        self.assertIsNone(reg.lookup("nope", "P43506"))

    def test_lookup_missing_uniprot(self) -> None:
        _write_jsonl(self.fam_dir / "predictions.jsonl", [{"uniprot_id": "P43506"}])
        reg = FamilyRegistry(self.root)
        self.assertIsNone(reg.lookup("tetr", "Q9KS52"))

    def test_lookup_returns_full_record(self) -> None:
        records = [
            {"uniprot_id": "P43506", "homologs": [{"uniprot_id": "x"}]},
            {"uniprot_id": "Q9KS52", "consensus_seq": "ACGT", "extra": list(range(50))},
            {"uniprot_id": "P0A6F5", "papers": [{"doi": "10.x", "year": 2018}]},
        ]
        _write_jsonl(self.fam_dir / "predictions.jsonl", records)
        reg = FamilyRegistry(self.root)
        self.assertEqual(reg.lookup("tetr", "P43506"), records[0])
        # Important: a lookup hits the right line even when offsets differ.
        self.assertEqual(reg.lookup("tetr", "Q9KS52"), records[1])
        self.assertEqual(reg.lookup("tetr", "P0A6F5"), records[2])

    def test_index_rebuilds_on_jsonl_mtime_change(self) -> None:
        # Initial state: one record.
        path = self.fam_dir / "predictions.jsonl"
        _write_jsonl(path, [{"uniprot_id": "P43506"}])
        reg = FamilyRegistry(self.root)
        self.assertIsNotNone(reg.lookup("tetr", "P43506"))
        self.assertIsNone(reg.lookup("tetr", "Q9KS52"))

        # Append a new record (precompute would do this).
        time.sleep(0.01)  # ensure mtime ticks on coarse-resolution filesystems
        with path.open("a") as fh:
            fh.write(json.dumps({"uniprot_id": "Q9KS52"}) + "\n")
        # Force an mtime change even when the OS rounds to seconds.
        import os

        st = path.stat()
        os.utime(path, (st.st_atime, st.st_mtime + 1))

        # Lookup should pick up the new record automatically.
        self.assertIsNotNone(reg.lookup("tetr", "Q9KS52"))

    def test_manifest_dict_includes_predictions_count(self) -> None:
        _write_jsonl(
            self.fam_dir / "predictions.jsonl",
            [{"uniprot_id": "P43506"}, {"uniprot_id": "Q9KS52"}],
        )
        reg = FamilyRegistry(self.root)
        manifest = reg.list()[0].manifest_dict()
        self.assertEqual(manifest["predictions_count"], 2)
        # No path internals leaked into the public manifest.
        self.assertNotIn("dir", manifest)
        self.assertNotIn("predictions_path", manifest)

    def test_malformed_lines_are_skipped_in_index(self) -> None:
        path = self.fam_dir / "predictions.jsonl"
        path.write_text(
            json.dumps({"uniprot_id": "P43506"}) + "\n"
            + "not json at all\n"
            + json.dumps({"uniprot_id": "Q9KS52"}) + "\n"
        )
        index = _build_index(path)
        self.assertIn("P43506", index.offsets)
        self.assertIn("Q9KS52", index.offsets)


if __name__ == "__main__":
    unittest.main()
