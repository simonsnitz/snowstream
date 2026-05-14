"""Tests for backend/smart_lookup.py.

Run: .venv/bin/python -m unittest backend.test_smart_lookup
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from .families import FamilyRegistry
from .schemas import PredictRequest
from . import smart_lookup as sl


def _write_family(root: Path, key: str = "tetr", *, with_dmnd: bool = True) -> Path:
    fam_dir = root / key
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
    if with_dmnd:
        # The smart_lookup code only checks `dmnd_path.exists()`; the file's
        # contents are irrelevant when we mock _diamond_top_hit.
        (fam_dir / "representatives.dmnd").write_text("fake-dmnd")
    return fam_dir


def _make_req(input_value: str = "MQRTAA", input_method: str = "Protein sequence", **overrides) -> PredictRequest:
    payload = {
        "input_method": input_method,
        "input_value": input_value,
        **overrides,
    }
    return PredictRequest(**payload)


class LookupTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name) / "families"
        self.fam_dir = _write_family(self.root)
        # Tiny JSONL with one cached prediction.
        (self.fam_dir / "predictions.jsonl").write_text(
            json.dumps(
                {
                    "uniprot_id": "P43506",
                    "evidence": {
                        "source": "groovdb",
                        "url": "https://www.groov.bio/entry/TetR/P43506",
                    },
                    "consensus_seq": "ACGT",
                }
            )
            + "\n"
        )
        self.registry = FamilyRegistry(self.root)

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_above_thresholds_returns_hit(self) -> None:
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
            return_value=[{"uniprot_id": "P43506", "identity_pct": 87.0, "coverage_pct": 99.0}],
        ):
            hit = sl.lookup(req, registry=self.registry)
        self.assertIsNotNone(hit)
        self.assertEqual(hit.family_key, "tetr")
        self.assertEqual(hit.uniprot_id, "P43506")
        self.assertEqual(hit.identity_pct, 87.0)
        self.assertEqual(hit.record["consensus_seq"], "ACGT")

    def test_below_identity_threshold_returns_none(self) -> None:
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
            # 49% identity < 50% threshold
            return_value=[{"uniprot_id": "P43506", "identity_pct": 49.0, "coverage_pct": 99.0}],
        ):
            hit = sl.lookup(req, registry=self.registry)
        self.assertIsNone(hit)

    def test_below_coverage_threshold_returns_none(self) -> None:
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
            # 90% coverage < 95% threshold
            return_value=[{"uniprot_id": "P43506", "identity_pct": 87.0, "coverage_pct": 90.0}],
        ):
            hit = sl.lookup(req, registry=self.registry)
        self.assertIsNone(hit)

    def test_no_diamond_hit_returns_none(self) -> None:
        req = _make_req()
        with patch.object(sl, "_diamond_hits", return_value=[]):
            hit = sl.lookup(req, registry=self.registry)
        self.assertIsNone(hit)

    def test_force_skips_smart_lookup(self) -> None:
        req = _make_req(force=True)
        with patch.object(
            sl,
            "_diamond_hits",
            return_value=[{"uniprot_id": "P43506", "identity_pct": 87.0, "coverage_pct": 99.0}],
        ) as mock_hit:
            hit = sl.lookup(req, registry=self.registry)
        self.assertIsNone(hit)
        # The diamond runner should never be called when force is set.
        mock_hit.assert_not_called()

    def test_missing_dmnd_skips_family(self) -> None:
        # Family exists in the manifest, but no representatives.dmnd was built.
        (self.fam_dir / "representatives.dmnd").unlink()
        registry = FamilyRegistry(self.root)
        req = _make_req()
        with patch.object(sl, "_diamond_hits") as mock_hit:
            hit = sl.lookup(req, registry=registry)
        self.assertIsNone(hit)
        mock_hit.assert_not_called()

    def test_hit_uniprot_not_in_predictions_jsonl(self) -> None:
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
            # Diamond claims a great match for an ID that *isn't* in our JSONL.
            return_value=[{"uniprot_id": "ORPHAN", "identity_pct": 99.0, "coverage_pct": 99.0}],
        ):
            hit = sl.lookup(req, registry=self.registry)
        self.assertIsNone(hit)

    def test_protein_sequence_input_is_used_directly(self) -> None:
        req = _make_req(input_method="Protein sequence", input_value="MNNNN")
        with patch.object(
            sl,
            "_diamond_hits",
            return_value=[{"uniprot_id": "P43506", "identity_pct": 99.0, "coverage_pct": 99.0}],
        ) as mock_hit:
            sl.lookup(req, registry=self.registry)
        # The query passed to diamond should be the raw input sequence —
        # no UniProt/RefSeq fetch when input_method is "Protein sequence".
        called_seq = mock_hit.call_args.args[0]
        self.assertEqual(called_seq, "MNNNN")

    def test_first_hit_below_coverage_second_above_returns_second(self) -> None:
        # Regression for the "max-target-seqs 1" bug: diamond's top hit by
        # bit score can be a partial-coverage match that fails the coverage
        # threshold, while a lower-bit-score hit clears both thresholds and
        # is the right cluster to use.
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
            return_value=[
                # Top by score: 99% ID but 85% coverage — fails 95% coverage gate.
                {"uniprot_id": "OTHER", "identity_pct": 99.0, "coverage_pct": 85.0},
                # Lower-ranked but clears both thresholds — should win.
                {"uniprot_id": "P43506", "identity_pct": 66.0, "coverage_pct": 97.0},
            ],
        ):
            hit = sl.lookup(req, registry=self.registry)
        self.assertIsNotNone(hit)
        self.assertEqual(hit.uniprot_id, "P43506")

    def test_multiple_families_first_above_threshold_wins(self) -> None:
        # Add a second family alphabetically before tetr, with a hit.
        other_dir = self.root / "araC"
        _write_family(self.root, "araC")
        (other_dir / "predictions.jsonl").write_text(
            json.dumps({"uniprot_id": "X11111", "consensus_seq": "TTTT"}) + "\n"
        )
        registry = FamilyRegistry(self.root)
        req = _make_req()

        # Both family lookups would hit; we expect the first registered family.
        with patch.object(
            sl,
            "_diamond_hits",
            side_effect=[
                [{"uniprot_id": "X11111", "identity_pct": 99.0, "coverage_pct": 99.0}],
                [{"uniprot_id": "P43506", "identity_pct": 99.0, "coverage_pct": 99.0}],
            ],
        ):
            hit = sl.lookup(req, registry=registry)
        self.assertIsNotNone(hit)
        self.assertEqual(hit.family_key, "araC")  # alphabetically first


class EventPayloadTests(unittest.TestCase):
    def test_payload_shape(self) -> None:
        hit = sl.SmartLookupHit(
            family_key="tetr",
            uniprot_id="P43506",
            identity_pct=87.5,
            coverage_pct=99.0,
            record={"uniprot_id": "P43506", "consensus_seq": "ACGT"},
        )
        payload = sl.to_event_payload(hit)
        self.assertEqual(payload["matched_via"]["family"], "tetr")
        self.assertEqual(payload["matched_via"]["uniprot_id"], "P43506")
        self.assertEqual(payload["matched_via"]["identity"], 0.875)
        self.assertEqual(payload["matched_via"]["coverage"], 0.99)
        self.assertEqual(payload["result"]["consensus_seq"], "ACGT")


if __name__ == "__main__":
    unittest.main()
