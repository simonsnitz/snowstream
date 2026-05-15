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


def _write_family(
    root: Path,
    key: str = "tetr",
    *,
    with_dmnd: bool = True,
    members: list[dict] | None = None,
    smart_lookup: dict | None = None,
) -> Path:
    """Materialise a fake family directory with all the files smart-lookup expects.

    `members` is a list of per-member records to write into members_predictions.jsonl.
    Each record needs at least a `uniprot_id`; promoter/operon/etc. are passed through
    to the synthesized result record.
    """
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
                "smart_lookup": smart_lookup
                or {"min_identity_pct": 50, "min_coverage_pct": 95, "max_homologs": 30},
            }
        )
    )
    if with_dmnd:
        # smart_lookup only checks `members_with_promoters_dmnd.exists()`; contents
        # are irrelevant when we mock _diamond_hits.
        (fam_dir / "members_with_promoters.dmnd").write_text("fake-dmnd")
    if members is None:
        members = [{"uniprot_id": "P43506", "promoter": "ACGT", "operon": {"genes": []}}]
    with (fam_dir / "members_predictions.jsonl").open("w") as fh:
        for m in members:
            fh.write(json.dumps(m) + "\n")
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
        # Synthesized record carries the member's cached promoter through homologs.
        self.assertEqual(len(hit.record["homologs"]), 1)
        self.assertEqual(hit.record["homologs"][0]["uniprot_id"], "P43506")
        self.assertEqual(hit.record["homologs"][0]["promoter"], "ACGT")
        self.assertEqual(hit.record["homologs"][0]["identity"], 87.0)

    def test_below_identity_threshold_returns_none(self) -> None:
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
            return_value=[{"uniprot_id": "P43506", "identity_pct": 49.0, "coverage_pct": 99.0}],
        ):
            hit = sl.lookup(req, registry=self.registry)
        self.assertIsNone(hit)

    def test_below_coverage_threshold_returns_none(self) -> None:
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
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
        mock_hit.assert_not_called()

    def test_missing_dmnd_skips_family(self) -> None:
        (self.fam_dir / "members_with_promoters.dmnd").unlink()
        registry = FamilyRegistry(self.root)
        req = _make_req()
        with patch.object(sl, "_diamond_hits") as mock_hit:
            hit = sl.lookup(req, registry=registry)
        self.assertIsNone(hit)
        mock_hit.assert_not_called()

    def test_hit_uniprot_not_in_members_jsonl(self) -> None:
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
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
        called_seq = mock_hit.call_args.args[0]
        self.assertEqual(called_seq, "MNNNN")

    def test_first_hit_below_coverage_second_above_returns_second(self) -> None:
        """Regression: diamond's top hit by bit score can fail coverage while
        a lower-ranked hit clears both thresholds. The lower hit should win
        as the anchor (and the only homolog, since the first didn't qualify).
        """
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
            return_value=[
                {"uniprot_id": "OTHER", "identity_pct": 99.0, "coverage_pct": 85.0},
                {"uniprot_id": "P43506", "identity_pct": 66.0, "coverage_pct": 97.0},
            ],
        ):
            hit = sl.lookup(req, registry=self.registry)
        self.assertIsNotNone(hit)
        self.assertEqual(hit.uniprot_id, "P43506")
        self.assertEqual(len(hit.record["homologs"]), 1)

    def test_multiple_qualifying_hits_become_homologs(self) -> None:
        # Family has 3 members with promoters; diamond ranks them.
        fam_dir = _write_family(
            self.root,
            key="multi",
            members=[
                {"uniprot_id": "M1", "promoter": "AAAA", "operon": {"genes": ["g1"]}},
                {"uniprot_id": "M2", "promoter": "CCCC", "operon": {"genes": ["g2"]}},
                {"uniprot_id": "M3", "promoter": "GGGG", "operon": {"genes": ["g3"]}},
            ],
        )
        registry = FamilyRegistry(self.root)
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
            side_effect=[
                # The 'multi' family is alphabetically first; we only need its
                # results. side_effect returns the next list on each call —
                # smart_lookup short-circuits on the first family that yields
                # at least one qualifying hit.
                [
                    {"uniprot_id": "M1", "identity_pct": 99.0, "coverage_pct": 99.0},
                    {"uniprot_id": "M2", "identity_pct": 85.0, "coverage_pct": 96.0},
                    {"uniprot_id": "M3", "identity_pct": 60.0, "coverage_pct": 95.0},
                ],
            ],
        ):
            hit = sl.lookup(req, registry=registry)
        self.assertIsNotNone(hit)
        self.assertEqual(hit.family_key, "multi")
        self.assertEqual(hit.uniprot_id, "M1")  # top hit becomes the anchor
        self.assertEqual([h["uniprot_id"] for h in hit.record["homologs"]], ["M1", "M2", "M3"])
        # Each homolog carries the per-query identity, not anything from the cache.
        self.assertEqual(hit.record["homologs"][1]["identity"], 85.0)
        self.assertEqual(hit.record["homologs"][1]["promoter"], "CCCC")

    def test_max_homologs_caps_synthesized_list(self) -> None:
        fam_dir = _write_family(
            self.root,
            key="bigfam",
            members=[{"uniprot_id": f"M{i}", "promoter": "X"} for i in range(10)],
            smart_lookup={"min_identity_pct": 50, "min_coverage_pct": 95, "max_homologs": 3},
        )
        registry = FamilyRegistry(self.root)
        req = _make_req()
        with patch.object(
            sl,
            "_diamond_hits",
            return_value=[
                {"uniprot_id": f"M{i}", "identity_pct": 90.0, "coverage_pct": 99.0} for i in range(10)
            ],
        ):
            hit = sl.lookup(req, registry=registry)
        self.assertEqual(len(hit.record["homologs"]), 3)

    def test_multiple_families_first_above_threshold_wins(self) -> None:
        # Add a second family alphabetically before tetr, with a hit of its own.
        _write_family(
            self.root,
            key="araC",
            members=[{"uniprot_id": "X11111", "promoter": "TTTT"}],
        )
        registry = FamilyRegistry(self.root)
        req = _make_req()

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
        self.assertEqual(hit.family_key, "araC")


class EventPayloadTests(unittest.TestCase):
    def test_payload_shape(self) -> None:
        hit = sl.SmartLookupHit(
            family_key="tetr",
            uniprot_id="P43506",
            identity_pct=87.5,
            coverage_pct=99.1,
            record={"homologs": []},
        )
        payload = sl.to_event_payload(hit)
        self.assertEqual(payload["matched_via"]["family"], "tetr")
        self.assertEqual(payload["matched_via"]["uniprot_id"], "P43506")
        self.assertEqual(payload["matched_via"]["identity"], 0.875)
        self.assertEqual(payload["matched_via"]["coverage"], 0.991)
        self.assertEqual(payload["result"], {"homologs": []})


if __name__ == "__main__":
    unittest.main()
