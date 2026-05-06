"""Smaller unit tests for the supporting precompute modules.

Run: .venv/bin/python -m unittest scripts.precompute.test_other_stages
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock

from . import cluster as cluster_mod
from . import dmnd as dmnd_mod
from . import groovdb as groovdb_mod
from . import interpro as interpro_mod
from . import paperblast as paperblast_mod
from . import sequences as sequences_mod


def _mock_response(payload, status_code: int = 200):
    """Build a thing that quacks like requests.Response for our needs."""
    r = MagicMock()
    r.status_code = status_code
    r.ok = 200 <= status_code < 300
    r.json.return_value = payload
    r.text = payload if isinstance(payload, str) else json.dumps(payload)
    r.raise_for_status = MagicMock()
    return r


class InterProTests(unittest.TestCase):
    def test_pagination_collects_all_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "members.txt"
            page1 = {
                "next": "https://example/api?page=2",
                "results": [{"metadata": {"accession": "P43506"}}, {"metadata": {"accession": "Q9KS52"}}],
            }
            page2 = {
                "next": None,
                "results": [{"metadata": {"accession": "P0A6F5"}}],
            }
            http = MagicMock()
            http.get.side_effect = [_mock_response(page1), _mock_response(page2)]
            n = interpro_mod.fetch_members("IPR001647", output, http=http)
            self.assertEqual(n, 3)
            self.assertEqual(output.read_text().splitlines(), ["P43506", "Q9KS52", "P0A6F5"])

    def test_resume_skips_when_output_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "members.txt"
            output.write_text("ALREADY\nHERE\n")
            http = MagicMock()
            n = interpro_mod.fetch_members("IPR001647", output, http=http)
            self.assertEqual(n, 2)
            http.get.assert_not_called()


class SequencesTests(unittest.TestCase):
    def test_writes_concatenated_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            members = Path(tmp) / "members.txt"
            members.write_text("P43506\nQ9KS52\n")
            output = Path(tmp) / "members.fasta"
            http = MagicMock()
            http.get.return_value = _mock_response(">P43506\nMQRT\n>Q9KS52\nMQQQ\n")
            n = sequences_mod.fetch_sequences(members, output, http=http)
            self.assertEqual(n, 2)
            self.assertIn(">P43506", output.read_text())
            self.assertIn(">Q9KS52", output.read_text())


class GroovDBTests(unittest.TestCase):
    def test_load_list_of_records(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "groovdb.json"
            path.write_text(json.dumps([{"uniprot_id": "P43506", "name": "TetR-X"}]))
            dump = groovdb_mod.load(path)
            self.assertIn("P43506", dump)

    def test_load_dict_of_categories(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "groovdb.json"
            path.write_text(json.dumps({"TetR": [{"uniprotID": "Q9KS52"}], "AraC": [{"accession": "X11111"}]}))
            dump = groovdb_mod.load(path)
            self.assertIn("Q9KS52", dump)
            self.assertIn("X11111", dump)

    def test_load_missing_file_returns_empty(self) -> None:
        dump = groovdb_mod.load(Path("/nonexistent/path"))
        self.assertEqual(dump, {})

    def test_make_url(self) -> None:
        url = groovdb_mod.make_url("https://www.groov.bio/entry/TetR/", "P43506")
        self.assertEqual(url, "https://www.groov.bio/entry/TetR/P43506")
        # Trailing-slash-insensitive.
        url2 = groovdb_mod.make_url("https://www.groov.bio/entry/TetR", "P43506")
        self.assertEqual(url2, "https://www.groov.bio/entry/TetR/P43506")


class PaperBlastTests(unittest.TestCase):
    def test_cloudflare_challenge_returns_no_hits(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp)
            http = MagicMock()
            http.get.return_value = _mock_response(
                "<html>Just a moment... please Enable JavaScript</html>"
            )
            hits = paperblast_mod.query("MQRT", cache, http=http)
            self.assertEqual(hits, [])

    def test_request_failure_returns_no_hits(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp)
            http = MagicMock()
            http.get.side_effect = paperblast_mod.requests.ConnectionError("nope")
            hits = paperblast_mod.query("MQRT", cache, http=http)
            self.assertEqual(hits, [])


class ClusterParseTests(unittest.TestCase):
    def test_parse_clusters_strips_uniprot_headers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "clusters.tsv"
            path.write_text("sp|P43506|FOO\tsp|P43506|FOO\nsp|P43506|FOO\ttr|Q9KS52|BAR\n")
            clusters = cluster_mod.parse_clusters(path)
            self.assertEqual(clusters, {"P43506": ["P43506", "Q9KS52"]})


class DmndFastaTests(unittest.TestCase):
    def test_extracts_only_representatives(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            members = Path(tmp) / "members.fasta"
            members.write_text(
                ">sp|P43506|FOO\nMQRTAA\n>sp|Q9KS52|BAR\nMQRTBB\n>sp|P0A6F5|BAZ\nMQRTCC\n"
            )
            output = Path(tmp) / "reps.fasta"
            reps = [{"uniprot_id": "P43506"}, {"uniprot_id": "P0A6F5"}]
            n = dmnd_mod.write_representatives_fasta(members, reps, output)
            self.assertEqual(n, 2)
            text = output.read_text()
            self.assertIn(">P43506\n", text)
            self.assertIn(">P0A6F5\n", text)
            self.assertNotIn(">Q9KS52", text)
            # Headers should be normalised to bare uniprot IDs.
            self.assertNotIn("|FOO", text)


if __name__ == "__main__":
    unittest.main()
