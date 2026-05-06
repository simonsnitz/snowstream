"""Tests for backend/blast_remote.py.

Run: .venv/bin/python -m unittest backend.test_blast_remote
"""

import json
import unittest
from unittest.mock import MagicMock

from . import blast_remote


def _resp(text: str = "", json_payload=None, status: int = 200):
    """Lightweight requests.Response stand-in."""
    r = MagicMock()
    r.status_code = status
    r.ok = 200 <= status < 300
    r.text = text
    if json_payload is not None:
        r.json.return_value = json_payload
    r.raise_for_status = MagicMock()
    return r


class SubmitTests(unittest.TestCase):
    def test_extracts_rid_and_rtoe(self) -> None:
        http = MagicMock()
        http.post.return_value = _resp("blah blah\n    RID = ABC123\n    RTOE = 60\nblah")
        rid, rtoe = blast_remote.submit("MQRTAA", http=http)
        self.assertEqual(rid, "ABC123")
        self.assertEqual(rtoe, 60)
        # Sent the right form data.
        kwargs = http.post.call_args.kwargs
        self.assertEqual(kwargs["data"]["DATABASE"], "nr")
        self.assertEqual(kwargs["data"]["PROGRAM"], "blastp")
        self.assertEqual(kwargs["data"]["QUERY"], "MQRTAA")

    def test_missing_rid_raises(self) -> None:
        http = MagicMock()
        http.post.return_value = _resp("no rid here")
        with self.assertRaises(RuntimeError):
            blast_remote.submit("MQRTAA", http=http)


class PollTests(unittest.TestCase):
    def test_status_codes(self) -> None:
        http = MagicMock()
        http.get.return_value = _resp("    Status=READY\n    ")
        self.assertEqual(blast_remote.poll("RID1", http=http), "READY")
        http.get.return_value = _resp("    Status=WAITING\n    ")
        self.assertEqual(blast_remote.poll("RID1", http=http), "WAITING")
        http.get.return_value = _resp("    Status=FAILED\n    ")
        self.assertEqual(blast_remote.poll("RID1", http=http), "FAILED")
        http.get.return_value = _resp("garbage")
        self.assertEqual(blast_remote.poll("RID1", http=http), "UNKNOWN")


class ParseXmlTests(unittest.TestCase):
    # Minimal BLAST XML that NCBIXML.read accepts. Two hits: one strong
    # (passes default ident/cov), one weak (filtered out).
    SAMPLE_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.13</BlastOutput_version>
  <BlastOutput_reference>...</BlastOutput_reference>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>query</BlastOutput_query-def>
  <BlastOutput_query-len>200</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>1</Parameters_sc-match>
      <Parameters_sc-mismatch>-3</Parameters_sc-mismatch>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>L</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>query</Iteration_query-def>
      <Iteration_query-len>200</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|123|ref|WP_013083972.1|</Hit_id>
          <Hit_def>HTH-type repressor</Hit_def>
          <Hit_accession>WP_013083972.1</Hit_accession>
          <Hit_len>200</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>400</Hsp_bit-score>
              <Hsp_score>1000</Hsp_score>
              <Hsp_evalue>1e-100</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>200</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>200</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>200</Hsp_identity>
              <Hsp_positive>200</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>200</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gi|456|ref|WP_999.1|</Hit_id>
          <Hit_def>poor hit</Hit_def>
          <Hit_accession>WP_999.1</Hit_accession>
          <Hit_len>200</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>20</Hsp_bit-score>
              <Hsp_score>40</Hsp_score>
              <Hsp_evalue>1.0</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>50</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>50</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>10</Hsp_identity>
              <Hsp_positive>10</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>50</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1</Statistics_db-num>
          <Statistics_db-len>200</Statistics_db-len>
          <Statistics_hsp-len>10</Statistics_hsp-len>
          <Statistics_eff-space>1000</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

    def test_filters_by_thresholds(self) -> None:
        df = blast_remote.parse_xml(self.SAMPLE_XML, ident_cutoff=40, cov_cutoff=90)
        # First hit passes (100% identity, 100% coverage). Second is dropped
        # (20% identity, 25% coverage).
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["Uniprot Id"], "WP_013083972.1")
        self.assertEqual(df.iloc[0]["Identity"], 100.0)
        self.assertEqual(df.iloc[0]["Coverage"], 100.0)


class ExtractAccessionTests(unittest.TestCase):
    def test_pulls_refseq_from_pipe_format(self) -> None:
        self.assertEqual(
            blast_remote._extract_accession_from_id("gi|123|ref|WP_X.1|"),
            "WP_X.1",
        )

    def test_falls_back_to_last_segment(self) -> None:
        self.assertEqual(blast_remote._extract_accession_from_id("BARE_ACC"), "BARE_ACC")


class TranslateToUniprotTests(unittest.TestCase):
    def test_maps_refseq_to_uniprot(self) -> None:
        http = MagicMock()
        http.get.return_value = _resp(
            json_payload={
                "results": [
                    {
                        "primaryAccession": "P43506",
                        "uniProtKBCrossReferences": [
                            {"database": "RefSeq", "id": "WP_013083972.1"}
                        ],
                    },
                    {
                        "primaryAccession": "Q9KS52",
                        "uniProtKBCrossReferences": [
                            {"database": "RefSeq", "id": "WP_OTHER.2"}
                        ],
                    },
                ]
            }
        )
        result = blast_remote.translate_to_uniprot(["WP_013083972.1", "WP_OTHER.2", "WP_NOMATCH.1"], http=http)
        self.assertEqual(result["WP_013083972.1"], "P43506")
        self.assertEqual(result["WP_OTHER.2"], "Q9KS52")
        self.assertNotIn("WP_NOMATCH.1", result)

    def test_strips_version_suffix_when_matching(self) -> None:
        http = MagicMock()
        # NCBI returned "WP_X.1" but UniProt's xref says "WP_X" (no version).
        http.get.return_value = _resp(
            json_payload={
                "results": [
                    {
                        "primaryAccession": "P43506",
                        "uniProtKBCrossReferences": [{"database": "RefSeq", "id": "WP_X"}],
                    },
                ]
            }
        )
        result = blast_remote.translate_to_uniprot(["WP_X.1"], http=http)
        self.assertEqual(result["WP_X.1"], "P43506")

    def test_request_failure_returns_empty(self) -> None:
        http = MagicMock()
        http.get.side_effect = blast_remote.requests.ConnectionError("nope")
        result = blast_remote.translate_to_uniprot(["WP_X.1"], http=http)
        self.assertEqual(result, {})


class BlastAgainstNrIntegrationTests(unittest.TestCase):
    def test_end_to_end_with_mocked_http(self) -> None:
        """Submit → poll → fetch → parse → translate, all mocked."""
        http = MagicMock()

        def get_side_effect(*args, **kwargs):
            params = kwargs.get("params", {}) or {}
            cmd = params.get("CMD")
            if cmd == "Get" and params.get("FORMAT_OBJECT") == "SearchInfo":
                return _resp("Status=READY")
            if cmd == "Get" and params.get("FORMAT_TYPE") == "XML":
                return _resp(ParseXmlTests.SAMPLE_XML)
            # UniProt translation:
            return _resp(
                json_payload={
                    "results": [
                        {
                            "primaryAccession": "P43506",
                            "uniProtKBCrossReferences": [
                                {"database": "RefSeq", "id": "WP_013083972.1"}
                            ],
                        }
                    ]
                }
            )

        http.post.return_value = _resp("RID = ABC123\nRTOE = 0\n")
        http.get.side_effect = get_side_effect

        # Patch time.sleep so the test runs instantly.
        with unittest.mock.patch.object(blast_remote.time, "sleep"):
            df = blast_remote.blast_against_nr(
                "MQRTAA",
                {"ident_cutoff": 40, "cov_cutoff": 90},
                http=http,
                poll_interval=0.0,
            )
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["Uniprot Id"], "P43506")
        self.assertEqual(df.iloc[0]["Identity"], 100.0)


import unittest.mock  # noqa: E402  (used in patch.object above)


if __name__ == "__main__":
    unittest.main()
