"""operon_fetch v0 — original acc2operon algorithm.

Port of src/accID2operon.py. Given IPG coordinates for a regulator, fetches
a ±10 kb genome fragment, parses every CDS in it, locates the regulator by
matching its start/stop positions in the CDS headers, then walks neighbouring
CDSs to build the operon.

Returns {"operon": [...], "protein_index": int, "genome": str} or None.

Each gene in the operon list is a dict with keys: alias, description,
accession, direction ('+' or '-'), start, stop.
"""

from __future__ import annotations

import re
from typing import Optional

from .._shared import ncbi_efetch_nuccore


# --- helpers (kept private; identical semantics to legacy code) ----------

def _get_genes(genome_id: str, start_pos: int, stop_pos: int) -> tuple[Optional[list[str]], Optional[int]]:
    """Fetch genome fragment around the regulator and locate the CDS header
    that contains both start and stop coordinates."""
    windows = [
        (start_pos - 10000, stop_pos + 10000),
        (start_pos - 5000,  stop_pos + 5000),
        (start_pos,         stop_pos + 5000),
        (start_pos - 5000,  stop_pos),
    ]
    genome_lines = None
    for s, e in windows:
        try:
            url_payload = ncbi_efetch_nuccore(genome_id, s, e, strand=1,
                                              rettype="fasta_cds_aa")
            if url_payload:
                genome_lines = url_payload.split("\n")
                break
        except Exception:
            continue
    if not genome_lines:
        return None, None

    re_start = re.compile(str(start_pos))
    re_stop  = re.compile(str(stop_pos))
    gene_index = 0
    reg_index: Optional[int] = None
    genes: list[str] = []
    for line in genome_lines:
        if not line:
            continue
        if line[0] != ">":
            continue
        if re_start.search(line) and re_stop.search(line):
            reg_index = gene_index
        gene_index += 1
        genes.append(line)
    return genes, reg_index


def _parse_header(fasta_header: str) -> dict:
    """Parse a CDS FASTA header line into a structured gene dict."""
    metadata: dict = {}
    for chunk in fasta_header.split(" ["):
        if chunk[:10] == "locus_tag=":
            metadata["alias"] = chunk[10:-1]
        elif chunk[:8] == "protein=":
            metadata["description"] = chunk[8:-1].replace("'", "")
        elif chunk[:11] == "protein_id=":
            metadata["accession"] = chunk[11:-1]
        elif chunk[:9] == "location=":
            if chunk[9:20] == "complement(":
                metadata["direction"] = "-"
                location = chunk[20:-2]
            else:
                metadata["direction"] = "+"
                location = chunk[9:-1]
            parts = location.split("..")
            metadata["start"] = int(re.sub(r"\D", "", parts[0]))
            metadata["stop"]  = int(re.sub(r"\D", "", parts[1]))
    metadata.setdefault("accession", "")
    return metadata


def _walk_operon(all_genes: list[str], reg_idx: int, seq_start: int,
                 strand: str) -> tuple[list[dict], int]:
    """Walk neighbouring CDSs to build the operon. Same inclusion rules as
    the legacy code: always include immediate neighbours; include further
    neighbours that share direction; stop on divergence or when distance
    exceeds 8 kb."""

    def _walk(geneStrand, direction, nextGene, geneList, index):
        while geneStrand == nextGene["direction"]:
            if direction == "+":
                next_index = index + 1
            elif direction == "-":
                next_index = index - 1
            else:
                break
            try:
                nextGene = _parse_header(all_genes[next_index])
                if abs(seq_start - nextGene["start"]) > 8000:
                    break
                if geneStrand == "-" and nextGene["direction"] == "+" and direction == "+":
                    geneList.append(nextGene)
                elif geneStrand == "+" and nextGene["direction"] == "-" and direction == "-":
                    geneList.append(nextGene)
                elif geneStrand == nextGene["direction"]:
                    geneList.append(nextGene)
                index = next_index
            except Exception:
                break

    gene_strand = strand
    try:
        idx_down = reg_idx - 1
        down_gene = _parse_header(all_genes[idx_down])
        if strand == "+" and down_gene["direction"] == "-":
            gene_strand = down_gene["direction"]
        down_genes = [down_gene]
        _walk(gene_strand, "-", down_gene, down_genes, idx_down)
        gene_array = list(reversed(down_genes))
    except Exception:
        gene_array = []

    gene_array.append(_parse_header(all_genes[reg_idx]))
    regulator_index = len(gene_array) - 1
    gene_strand = strand

    try:
        idx_up = reg_idx + 1
        up_gene = _parse_header(all_genes[idx_up])
        if strand == "-" and up_gene["direction"] == "+":
            gene_strand = up_gene["direction"]
        gene_array.append(up_gene)
        _walk(gene_strand, "+", up_gene, gene_array, idx_up)
    except Exception:
        return gene_array, regulator_index

    return gene_array, regulator_index


# --- public interface ----------------------------------------------------

def fetch(homolog_dict: dict) -> Optional[dict]:
    """Original acc2operon algorithm. See module docstring for I/O contract."""
    if "Genome" not in homolog_dict:
        return None
    genes, reg_idx = _get_genes(
        homolog_dict["Genome"],
        int(homolog_dict["Start"]),
        int(homolog_dict["Stop"]),
    )
    if reg_idx is None:
        return None
    reg = _parse_header(genes[reg_idx])
    operon, regulator_index = _walk_operon(genes, reg_idx, reg["start"], reg["direction"])
    return {
        "operon": operon,
        "protein_index": regulator_index,
        "genome": homolog_dict["Genome"],
    }
