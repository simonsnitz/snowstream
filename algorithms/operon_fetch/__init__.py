"""operon_fetch — given a regulator accession, return its surrounding operon.

Interface (all versions):
    fetch(homolog_dict: dict) -> dict | None
where homolog_dict has at minimum keys ``Genome``, ``Start``, ``Stop``,
``Strand`` (the IPG-derived CDS coordinates), and the return is
    {"operon": [...gene dicts...], "protein_index": int, "genome": str}
or None on failure.
"""

from . import v0  # noqa: F401
