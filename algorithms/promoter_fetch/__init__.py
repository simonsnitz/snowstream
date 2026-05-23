"""promoter_fetch — given an operon + regulator index, return the intergenic
region(s) that may contain the operator.

Interface:
    fetch(operon_data: dict, params: dict) -> str | list[str]
where operon_data = {"operon": [...], "protein_index": int, "genome": str}
and the return is either a single intergenic sequence (V0) or a list of
candidate sequences (V1+).

The pipeline runner in `algorithms/pipeline.py` is aware of both shapes.
"""

from . import v0  # noqa: F401
from . import v1  # noqa: F401
