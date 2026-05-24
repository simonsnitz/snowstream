"""operator_fetch — given a list of homologs each with a promoter, find the
inverted-repeat operator and return a consensus motif.

Interface:
    fetch(homologs: list[dict], params: dict) -> dict
where each homolog has at minimum {"Uniprot Id": str, "promoter": str}
and the return is:
    {
      "consensus_score": float,
      "num_seqs": int,
      "motif": list[{"base": str, "score": float}],
      "aligned_seqs": [...],
      "frequency_matrix": [...],
      "intergenic": str,
      ...
    }
"""

from . import v0  # noqa: F401
from . import v1  # noqa: F401
