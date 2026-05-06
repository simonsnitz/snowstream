from typing import Literal, Optional, Any
from pydantic import BaseModel, Field


InputMethod = Literal["RefSeq", "Uniprot", "Protein sequence"]
CoordinatesMethod = Literal["batch", "individually"]
Database = Literal["local_diamond", "nr_remote"]


class BlastParams(BaseModel):
    ident_cutoff: int = Field(40, ge=30, le=90)
    cov_cutoff: int = Field(90, ge=50, le=100)
    max_homologs: int = Field(30, ge=10, le=100)
    filter_redundant: bool = True


class PromoterParams(BaseModel):
    min_length: int = Field(80, ge=1, le=500)
    max_length: int = Field(800, ge=20, le=9000)


class PredictRequest(BaseModel):
    input_method: InputMethod
    input_value: str
    blast_params: BlastParams = BlastParams()
    promoter_params: PromoterParams = PromoterParams()
    get_coordinates_method: CoordinatesMethod = "batch"
    database: Database = "local_diamond"
    force: bool = False


class ProteinInfo(BaseModel):
    annotation: Optional[str] = None
    organism: Optional[str] = None
    lineage: Optional[list[str]] = None


class HomologResult(BaseModel):
    uniprot_id: str
    identity: float
    coverage: float
    genome: Optional[str] = None
    start: Optional[int] = None
    stop: Optional[int] = None
    strand: Optional[str] = None
    operon: Optional[Any] = None
    promoter: Optional[str] = None


class PredictResult(BaseModel):
    cache_key: str
    cached_at: str
    input: PredictRequest
    protein_info: Optional[ProteinInfo] = None
    homologs: list[HomologResult]
