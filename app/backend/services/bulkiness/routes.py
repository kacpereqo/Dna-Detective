from fastapi import APIRouter
from .schemas import Frame
from .service import ProteinBulkiness

router = APIRouter()


@router.post("/api/bulkiness", tags=["properties"], description="Returns weight of protein in daltons")
def get_weight_of_protein(frame: Frame):
    """Returns weight of protein in units"""
    bulkiness = ProteinBulkiness(frame.frame).get_polarity()
    return {"bulkiness": bulkiness}
