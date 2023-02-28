from fastapi import APIRouter
from .schemas import Frame
from .service import ProteinRecognition

router = APIRouter()


@router.post("/api/recognition", tags=["properties"], description="Returns weight of protein in daltons")
def get_weight_of_protein(frame: Frame):
    """Returns weight of protein in units"""
    bulkiness = ProteinRecognition(frame.frame).get_polarity()
    return {"recognition": bulkiness}
