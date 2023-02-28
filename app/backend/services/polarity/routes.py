from fastapi import APIRouter
from .schemas import Frame
from .service import ProteinPolarity

router = APIRouter()


@router.post("/api/polarity", tags=["properties"], description="Returns weight of protein in daltons")
def get_weight_of_protein(frame: Frame):
    """Returns weight of protein in units"""
    weight = ProteinPolarity(frame.frame).get_polarity()
    return {"polarity": weight}
