from fastapi import APIRouter
from services.database.db import DB
from .service import ProteinPolarity

router = APIRouter()


@router.get("/api/polarity/{_id}", tags=["properties"], description="Returns weight of protein in daltons")
def get_weight_of_protein(_id: int):
    """Returns weight of protein in units"""
    weight = ProteinPolarity(DB().get_frame(_id)).get_polarity()
    return {"polarity": weight}
