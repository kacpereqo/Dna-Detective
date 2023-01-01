from fastapi import APIRouter
from .schemas import Protein
from .service import ProteinAromacity
from typing import Optional


#|------------------------------------------------------------------------------|#

router = APIRouter()

#|------------------------------------------------------------------------------|#

@router.post("/api/aromacity", tags=["Aromacity"], description="Returns aromacity of protein")
def get_flexibility_of_protein(protein: Protein):
    """Returns flexibility of protein in units"""
    
    sequence = ProteinAromacity(protein.sequence)
    aromacity = sequence.aromacity()
    return {"aromacity":aromacity}