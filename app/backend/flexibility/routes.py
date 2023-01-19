from fastapi import APIRouter
from .schemas import Protein
from .service import ProteinFlexibility
from typing import Optional


#|------------------------------------------------------------------------------|#

router = APIRouter()

#|------------------------------------------------------------------------------|#

@router.post("/api/flexibility", tags=["Flexibility"], description="Returns flexibility of protein")
def get_flexibility_of_protein(protein: Protein, window_size: Optional[int] = 9):
    """Returns flexibility of protein in units"""
    
    sequence = ProteinFlexibility(protein.sequence)
    flexibility = sequence.flexibility(window_size = window_size)
    return {"flexibility":flexibility}