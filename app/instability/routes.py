from fastapi import APIRouter
from .schemas import Protein
from .service import ProteinInstability
from typing import Optional


#|------------------------------------------------------------------------------|#

router = APIRouter()

#|------------------------------------------------------------------------------|#

@router.post("/api/instability", tags=["instability"], description="Returns instability of protein")
def get_instability_of_protein(protein: Protein):
    """Returns flexibility of protein in units"""
    
    sequence = ProteinInstability(protein.sequence)
    instability = sequence.instability()
    return {"instability":instability}