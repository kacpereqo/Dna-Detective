from fastapi import APIRouter
from .schemas import Protein
from .service import ProteinStatistics
from typing import Optional


#|------------------------------------------------------------------------------|#

router = APIRouter()

#|------------------------------------------------------------------------------|#

@router.post("/api/percentage/{amine_acid}", tags=["Statistics"], description="Returns flexibility of protein")
def get_percentage_of_amine_acid(protein: Protein, amine_acid: str):
    """Returns flexibility of protein in units""" 
    protein  = ProteinStatistics(protein.sequence)
    percentage = protein.perctage_of(amine_acid)
    return {amine_acid:percentage}

#|------------------------------------------------------------------------------|#

@router.post("/api/percentage", tags=["Statistics"], description="Returns flexibility of protein")
def get_percentage_of_all_amine_acids(protein: Protein):
    """Returns flexibility of protein in units""" 
    protein  = ProteinStatistics(protein.sequence)
    percentage = protein.percentage_of_all()
    return percentage
