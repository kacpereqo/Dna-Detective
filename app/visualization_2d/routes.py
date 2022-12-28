from fastapi import APIRouter
from .schemas import Protein

router = APIRouter()

@router.post("/api/tosmiles", tags=["converting"], description="Returns SMILES for a given sequence")
def get_2d_visualization(sequence: Protein):
    return {"sequence": sequence}