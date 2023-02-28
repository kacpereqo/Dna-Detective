from fastapi import APIRouter, Depends
from .schemas import Protein, HydrophobicityScale, Frame
from peptides import Peptide
from .service import ProteinHydrophobicity
from services.database.db import DB


# |------------------------------------------------------------------------------|#

router = APIRouter()

# |------------------------------------------------------------------------------|#


@router.post("/api/hydrophobicity", tags=["Hydrophobicity"], description="Returns hydrophobicity of each amino acid in protein")
def get_hydrophobicity_of_protein(frame: Frame, hydrophobicity_scale: HydrophobicityScale = Depends()):
    """Returns hydrophobicity of protein in units"""

    peptide = Peptide(frame.frame)
    hydrophobicity = peptide .hydrophobicity_profile(window=3, scale="KyteDoolittle")
    return {"hydrophobicity": [round(x, 3) for x in hydrophobicity]}


# |------------------------------------------------------------------------------|#


@router.post("/api/avghydrophobicity", tags=["Hydrophobicity"], description="Returns average hydrophobicity of protein")
def get_hydrophobicity_of_protein(frame: Frame, hydrophobicity_scale: HydrophobicityScale = Depends()):
    """Returns hydrophobicity of protein in units"""

    protein = ProteinHydrophobicity(frame.frame)
    return {"hydrophobicity": round(protein.hydrophobicity(scale=hydrophobicity_scale.scale), 2)}
