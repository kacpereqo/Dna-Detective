from fastapi import APIRouter, Depends
from .schemas import Protein, HydrophobicityScale
from localcider.sequenceParameters import SequenceParameters
from peptides import Peptide
from .service import ProteinHydrophobicity
from backend.services.database.db import DB


# |------------------------------------------------------------------------------|#

router = APIRouter()

# |------------------------------------------------------------------------------|#


@router.post("/api/avghydrophobicity", tags=["Hydrophobicity"], description="Returns average hydrophobicity of protein")
def get_hydrophobicity_of_protein(protein: Protein, hydrophobicity_scale: HydrophobicityScale = Depends()):
    """Returns hydrophobicity of protein in units"""

    protein = ProteinHydrophobicity(protein.sequence)
    return {"hydrophobicity": round(protein.hydrophobicity(scale=hydrophobicity_scale.scale), 2)}

# |------------------------------------------------------------------------------|#


@router.post("/api/hydrophobicity", tags=["Hydrophobicity"], description="Returns hydrophobicity of each amino acid in protein")
def get_hydrophobicity_of_protein(protein: Protein, hydrophobicity_scale: HydrophobicityScale = Depends()):
    """Returns hydrophobicity of protein in units"""

    peptide = Peptide(protein.sequence)
    hydrophobicity = peptide .hydrophobicity_profile(window=3, scale="KyteDoolittle")
    return {"hydrophobicity": [round(x, 3) for x in hydrophobicity]}

# |------------------------------------------------------------------------------|#


@router.get("/api/hydrophobicity/{_id}", tags=["Hydrophobicity"], description="Returns hydrophobicity of each amino acid in protein")
def get_hydrophobicity_of_protein(_id: int, hydrophobicity_scale: HydrophobicityScale = Depends()):
    """Returns hydrophobicity of protein in units"""

    sequence = DB().get_frame(_id)

    peptide = Peptide(sequence)
    hydrophobicity = peptide .hydrophobicity_profile(window=3, scale="KyteDoolittle")
    return {"hydrophobicity": [round(x, 3) for x in hydrophobicity]}


@router.get("/api/avghydrophobicity/{_id}", tags=["Hydrophobicity"], description="Returns average hydrophobicity of protein")
def get_hydrophobicity_of_protein(_id: int, hydrophobicity_scale: HydrophobicityScale = Depends()):
    """Returns hydrophobicity of protein in units"""

    sequence = DB().get_frame(_id)

    protein = ProteinHydrophobicity(sequence)
    return {"hydrophobicity": round(protein.hydrophobicity(scale=hydrophobicity_scale.scale), 2)}
