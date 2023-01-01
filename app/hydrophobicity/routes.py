from fastapi import APIRouter, Depends
from .schemas import Protein, HydrophobicityScale
from localcider.sequenceParameters import SequenceParameters
from typing import Optional
from .service import ProteinHydrophobicity


#|------------------------------------------------------------------------------|#

router = APIRouter()

#|------------------------------------------------------------------------------|#

@router.post("/api/averagehydrophobicity", tags=["Hydrophobicity"], description="Returns average hydrophobicity of protein")
def get_hydrophobicity_of_protein(protein: Protein, hydrophobicity_scale: HydrophobicityScale = Depends()):
    """Returns hydrophobicity of protein in units"""
    
    protein = ProteinHydrophobicity(protein.sequence)

    return {"hydrophobicity":round(protein.hydrophobicity(scale = hydrophobicity_scale.scale), 2)}

#|------------------------------------------------------------------------------|#

@router.post("/api/hydrophobicity", tags=["Hydrophobicity"], description="Returns hydrophobicity of each amino acid in protein")
def get_hydrophobicity_of_protein(protein: Protein):
    """Returns hydrophobicity of protein in units"""

    sequence = SequenceParameters(protein.sequence)
    hydrophobicity = sequence.get_linear_hydropathy()
    response = dict(enumerate(hydrophobicity.flatten(), 1))
    return {"hydrophobicity":response}

#|------------------------------------------------------------------------------|#