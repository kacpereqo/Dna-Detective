from fastapi import APIRouter
from .schemas import Protein
from localcider.sequenceParameters import SequenceParameters
from typing import Optional
from .service import ProteinPropeties


#|------------------------------------------------------------------------------|#

def renormalize(n, range1, range2):
    delta1 = range1[1] - range1[0]
    delta2 = range2[1] - range2[0]
    return (delta2 * (n - range1[0]) / delta1) + range2[0]

router = APIRouter()

@router.post("/api/weight", tags=["properties"], description="Returns weight of protein in daltons")
def get_weight_of_protein(protein: Protein):
    """Returns weight of protein in units"""
    sequence = SequenceParameters(protein.sequence)
    return {"weight":round(sequence.get_molecular_weight(),3)}
    
#|------------------------------------------------------------------------------|#

@router.post("/api/lenght", tags=["properties"], description="Returns lenght of protein in amino acids")
def get_lenght_of_protein(protein: Protein):
    """Returns charge of protein in units"""
    sequence = SequenceParameters(protein.sequence)
    return {"lenght":sequence.get_sequence_length()}

#|------------------------------------------------------------------------------|#

@router.post("/api/netcharge", tags=["properties"], description="Returns charge of protein in units")
def get_netcharge_of_protein(protein: Protein, pH: Optional[float] = 7.0):
    """Returns charge of protein in units"""
    sequence = SequenceParameters(protein.sequence)
    return {"netcharge":sequence.get_NCPR(pH=pH)}
    
#|------------------------------------------------------------------------------|#

@router.post("/api/isoelectricpoint", tags=["properties"], description="Returns isoelectric point of protein in units")
def get_isoelectricpoint_of_protein(protein: Protein):
    """Returns isoelectric point of protein in units"""
    pI = ProteinPropeties(protein.sequence).isoelectric_point()
    return {"isoelectricpoint":pI}

#|------------------------------------------------------------------------------|#

@router.post("/api/hydrophobicity-index", tags=["properties"], description="Returns hydrophobicity of protein in units")
def get_hydrophobicity_of_protein(protein: Protein):
    """Returns hydrophobicity of protein in units"""
    return {"hydrophobicity":round(ProteinPropeties(protein.sequence).hydrophobicity(), 3)}

#|------------------------------------------------------------------------------|#

@router.post("/api/hydrophobicity", tags=["properties"], description="Returns hydrophobicity of protein in units")
def get_hydrophobicity_of_protein(protein: Protein):
    """Returns hydrophobicity of protein in units"""
    sequence = SequenceParameters(protein.sequence)
    hydrophobicity = sequence.get_linear_hydropathy()
    response = dict(enumerate(hydrophobicity.flatten(), 1))
    return {"hydrophobicity":response}
