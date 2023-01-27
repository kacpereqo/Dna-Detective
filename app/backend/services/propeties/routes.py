from fastapi import APIRouter
from .schemas import Protein
from localcider.sequenceParameters import SequenceParameters
from typing import Optional
from backend.services.database.db import DB


# |------------------------------------------------------------------------------|#

router = APIRouter()

# |------------------------------------------------------------------------------|#


@router.post("/api/weight", tags=["properties"], description="Returns weight of protein in daltons")
def get_weight_of_protein(protein: Protein):
    """Returns weight of protein in units"""
    sequence = SequenceParameters(protein.sequence)
    return {"weight": round(sequence.get_molecular_weight(), 3)}

# |------------------------------------------------------------------------------|#


@router.post("/api/lenght", tags=["properties"], description="Returns lenght of protein in amino acids")
def get_lenght_of_protein(protein: Protein):
    """Returns charge of protein in units"""
    sequence = SequenceParameters(protein.sequence)
    return {"lenght": sequence.get_sequence_length()}

# |------------------------------------------------------------------------------|#


@router.get("/api/lenght/{_id}", tags=["properties"], description="Returns charge of protein in units")
def get_charge_of_protein(_id: int):
    """Returns charge of protein in units"""
    sequence = SequenceParameters(DB().get_frame(_id))
    return {"charge": sequence.get_sequence_length()}

# |------------------------------------------------------------------------------|#


@router.get("/api/weight/{_id}", tags=["properties"], description="Returns weight of protein in daltons")
def get_weight_of_protein(_id: int):
    """Returns weight of protein in units"""
    sequence = SequenceParameters(DB().get_frame(_id))
    return {"weight": round(sequence.get_molecular_weight(), 3)}
