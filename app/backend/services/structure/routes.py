from .constants import ALPHA_HELIX_SCALE, MUTABILITY_SCALE
from services.common.getproperty import get_property
from services.common.schemas import Frame
from fastapi import APIRouter

router = APIRouter()


@router.post("/api/refractivity", tags=["AlphaHelix"], description="Returns alpha helix of protein")
def get_alpha_helix_of_protein(frame: Frame):
    """Returns alpha helix of protein in units"""
    return {"refractivity": get_property(ALPHA_HELIX_SCALE, frame.frame)}


@router.post("/api/mutability", tags=["mutability"], description="Returns alpha helix of protein in percentage")
def get_mutability_of_protein(frame: Frame):
    """Returns alpha helix of protein in percentage"""
    return {"mutability": get_property(MUTABILITY_SCALE, frame.frame)}
