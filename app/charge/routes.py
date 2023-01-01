from fastapi import APIRouter
from .schemas import Protein, PH_Range, PH, Pka_scale
from typing import Optional
from .service import ProteinCharge


#|------------------------------------------------------------------------------|#

router = APIRouter()

#|------------------------------------------------------------------------------|#

@router.post("/api/netcharge", tags=["Charge"], description="Returns charge of protein in units")
def get_netcharge_of_protein(protein: Protein, start: Optional[float] = 0.0, end: Optional[float] = 14.0, step: Optional[float] = 0.5):
    """Returns charge of protein in units"""

    pH_range = PH_Range(start = start, end = end, step = step)
    protein = ProteinCharge(protein.sequence)
    return {"netcharge":protein.net_charge(start = pH_range.start,end = pH_range.end,step = pH_range.step)}
    
#|------------------------------------------------------------------------------|#

@router.post("/api/isoelectricpoint", tags=["Charge"], description="Returns isoelectric point of protein in units")
def get_isoelectricpoint_of_protein(protein: Protein):
    """Returns isoelectric point of protein in units"""

    isoelectric_point = ProteinCharge(protein.sequence).isoelectric_point()
    return {"isoelectricpoint": isoelectric_point}

#|------------------------------------------------------------------------------|#

@router.post("/api/charge", tags=["Charge"], description="Returns charge of protein at certain pH")
def get_charge_of_protein(protein: Protein, pka_scale: Optional[str] = "Rodwell", pH: Optional[float] = 7.0):
    """Returns charge of protein in units"""

    pka_scale = Pka_scale(scale = pka_scale) 
    pH = PH(pH = pH)

    sequence = ProteinCharge(protein.sequence)
    return {"charge":sequence.charge_at_ph(scale = pka_scale.scale, pH = pH.pH)}
