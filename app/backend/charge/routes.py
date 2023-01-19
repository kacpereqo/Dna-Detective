from fastapi import APIRouter, Depends
from .schemas import Protein, PH_Range, PH, Pka_scale
from typing import Optional
from .service import ProteinCharge


#|------------------------------------------------------------------------------|#

router = APIRouter()

#|------------------------------------------------------------------------------|#

@router.post("/api/netcharge", tags=["Charge"], description="Returns charge of protein in units")
def get_netcharge_of_protein(protein: Protein, pH_range: PH_Range = Depends()):
    """Returns charge of protein in units"""

    protein = ProteinCharge(protein.sequence)
    return {"netcharge":protein.net_charge(start = pH_range.start,end = pH_range.end,step = pH_range.step)}
    
#|------------------------------------------------------------------------------|#

@router.post("/api/isoelectricpoint", tags=["Charge"], description="Returns isoelectric point of protein in units")
def get_isoelectricpoint_of_protein(protein: Protein, pka_scale: Pka_scale = Depends()):
    """Returns isoelectric point of protein in units"""

    isoelectric_point = ProteinCharge(protein.sequence).isoelectric_point()
    return {"isoelectricpoint": isoelectric_point}

#|------------------------------------------------------------------------------|#

@router.post("/api/chargeatph", tags=["Charge"], description="Returns charge of protein at certain pH")
def get_charge_of_protein(protein: Protein, pka_scale: Pka_scale = Depends(), pH: PH = Depends()):
    """Returns charge of protein in units"""

    sequence = ProteinCharge(protein.sequence)
    return {"charge":sequence.charge_at_ph(scale = pka_scale.scale, pH = pH.pH)}
