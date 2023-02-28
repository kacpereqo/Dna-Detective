from fastapi import APIRouter, Depends
from .schemas import PH_Range, PH, Pka_scale, Frame
from .service import ProteinCharge


# |------------------------------------------------------------------------------|#

router = APIRouter()


@router.post("/api/netcharge/", tags=["Charge"], description="Returns charge of protein in units")
def get_netcharge_of_protein(frame: Frame, pH_range: PH_Range = Depends()):
    """Returns charge of protein in units"""

    protein = ProteinCharge(frame.frame)
    return {"netcharge": protein.net_charge(start=pH_range.start, end=pH_range.end, step=pH_range.step)}

# |------------------------------------------------------------------------------|#


@router.post("/api/isoelectricpoint", tags=["Charge"], description="Returns isoelectric point of protein in units")
async def get_isoelectricpoint_of_protein(frame: Frame, pka_scale: Pka_scale = Depends()):
    """Returns isoelectric point of protein in units"""

    isoelectric_point = ProteinCharge(frame.frame).isoelectric_point()
    return {"isoelectricpoint": round(isoelectric_point, 2)}
