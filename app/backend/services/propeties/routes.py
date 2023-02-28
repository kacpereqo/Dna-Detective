from fastapi import APIRouter
from localcider.sequenceParameters import SequenceParameters
from services.database.db import DB
from .service import SequncePropeties
from .schemas import Frame


# |------------------------------------------------------------------------------|#

router = APIRouter()

# |------------------------------------------------------------------------------|#


@router.get("/api/lenght/{_id}", tags=["properties"], description="Returns charge of protein in units")
def get_charge_of_protein(_id: str):
    """Returns charge of protein in units"""
    sequence = SequenceParameters(DB().get_frame(_id))
    return {"charge": sequence.get_sequence_length()}

# |------------------------------------------------------------------------------|#


@router.post("/api/weight/", tags=["properties"], description="Returns weight of protein in daltons")
def get_weight_of_protein(frame : Frame):
    """Returns weight of protein in units"""
    weight = SequncePropeties(frame.frame).get_weight()
    return {"weight": weight}
