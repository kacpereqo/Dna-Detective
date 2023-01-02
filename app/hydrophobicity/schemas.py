from pydantic import BaseModel, validator
from fastapi import HTTPException
from .constants import HYDROPHOBICITY_SCALE

#|------------------------------------------------------------------------------|#

class Protein(BaseModel):
    sequence: str

    @validator('sequence')
    def validate_sequence(cls, v):
        v = v.strip()
        if v.islower():
            v = v.upper()
        if any(char not in 'ARNDCQEGHILKMFPSTWYV' for char in v):
            raise HTTPException(status=422, detail='Sequence contains invalid characters that are not amino acids')
        return v

#|------------------------------------------------------------------------------|#

class HydrophobicityScale(BaseModel):
    scale: str
    @validator('scale')
    def validate_scale(cls, v):
        if HYDROPHOBICITY_SCALE.get(v) is None:
            raise HTTPException(status_code=422, detail="Scale not found")
        return v 