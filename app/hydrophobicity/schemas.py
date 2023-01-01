from pydantic import BaseModel, validator
from .constants import HYDROPHOBICITY_SCALE

#|------------------------------------------------------------------------------|#

class Protein(BaseModel):
    sequence: str

    @validator('sequence')
    def validate_sequence(cls, v):
        if v.islower():
            v = v.upper()
        if any(char not in 'ARNDCQEGHILKMFPSTWYV' for char in v):
            raise ValueError('Sequence contains invalid amino acid')
        return v

#|------------------------------------------------------------------------------|#

class HydrophobicityScale(BaseModel):
    scale: str

    @validator('scale')
    def validate_scale(cls, v):
        if v not in HYDROPHOBICITY_SCALE:
            raise ValueError('Scale not found')
        return v