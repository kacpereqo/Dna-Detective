from pydantic import BaseModel, validator
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
            raise ValueError('Sequence contains invalid characters that are not amino acids')
        return v

#|------------------------------------------------------------------------------|#

class HydrophobicityScale(BaseModel):
    scale: str = 'Kyte-Doolittle'
    @validator('scale')
    def validate_scale(cls, v):
        if HYDROPHOBICITY_SCALE.get(v) is None:
            raise ValueError("Scale not found")
        return v 