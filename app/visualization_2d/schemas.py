from pydantic import BaseModel, validator
from fastapi import HTTPException

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