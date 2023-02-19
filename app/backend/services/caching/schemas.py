from pydantic import BaseModel, validator
from fastapi import HTTPException
from typing import Dict

# |------------------------------------------------------------------------------|#W


class Sequence(BaseModel):
    sequence: str

    @validator("sequence")
    def rna_is_valid(cls, v):
        if len(v) < 3:
            raise HTTPException(
                status_code=422, detail="RNA must be at least 3 nucleotides long")
        if not v.isupper():
            v = v.upper()

        v = v.strip()

        if ("T" in v):
            v = v.replace("T", "U")

        if any(x not in ["A", "U", "G", "C"] for x in v):
            raise HTTPException(
                status_code=422, detail="RNA must contain only A, U, G, C")
        return v


class Frame(BaseModel):
    frame: str


class Data(BaseModel):
    data: Dict[str, str]
