from typing import Optional, List, Dict
from pydantic import BaseModel, validator
from fastapi import HTTPException

# |------------------------------------------------------------------------------|#W


class Rna_to_translate(BaseModel):
    rna: str

    @validator("rna")
    def rna_is_valid(cls, v):
        if not v.isupper():
            v = v.upper()

        v.replace(" ", "")

        if ("T" in v):
            v = v.replace("T", "U")

        if any(x not in ["A", "U", "G", "C"] for x in v):
            raise HTTPException(
                status_code=422, detail="RNA must contain only A, U, G, C")
        return v

# |------------------------------------------------------------------------------|#


class Rna_translated(BaseModel):
    frames: Optional[Dict[str, List[str]]] = None
    translated_frames: Optional[Dict[str, List[str]]] = None
    open_reading_frames: Dict[str, Dict[str, List[str]]]
