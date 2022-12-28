from typing import  Optional, List, Dict
from pydantic import BaseModel, validator

class Rna_to_translate(BaseModel):
    rna: str

    @validator("rna")
    def rna_is_valid(cls, v):
        if not v.isupper():
            v = v.upper()

        if ("T" in v):
            v = v.replace("T", "U")
            
        if any(x not in ["A", "U", "G", "C"] for x in v):
            raise ValueError("RNA must contain only A, U, G, C")
        return v


class Rna_translated(BaseModel):
    frames: Optional[List[str]] = None
    translated_frames: Optional[List[str]] = None
    open_reading_frames: Dict[int, List[str]]

