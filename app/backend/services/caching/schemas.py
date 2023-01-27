from pydantic import BaseModel, validator
from fastapi import HTTPException

# |------------------------------------------------------------------------------|#W


class Rna_to_translate(BaseModel):
    sequence: str

    @validator("sequence")
    def rna_is_valid(cls, v):
        if not v.isupper():
            v = v.upper()

        if ("T" in v):
            v = v.replace("T", "U")

        if any(x not in ["A", "U", "G", "C"] for x in v):
            raise HTTPException(
                status_code=422, detail="RNA must contain only A, U, G, C")
        return v


class Frame(BaseModel):
    frame: str
