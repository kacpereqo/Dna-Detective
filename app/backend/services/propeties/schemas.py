from pydantic import BaseModel, validator

# |------------------------------------------------------------------------------|#


class Protein(BaseModel):
    sequence: str

    @validator('sequence')
    def validate_sequence(cls, v):
        if v.islower():
            v = v.upper()
        if any(char not in 'ARNDCQEGHILKMFPSTWYV' for char in v):
            raise ValueError('Sequence contains invalid amino acid')
        return v

# |------------------------------------------------------------------------------#


class Frame(BaseModel):
    frame: str

    @validator('frame')
    def validate_frame(cls, v):
        if v.islower():
            v = v.upper()
        if any(char not in 'ARNDCQEGHILKMFPSTWYV' for char in v):
            raise ValueError('Sequence contains invalid amino acid')
        return v
