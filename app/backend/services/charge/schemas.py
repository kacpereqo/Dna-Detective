from pydantic import BaseModel, validator
from .constants import PKA_SCALE
from typing import Optional

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

class PH_Range(BaseModel):
    start: float= 0.0 
    end: float = 14.0  
    step: float = 0.5  

    @validator('start')
    def validate_start(cls,v):
        if v < 0.0:
            raise ValueError('Start pH must be greater than 0.0')
        if v > 14.0:
            raise ValueError('Start pH must be less than 14.0')
        return v
    
    @validator('end')
    def validate_end(cls, v, values):
        if v > 14.0:
            raise ValueError('End pH must be less than 14.0')
        if v < 0:
            raise ValueError('End pH must be greater than 0.0')
        if v < values['start']:
            raise ValueError('End pH must be greater than start pH')
        return v

    @validator('step')
    def validate_step(cls, v, values):
        if v <= 0.0:
            raise ValueError('Step must be greater than 0.0')
        if v > 1.0:
            raise ValueError('Step must be less than 1.0')
        if v > values['end'] - values['start']:
            raise ValueError('Step must be less than the range')
        print(values)
        if v < 0.1:
            raise ValueError('Step must be greater than 0.1')
        return v


#|------------------------------------------------------------------------------|#

class PH(BaseModel):
    pH: float = 7.0 

    @validator('pH')
    def validate_pH(cls, v):
        if v < 0.0:
            raise ValueError('pH must be greater than 0.0')
        if v > 14.0:
            raise ValueError('pH must be less than 14.0')
        return v

#|------------------------------------------------------------------------------|#

class Pka_scale(BaseModel):
    scale: str = 'Rodwell'

    @validator('scale')
    def validate_scale(cls, v):
        if v not in PKA_SCALE:
            raise ValueError('pKa scale not found')
        return v