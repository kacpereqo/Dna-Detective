from .constants import RNA_CODON_TABLE
from typing import List
import re

#|------------------------------------------------------------------------------|#

class Translator():
    def __init__(self, rna: str, is_reversed: bool = False, is_forward: bool = True):
        self.rna = rna
        self.frames = self.get_frames(is_reversed, is_forward)
        self.translated_frames = self.translate_frames()
        self.open_reading_frames = self.find_open_reading_frames()

    #|------------------------------------------------------------------------------|#
     
    
    def get_frames(self, is_reversed: bool, is_forward: bool) -> List[str]:
        _frames = [] 

        if is_forward:
            _frames.extend(self.rna[i:] for i in range(3))

        if is_reversed:
            _frames.extend(self.rna[i:][::-1] for i in range(3))
        
        return _frames

    #|------------------------------------------------------------------------------|#

    def translate_frames(self):
        _proteins = []

        for frame in self.frames:
            protein = ""
            for i in range(0, len(frame), 3):
                codon = frame[i:i+3]
                if codon in RNA_CODON_TABLE:
                    protein += RNA_CODON_TABLE[codon]
                else:
                    break
            
            _proteins.append(protein)

        return _proteins        

    #|------------------------------------------------------------------------------|#

    def find_open_reading_frames(self):
        _proteins = {}
        for i,frame in enumerate(self.translated_frames):
            _proteins[i] = []
            for protein in frame.split("-"):
                if len(protein) > 0:
                    if protein[0] != "M":
                        _proteins[i].append(protein[0])
                    _proteins[i].extend(re.findall(r"M[A-Z]+", protein))

        return _proteins

#|------------------------------------------------------------------------------|#