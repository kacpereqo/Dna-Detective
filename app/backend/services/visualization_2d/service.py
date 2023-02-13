from .constants import AMINE_ACIDS_SMILES
from fastapi import BackgroundTasks
import asyncio
import concurrent.futures
from rdkit import Chem
from rdkit.Chem import rdDepictor, Draw
from threading import Thread
import multiprocessing


# |------------------------------------------------------------------------------|#


class Visualization2DService:
    def __init__(self, sequence):
        self.sequence = sequence
        self.smiles = self.protein_to_smiles()

    # |------------------------------------------------------------------------------|#

    def protein_to_smiles(self):
        """Returns SMILES for a given sequence"""
        return (
            "".join(AMINE_ACIDS_SMILES[amino_acid] for amino_acid in self.sequence)
            + "O"
        )
    # |------------------------------------------------------------------------------|#

    def draw(self):
        mol = Chem.MolFromSmiles(self.smiles)
        rdDepictor.SetPreferCoordGen(True)
        Draw.MolToFile(mol, "test.png",
                       (50 * len(self.sequence) , 200))

    def protein_to_svg(self):
        """Returns 2D visualization for a given sequence"""
        self.draw()
        return "test.png"

    # |------------------------------------------------------------------------------|#
