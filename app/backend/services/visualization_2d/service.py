from .constants import AMINE_ACIDS_SMILES
from fastapi import BackgroundTasks
import asyncio
import concurrent.futures
from rdkit import Chem
from rdkit.Chem import rdDepictor, Draw


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
        rdDepictor.SetPreferCoordGen(True)
        mol = Chem.MolFromSmiles(self.smiles)
        Draw.MolToFile(mol, "test.png",
                       (50 * len(self.sequence) , 200))

    # |------------------------------------------------------------------------------|#

    async def protein_to_svg(self):
        """Returns 2D visualization for a given sequence"""
        loop = asyncio.get_running_loop()
        with concurrent.futures.ProcessPoolExecutor() as pool:
            await loop.run_in_executor(pool, self.draw)

        return "test.png"

    # |------------------------------------------------------------------------------|#
