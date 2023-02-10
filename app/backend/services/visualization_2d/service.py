from .constants import AMINE_ACIDS_SMILES
import asyncio
import concurrent.futures
from rdkit import Chem
from rdkit.Chem import rdDepictor, Draw
from threading import Thread


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
        Draw.MolToFile(mol, "backend/visualizations/test.png",
                       (50 * len(self.sequence), 200))

    async def protein_to_svg(self):
        """Returns visualization of protein in SVG format"""

        loop = asyncio.get_running_loop()
        with concurrent.futures.ProcessPoolExecutor() as pool:
            await loop.run_in_executor(pool, self.draw)

        return "backend/visualizations/test.png"

    # |------------------------------------------------------------------------------|#
