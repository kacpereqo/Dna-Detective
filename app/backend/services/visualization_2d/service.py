from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor, rdMolDescriptors, Draw
import multiprocessing
import concurrent.futures
import asyncio

from .constants import AMINE_ACIDS_SMILES

# |------------------------------------------------------------------------------|#


class Visualization2DService:
    def __init__(self, sequence):
        self.sequence = sequence

    # |------------------------------------------------------------------------------|#

    def protein_to_smiles(self):
        """Returns SMILES for a given sequence"""
        return (
            "".join(AMINE_ACIDS_SMILES[amino_acid] for amino_acid in self.sequence)
            + "O"
        )
    # |------------------------------------------------------------------------------|#

    def draw(self, mol):
        rdDepictor.SetPreferCoordGen(True)
        Draw.MolToFile(mol, "backend/visualizations/test.png",
                       ((100 * len(self.sequence)), 200))

    async def protein_to_svg(self):
        """Returns visualization of protein in SVG format"""

        smiles = self.protein_to_smiles()
        mol = Chem.MolFromSmiles(smiles)
        loop = asyncio.get_running_loop()
        with concurrent.futures.ProcessPoolExecutor() as pool:
            await loop.run_in_executor(pool, self.draw, mol)

        return "backend/visualizations/test.png"

        # with open('test.svg', 'w') as f:
        #     f.write(drawer.GetDrawingText())

        # return drawer.GetDrawingText()

    # |------------------------------------------------------------------------------|#
