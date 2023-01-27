from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor, rdMolDescriptors, Draw

from .constants import AMINE_ACIDS_SMILES

# |------------------------------------------------------------------------------|#


class Visualization2DService:
    def __init__(self, sequence):
        self.sequence = sequence
        rdDepictor.SetPreferCoordGen(True)

    # |------------------------------------------------------------------------------|#

    def protein_to_smiles(self):
        """Returns SMILES for a given sequence"""
        return (
            "".join(AMINE_ACIDS_SMILES[amino_acid] for amino_acid in self.sequence)
            + "O"
        )
    # |------------------------------------------------------------------------------|#

    def protein_to_svg(self):
        """Returns visualization of protein in SVG format"""

        smiles = self.protein_to_smiles()
        mol = Chem.MolFromSmiles(smiles)

        Draw.MolToFile(mol, "backend/visualizations/test.png",
                       size=((100 * len(self.sequence)), 200))

        # with open('test.svg', 'w') as f:
        #     f.write(drawer.GetDrawingText())

        # return drawer.GetDrawingText()

    # |------------------------------------------------------------------------------|#