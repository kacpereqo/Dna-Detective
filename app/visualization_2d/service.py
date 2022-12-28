from .constants import AMINE_ACIDS_SMILES

def protein_to_smiles(protein):
    smiles = ""
    for amino_acid in protein:
        smiles += AMINE_ACIDS_SMILES[amino_acid]
    return smiles