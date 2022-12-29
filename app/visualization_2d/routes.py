from fastapi import APIRouter
from .schemas import Protein
from .service import Visualization2DService

#|------------------------------------------------------------------------------|#

router = APIRouter()

@router.post("/api/tosmiles", tags=["converting"], description="Returns SMILES for a given sequence")
def get_2d_visualization(protein: Protein):
    visualizer = Visualization2DService(protein.sequence)

    return {"smiles": visualizer.protein_to_smiles()}

#|------------------------------------------------------------------------------|#

@router.post("/api/to2d", tags=["converting"], description="Returns 2D visualization for a given sequence")
def get_2d_visualization(protein: Protein):
    visualizer = Visualization2DService(protein.sequence)

    return {"svg": visualizer.protein_to_svg()}
