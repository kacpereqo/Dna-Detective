from fastapi import APIRouter
from .schemas import Protein
from .service import Visualization2DService
from fastapi.responses import FileResponse
from .schemas import Frame

# |------------------------------------------------------------------------------|#

router = APIRouter()


@router.post("/api/tosmiles", tags=["converting"], description="Returns SMILES for a given sequence")
def get_2d_visualization(protein: Protein):
    visualizer = Visualization2DService(protein.sequence)

    return {"smiles": visualizer.protein_to_smiles()}

# |------------------------------------------------------------------------------|#


@router.post("/api/to2d", tags=["converting"], description="Returns 2D visualization for a given sequence")
def get_2d_visualization(protein: Protein):
    visualizer = Visualization2DService(protein.sequence)
    visualizer.protein_to_svg()
    return FileResponse("test.png")


@router.post("/api/visualizaiton", tags=["converting"], description="Returns 2D visualization for a given sequence", response_class=FileResponse)
async def get_2d_visualization(frame : Frame):
    visualizer = Visualization2DService(frame.frame)
    await visualizer.protein_to_svg()
    return FileResponse("test.png", media_type="image/png", headers={'Access-Control-Expose-Headers': 'Content-Disposition'}, filename="test.png")
