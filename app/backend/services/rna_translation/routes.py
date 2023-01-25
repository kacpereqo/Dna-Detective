from fastapi import APIRouter, Body, Form
from .service import Translator

from .schemas import Rna_to_translate, Rna_translated

# |------------------------------------------------------------------------------|#

router = APIRouter()


@router.post("/api/translate", tags=["Rna translation"], description="Translates RNA to proteins")
async def translate_rna(rna: Rna_to_translate, is_reversed: bool = False, is_forward: bool = True):

    # rna = await db.get_sequence(_id)
    # print(rna)
    translator = Translator(rna.rna, is_reversed, is_forward)

    return Rna_translated(
        frames=translator.frames,
        open_reading_frames=translator.open_reading_frames,
        translated_frames=translator.translated_frames
    )

# |------------------------------------------------------------------------------|#
