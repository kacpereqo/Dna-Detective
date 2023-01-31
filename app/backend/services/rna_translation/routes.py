from fastapi import APIRouter, Body, Form

from backend.services.database.db import DB
from .service import Translator
from .schemas import Rna_to_translate, Rna_translated


# |------------------------------------------------------------------------------|#

router = APIRouter()


@router.post("/api/translate", tags=["Rna translation"], description="Translates RNA to proteins")
async def translate_rna(rna: Rna_to_translate, is_reversed: bool = False, is_forward: bool = True):

    translator = Translator(rna.rna, is_reversed, is_forward)

    return Rna_translated(
        frames=translator.frames,
        open_reading_frames=translator.open_reading_frames,
        translated_frames=translator.translated_frames
    )

# |------------------------------------------------------------------------------|#


@router.get("/api/{_id}/translate", tags=["Rna translation"], description="Translates RNA to proteins")
async def translate_rna(_id: int, is_reversed: bool = False, is_forward: bool = True):
    sequence = DB().get_sequence(_id)
    translator = Translator(sequence, is_reversed, is_forward)

    return translator.parse()

    # |------------------------------------------------------------------------------|#
