from fastapi import APIRouter, Request
from services.auth.service import get_current_user
from services.database.db import DB
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
async def translate_rna(req: Request, _id: str, is_reversed: bool = False, is_forward: bool = True):
    berear = req.headers.get("Authorization")
    if berear is not None:
        try:
            user = await get_current_user(berear.split(" ")[1])
            if user is not None:
                DB().add_sequence_to_user(_id, user['_id'])
        except Exception as e:
            print(e)

    data = DB().get_sequence(_id)
    if "translation" in data:
        return data["translation"]
    else:
        translator = Translator(data["sequence"], is_reversed, is_forward)
        return translator.parse()

    # |------------------------------------------------------------------------------|#
