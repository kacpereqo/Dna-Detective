from fastapi import APIRouter
from backend.services.database.db import DB
from .schemas import Rna_to_translate, Frame

router = APIRouter(tags=["caching"])


@router.post("/api/sequence")
def cache_sequence(sequence: Rna_to_translate):
    """
    Caches a sequence in the database.

    :param sequence: The sequence to cache
    :return: The inserted row id
    """
    return DB().post_sequence(sequence.sequence)


@router.post("/api/frame")
def cache_sequence(frame: Frame):
    return DB().post_frame(frame.frame)
