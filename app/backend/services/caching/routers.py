from fastapi import APIRouter
from services.database.db import DB
from .schemas import Sequence, Data, Frame

router = APIRouter()


@router.post("/api/sequence")
def cache_sequence(data : Data):
    """
    Caches a sequence in the database.

    :param sequence: The sequence to cache
    :return: The inserted row id
    """

    extension = data.data['extension']
    content = data.data['content']

    if extension == 'txt':
        sequence = Sequence(sequence=content)
        _id = DB().post_sequence(sequence.sequence)

    elif extension == 'fasta':
        sequence = Sequence(sequence=content)
        _id = DB().post_sequence(sequence.sequence)

    elif extension == 'csv':
        sequence = Sequence(sequence=content)
        _id = DB().post_sequence(sequence.sequence)

    elif extension == 'json':
        sequence = Sequence(sequence=content)
        _id = DB().post_sequence(sequence.sequence)

    elif extension == 'plaintext':
        sequence = Sequence(sequence=content)
        _id = DB().post_sequence(sequence.sequence)

    return {"id": _id}


@router.post("/api/frameid")
def cache_sequence(frame: Frame):
    return {"id": DB().post_frame(frame.frame)}
