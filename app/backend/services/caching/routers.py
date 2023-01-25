from fastapi import APIRouter
router = APIRouter(tags=["caching"])


# @router.post("/sequence")
# async def cache_sequence(sequence: str):
#     """
#     Caches a sequence in the database.

#     :param sequence: The sequence to cache
#     :return: The inserted row id
#     """

#     return await db.post_sequence(sequence)
