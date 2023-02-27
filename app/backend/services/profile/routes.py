from fastapi import APIRouter, Depends, status
from fastapi.security import OAuth2PasswordBearer
from services.database.db import DB
from services.auth.service import get_current_user
from typing import Dict

router = APIRouter()
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token")

# |------------------------------------------------------------------------------|#


@router.get("/api/profile", tags=["Profile"], description="Get user profile")
async def get_profile(token: str = Depends(oauth2_scheme)):
    user = await get_current_user(token)
    sequences = DB().get_user_sequences(user["_id"])
    return {"sequences": sequences}


@router.post("/api/preferences", tags=["Profile"], description="Update user profile")
async def update_preferences(token: str = Depends(oauth2_scheme), key : str = None, value: str = None):
    user = await get_current_user(token)
    DB().update_user_preferences(user["_id"], key, value)
    return status.HTTP_200_OK

# |------------------------------------------------------------------------------|#
