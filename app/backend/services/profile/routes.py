from fastapi import APIRouter, Depends
from fastapi.security import OAuth2PasswordBearer
from services.database.db import DB
from services.auth.service import get_current_user

router = APIRouter()
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token")

# |------------------------------------------------------------------------------|#


@router.get("/api/profile", tags=["Profile"], description="Get user profile")
async def get_profile(token: str = Depends(oauth2_scheme)):
    user = await get_current_user(token)
    sequences = DB().get_user_sequences(user["_id"])
    return {"sequences": sequences}

# |------------------------------------------------------------------------------|#
