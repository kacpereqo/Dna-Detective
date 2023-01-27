from fastapi.middleware.cors import CORSMiddleware
from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from backend.config.include_routes import include_routes
from backend.services.database.db import DB


def get_application():
    _app = FastAPI()

    include_routes(_app)

    origins = ["*"]

    _app.add_middleware(
        CORSMiddleware,
        allow_origins=origins,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    return _app


app = get_application()


@app.exception_handler(ValueError)
async def value_error_exception_handler(request: Request, exc: ValueError):
    return JSONResponse(
        status_code=422,
        content={"message": str(exc)},
    )


@app.on_event("startup")
async def database_connect():
    DB().migrate()


@app.on_event("shutdown")
async def database_disconnect():
    DB().migrate()
