from fastapi import FastAPI
from app.core.settings import DEFAULT_ROUTE_STR
from app.controllers import router
from starlette.middleware.gzip import GZipMiddleware
from starlette.middleware.cors import CORSMiddleware

import uvicorn

# initate FastAPI instance
app = FastAPI()

# add middleware
app.add_middleware(GZipMiddleware, minimum_size=1000)

app.add_middleware(
    CORSMiddleware,
    allow_origin_regex="(http://(0\\.0\\.0\\.0|localhost)(:\\d+)?|)",
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# add router

app.include_router(router, prefix=DEFAULT_ROUTE_STR)

@app.get("/")
def ping():
    return {"ping": "live"}


if __name__ == "__main__":
    uvicorn.run(app, log_level="debug", reload=True)

