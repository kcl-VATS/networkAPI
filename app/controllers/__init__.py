from fastapi import APIRouter
from app.controllers.data import data
from app.controllers.network import network

router = APIRouter()

router.include_router(data, prefix="/data")
router.include_router(network,prefix='/network')


