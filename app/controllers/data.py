from fastapi import HTTPException,APIRouter,File, UploadFile
from loguru import logger
from pathlib import Path
from app.core.settings import DATA_PATH
import pandas as pd
import glob
import os
from pydantic import BaseModel
import aiofiles

data = APIRouter()


class filesOut(BaseModel):
    fileList: list[str]


@data.get('/check/',response_model=filesOut)
async def check():
    """
    checks files uploaded to the system
    :return: list of file names , 'No file' if none available
    """   
    files:list[str] = [os.path.basename(x) 
                       for x in glob.glob(DATA_PATH+"/*.txt")]
    if not files:
        raise HTTPException(status_code=404,
              detail='no files uploaded') 
    else:
        return ({'fileList':files})


@data.post('/upload/')
async def upload(file:UploadFile=File(...)):
    """
    uploads file to the local server
    :return: success if uploaded , fail if file already exists
    """      
    if file.filename in  [os.path.basename(x) 
                       for x in glob.glob(DATA_PATH+"/*.txt")]:
        return 'File already exists'
    else:
        async with aiofiles.open(DATA_PATH+'/'+file.filename, 'wb') as out_file:
            content = await file.read()  # async read
            await out_file.write(content)  # async write
        return 'File added succesfully'
