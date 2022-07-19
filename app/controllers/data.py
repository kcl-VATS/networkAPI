from fastapi import HTTPException,APIRouter,File, UploadFile
from pathlib import Path
from app.core.settings import DATA_PATH
from pydantic import BaseModel
import aiofiles
from app.controllers.helpers import get_files_list


data = APIRouter()


class filesOut(BaseModel):
    fileList: list[str]



@data.get('/check/',response_model=filesOut)
async def check():
    """
    checks files uploaded to the system
    :return: list of file names , 'No file' if none available
    """   
    files:list[str] = get_files_list(DATA_PATH)
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
    if file.filename in get_files_list(DATA_PATH):
        return 'File already exists'
    else:
        async with aiofiles.open(DATA_PATH+'/'+file.filename, 'wb') as out_file:
            content = await file.read()  # async read
            await out_file.write(content)  # async write
        return 'File added succesfully'

@data.get('/remove/')
async def remove(filename:str):
    return 0 