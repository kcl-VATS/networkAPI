from fastapi import HTTPException,APIRouter,File, UploadFile
from loguru import logger
from pathlib import Path
from app.core.settings import DATA_PATH
import networkx as nx
import pandas as pd
from typing import Tuple

network = APIRouter()

cpgNet = nx.DiGraph()

@network.get('/load')
async def load(filename:str):
    if len(cpgNet)>0:
        return 'Network already loaded'
    data = pd.read_csv(DATA_PATH+'/'+filename,delimiter = "\t")
    cpgs:list[str] = data['CpG'].unique() # cpg nodes
    #snps:List[str]= data['Top SNP'].unique() # snp nodes
    lds:list[str] = data['LD clump'].unique() # ld clump nodes
    cpg_ld_edges:list[Tuple] = [(cpg,'ld_'+ld) for cpg,ld in zip(data['CpG'],data['LD clump'])]
    cpgNet.add_nodes_from(cpgs)
    cpgNet.add_nodes_from(lds)
    cpgNet.add_edges_from(cpg_ld_edges)
    return 'Success'

@network.get('/length')
async def length():
    return len(cpgNet)
