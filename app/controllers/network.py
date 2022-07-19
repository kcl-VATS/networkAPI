from fastapi import HTTPException,APIRouter,File, UploadFile
from loguru import logger
from app.core.settings import DATA_PATH
from pathlib import Path 
import networkx as nx
import pandas as pd
from typing import Tuple
from pydantic import BaseModel

network = APIRouter()

cpgNet = nx.DiGraph()




@network.get('/load')
async def load(filename:str):
    if len(cpgNet)>0:
        return 'Network already loaded'
    data = pd.read_csv(Path(DATA_PATH)/filename,delimiter = "\t")
    cpgs:list[str] = data['CpG'].unique() # cpg nodes
    lds:list[str] = data['LD clump'].unique() # ld clump nodes
    cpg_ld_edges:list[Tuple] = [(cpg,'ld_'+ld) for cpg,ld in zip(data['CpG'],data['LD clump'])]
    cpgNet.add_nodes_from(cpgs)
    cpgNet.add_nodes_from(lds)
    cpgNet.add_edges_from(cpg_ld_edges)
    return 'Network Loaded'

@network.get('/subgraph')
async def subgraph(cpg:str):
    
    if len(cpgNet) == 0:
       return 'Empty network!'
    
    if not cpgNet.has_node(cpg):
        return 'CpG not in network'

    undirected_net:nx.classes.graph.Graph = cpgNet.to_undirected() # turn network to undirected to locate all CpG associations by LD
    allowed_jumps =  2
    reachable_nodes = [node[0] for node in nx.shortest_path_length(undirected_net,source=cpg).items() 
                       if node[1] <= allowed_jumps] 

    subGraph =  cpgNet.subgraph(reachable_nodes)

    return nx.to_dict_of_dicts(subGraph)
 

@network.get('/length')
async def length():
    return len(cpgNet)
