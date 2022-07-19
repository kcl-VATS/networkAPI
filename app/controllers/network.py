from fastapi import APIRouter
from loguru import logger
from app.core.settings import DATA_PATH
from pathlib import Path 
import networkx as nx
import pandas as pd
from typing import Tuple
import os 
import glob 

network = APIRouter()

cpgNet = nx.DiGraph()

@network.get('/load')
async def load(filename:str):
    if len(cpgNet)>0:
        return 'Network already loaded'
    if not filename in  [os.path.basename(x) 
                       for x in glob.glob(DATA_PATH+"/*.txt")]:
        return 'No such file in the system'
                    
    data = pd.read_csv(Path(DATA_PATH)/filename,delimiter = "\t")
    cpgs:list[str] = data['CpG'].unique() # cpg nodes
    lds:list[str] = data['LD clump'].unique() # ld clump nodes
    cpg_ld_edges:list[Tuple] = [(cpg,'ld_'+ld) 
                                for cpg,ld in zip(data['CpG'],data['LD clump'])]
    cpgNet.add_nodes_from(cpgs)
    cpgNet.add_nodes_from(lds)
    cpgNet.add_edges_from(cpg_ld_edges)
    return 'Network loaded'


@network.get('/empty_network')
async def empty():
    if len(cpgNet) == 0: 
        return 'Net is already empty'
    cpgNet.remove_edges_from([edges for edges in cpgNet.edges]) # remove all edges
    cpgNet.remove_nodes_from([nodes for nodes in cpgNet.nodes]) # remove all nodes
    return 'Net is emptied'


@network.get('/subgraph')
async def subgraph(cpg:str):
    if len(cpgNet) == 0:
       return 'Empty network!'
    
    if not cpgNet.has_node(cpg):
        return 'CpG not in network'

    undirected_net:nx.classes.graph.Graph = cpgNet.to_undirected() 
    allowed_jumps:int =  2
    reachable_nodes:list[str] = [node[0] for node in 
                        nx.shortest_path_length(undirected_net,source=cpg).items() 
                        if node[1] <= allowed_jumps] 

    subGraph:nx.classes.digraph.DiGraph =  cpgNet.subgraph(reachable_nodes)

    return nx.to_dict_of_dicts(subGraph)


@network.get('/connectivity')
async def isConnected(cpgs:Tuple):
    return 'not connected'

