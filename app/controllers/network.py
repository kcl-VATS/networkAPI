from fastapi import APIRouter
from app.core.settings import DATA_PATH
from app.controllers.helpers import get_files_list
import networkx.algorithms.components as nc 
from pathlib import Path 
import networkx as nx
import pandas as pd
from typing import Tuple
import numpy as np

cpg_color = '#7fc97f'
ld_color = '#beaed4'
edge_color = '#fdc086'
selected_cpg_color = '#f0027f'
selected_edge_color = '#386cb0'

network = APIRouter()
cpgNet = nx.DiGraph()

@network.get('/load')
async def load(file:str):
    filename = file
    if len(cpgNet)>0:
        return 'Network already loaded'
    if not filename in  get_files_list(DATA_PATH):
        return 'No such file in the system'
                    
    data = pd.read_csv(Path(DATA_PATH)/filename,delimiter = "\t")
    cpgs:list[str] = data['CpG'].unique() # cpg nodes
    lds:list[str] = data['LD clump'].unique() # ld clump nodes
    lds = ['ld_'+ld for ld in lds]
    cpg_ld_edges:list[Tuple] = [(cpg,'ld_'+ld) 
                                for cpg,ld in zip(data['CpG'],data['LD clump'])]
    cpgNet.add_nodes_from(cpgs,color=cpg_color)
    cpgNet.add_nodes_from(lds,color=ld_color)
    cpgNet.add_edges_from(cpg_ld_edges,color=edge_color)
    weak_components:list = [list(weak_component) 
                    for weak_component in nc.weakly_connected_components(cpgNet)]
    cpg_ld_singles:list = np.array([cpg_ld 
                    for cpg_ld in weak_components if len(cpg_ld)<=2]).flatten()
    cpgNet.remove_nodes_from(cpg_ld_singles)
    return 'Network loaded'

@network.get('/empty_network')
async def empty():
    if len(cpgNet) == 0: 
        return 'Net is already empty'
    cpgNet.remove_edges_from([edges for edges in cpgNet.edges]) # remove all edges
    cpgNet.remove_nodes_from([nodes for nodes in cpgNet.nodes]) # remove all nodes
    return 'Net is emptied'


@network.get('/subgraph')
async def subgraph(cpg:str,numJumps:int):
    if len(cpgNet) == 0:
       return 'Empty network!'
    
    if not cpgNet.has_node(cpg):
        return 'CpG not in network'

    undirected_net:nx.classes.graph.Graph = cpgNet.to_undirected() 
    allowed_jumps:int =  numJumps
    reachable_nodes:list[str] = [node[0] for node in 
                        nx.shortest_path_length(undirected_net,source=cpg).items() 
                        if node[1] <= allowed_jumps] 

    subGraph:nx.classes.digraph.DiGraph =  cpgNet.subgraph(reachable_nodes)
    nx.set_node_attributes(subGraph,{cpg:selected_cpg_color},'color')
    edges_dict = {edges : selected_edge_color 
                for edges in subGraph.edges(cpg)}
    nx.set_edge_attributes(subGraph,edges_dict,'color')
    
    node_color =nx.get_node_attributes(subGraph, "color")
    edge_color =nx.get_edge_attributes(subGraph, "color")
    graphologyNodes:list =[{'key':node,'attributes':{
                                'color':node_color[node],
                                'size': subGraph.degree(node)/2,
                                'label':node
                                }} 
                            for node in list(subGraph.nodes())]
    graphologyEdges:list =[{'source':edge[0],'target':edge[1],
                            'attributes':{'color':edge_color[edge]}} 
                            for edge in list(subGraph.edges())]
    graphologyAttribute:dict = {'name': f'Network for {cpg}'} 
    graphologyObject:dict = {'attributes': graphologyAttribute,
    'nodes':graphologyNodes,
    'edges':graphologyEdges}

    return graphologyObject

@network.get('/connected_subgraph')
def connected_subgraph(cpg:str):
    undirected_net = cpgNet.to_undirected() 
    nodes = nc.node_connected_component(undirected_net,cpg)
    subGraph:nx.classes.digraph.DiGraph = cpgNet.subgraph(nodes)
    
     # change the color of source node
    nx.set_node_attributes(subGraph,{cpg:'#f0027f'},'color')

    # change the color of source edges

    edges_dict = {edges : '#386cb0' for edges in subGraph.edges(cpg)}
    
    nx.set_edge_attributes(subGraph,edges_dict,'color')

    node_color =nx.get_node_attributes(subGraph, "color")
    edge_color =nx.get_edge_attributes(subGraph, "color")

    graphologyNodes:list =[{'key':node,'attributes':{
                                'color':node_color[node],
                                'size': subGraph.degree(node)/2,
                                'label':node
                                }} 
                            for node in list(subGraph.nodes())]
    graphologyEdges:list =[{'source':edge[0],'target':edge[1],
                            'attributes':{'color':edge_color[edge]}} 
                            for edge in list(subGraph.edges())]
    graphologyAttribute:dict = {'name': f'Network for {cpg}'} 
    graphologyObject:dict = {'attributes': graphologyAttribute,
    'nodes':graphologyNodes,
    'edges':graphologyEdges}

    return graphologyObject

