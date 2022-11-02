from fastapi import APIRouter
from app.core.settings import DATA_PATH
from app.controllers.helpers import get_files_list
import networkx.algorithms.components as nc 
from pathlib import Path 
import networkx as nx
import pandas as pd
from typing import Tuple
import numpy as np
from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path
from networkx.algorithms.shortest_paths.generic import has_path
import requests 


cpg_color = '#7fc97f'
ld_color = '#beaed4'
edge_color = '#f0f0f0'
selected_cpg_color = '#f0027f'
selected_edge_color = '#386cb0'

network = APIRouter()
cpgNet = nx.DiGraph()
columns_to_drop = ['index','A1','A2','n','N','dist','MAF','CpG pos_abs','SNP pos_abs']

chromosome_distance = [0,
 249250621,
 492449994,
 690472424,
 881626700,
 1062541960,
 1233657027,
 1392795690,
 1539159712,
 1680373143,
 1815907890,
 1950914406,
 2084766301,
 2199936179,
 2307285719,
 2409817111,
 2500171864,
 2581367074,
 2659444322,
 2718573305,
 2781598825,
 2829728720,
 2881033286,
 3036303846]


@network.get('/load')
async def load(file:str,minDistance:int,minAssoc:int):
    filename = file
    if len(cpgNet)>0:
        cpgNet.clear()
    data = pd.read_csv(Path(DATA_PATH)/filename,delimiter = "\t")
    dummy = data.copy()
    data['CpG pos_abs'] = dummy['CpG chr'].apply(lambda cpg_chr: chromosome_distance[cpg_chr-1]) + dummy['CpG pos'] # calculate absolute distance of cpg
    data['SNP pos_abs'] = dummy['SNP chr'].apply(lambda snp_chr: chromosome_distance[snp_chr-1]) + dummy['SNP pos'] # calculate absolute distance of snp
    data['dist'] = abs(data['CpG pos_abs'] - data['SNP pos_abs']) # calculate distance between pairs
    df_g = data.groupby('CpG') # group data by cpgs
    df_g = df_g.filter(lambda x: len(x) >= minAssoc) # filter by number of associations per cpg
    df_g = df_g[df_g['dist']>minDistance]  # filter by minimum distance
    df_g['inter'] = (df_g['CpG chr'] == df_g['SNP chr'])
    df_g = df_g[df_g['inter']== False]
    df_g.reset_index(inplace=True)
    df_g.drop(columns=columns_to_drop,inplace=True)
    df_g['id'] = df_g.index
    #cpgs:list[str] = df_g['CpG'].unique() # cpg nodes
    #lds:list[str] = df_g['Top SNP'].unique() # ld clump nodes
    #lds = ['ld_'+ld for ld in lds]
    #print(cpgs)
    #cpg_ld_edges:list[Tuple] = [(cpg,'ld_'+ld) 
    #                            for cpg,ld in zip(df_g['CpG'],df_g['Top SNP'])]
    #cpgNet.add_nodes_from(cpgs,color=cpg_color)
    #cpgNet.add_nodes_from(lds,color=ld_color)
    #cpgNet.add_edges_from(cpg_ld_edges,color=edge_color)

    return df_g.to_dict('records')


@network.get('/empty_network')
async def empty():
    if len(cpgNet) == 0: 
        return 'Net is already empty'
    cpgNet.remove_edges_from([edges for edges in cpgNet.edges]) # remove all edges
    cpgNet.remove_nodes_from([nodes for nodes in cpgNet.nodes]) # remove all nodes
    return 'Net is emptied'


#@network.get('/subgraph')
#async def subgraph():
    
#    if len(cpgNet) == 0:
#        return 'Empty network!'

    #if not cpgNet.has_node(cpg):
    #    return 'CpG not in network'
    #undirected_net:nx.classes.graph.Graph = cpgNet.to_undirected() 
    #allowed_jumps:int =  numJumps
    #reachable_nodes:list[str] = [node[0] for node in 
    #                    nx.shortest_path_length(undirected_net,source=cpg).items() 
    #                    if node[1] <= allowed_jumps] 

    #subGraph:nx.classes.digraph.DiGraph =  cpgNet.subgraph(reachable_nodes)
    #subGraph = nx.Graph(subGraph)

    #nx.set_node_attributes(cpgNet,{cpg:selected_cpg_color},'color')
    #edges_dict = {edges : selected_edge_color 
    #            for edges in cpgNet.edges(cpg)}
    #nx.set_edge_attributes(cpgNet,edges_dict,'color')
    
#    node_color =nx.get_node_attributes(cpgNet, "color")
#    edge_color =nx.get_edge_attributes(cpgNet, "color")

#    graphologyNodes:list =[{'key':node,'attributes':{
#                                'color':node_color[node],
#                                'size': np.sqrt(cpgNet.degree(node)),
#                                'label':node
#                                }} 
#                            for node in list(cpgNet.nodes())]
#    graphologyEdges:list =[{'source':edge[0],'target':edge[1],
#                            'attributes':{'color':edge_color[edge]}} 
#                           for edge in list(cpgNet.edges())]
#    graphologyAttribute:dict = {'name': 'Network'} 
#    graphologyObject:dict = {'attributes': graphologyAttribute,
#    'nodes':graphologyNodes,
#    'edges':graphologyEdges}

#    return graphologyObject

@network.get('/connected_subgraph')
def connected_subgraph(cpg:str):
    undirected_net = cpgNet.to_undirected() 
    nodes = nc.node_connected_component(undirected_net,cpg)
    subGraph:nx.classes.digraph.DiGraph = cpgNet.subgraph(nodes)
    single_nodes = [nodes for nodes in subGraph.nodes() if subGraph.degree<=1]
     # change the color of source node
    subGraph.remove_nodes_from(single_nodes)
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

@network.get('/path')
def connected_subgraph(source:str,target:str):
    undirected_net = cpgNet.to_undirected() 
    if(has_path(undirected_net,source,target)):
        subGraph = cpgNet.subgraph(
                bidirectional_shortest_path(undirected_net,source,target))



        node_color =nx.get_node_attributes(subGraph, "color")
        edge_color =nx.get_edge_attributes(subGraph, "color")

        graphologyNodes:list =[{'key':node,'attributes':{
                                    'color':node_color[node],
                                    'label':node
                                    }} 
                                for node in list(subGraph.nodes())]
        graphologyEdges:list =[{'source':edge[0],'target':edge[1],
                                'attributes':{'color':edge_color[edge]}} 
                                for edge in list(subGraph.edges())]
        graphologyAttribute:dict = {'name': f'Network'} 
        graphologyObject:dict = {'attributes': graphologyAttribute,
        'nodes':graphologyNodes,
        'edges':graphologyEdges}

        return graphologyObject
    else:
        return 'no path'


@network.get('/process')
def get_data(file:str,minDistance:int,minAssoc:int,minChrom:int):
    
    data = pd.read_csv(Path(DATA_PATH)/file,delimiter = "\t")
    data['CpG pos_abs'] = data['CpG chr'].apply(lambda cpg_chr: chromosome_distance[cpg_chr-1]) + data['CpG pos'] # calculate absolute distance of cpg
    data['SNP pos_abs'] = data['SNP chr'].apply(lambda snp_chr: chromosome_distance[snp_chr-1]) + data['SNP pos'] # calculate absolute distance of snp
    data['dist'] = abs(data['CpG pos_abs'] - data['SNP pos_abs']) # calculate distance between pairs
    df_g = data.groupby('CpG') # group data by cpgs
    df_g = df_g.filter(lambda x: len(x) >= minAssoc) # filter by number of associations per cpg
    df_g = df_g[df_g['dist']>=minDistance] # filter by min distance
    
    num_chrom_unique = df_g.groupby('CpG')['SNP chr'].nunique()
    df_g = pd.merge(df_g,num_chrom_unique,on='CpG')
    df_g = df_g[df_g['SNP chr_y'] >= minChrom]
   
    #df_g['inter'] = (df_g['CpG chr'] == df_g['SNP chr'])
    #df_g = df_g[df_g['inter']== False]
    df_g.reset_index(inplace=True)
    df_g.drop(columns=columns_to_drop,inplace=True)
    df_g['id'] = df_g.index
    
  
    return df_g.to_dict('records')


@network.get('/ewas')
def ewas(cpg:str):
    try:
        cpg = requests.get(f'http://ewascatalog.org/api/?cpg={cpg}')
        return cpg.json()
    except:
        print('error')
        return 'error'


@network.get('/subgraph')
def get_data(file:str,targetCpg:str):
    
    data = pd.read_csv(Path(DATA_PATH)/file,delimiter = "\t")
    cpg_cons = data[data['CpG'] == targetCpg]
    snps = cpg_cons['Top SNP'].values
    snp_cons = data[data['Top SNP'].isin(snps)]
    out = pd.concat([snp_cons,cpg_cons],ignore_index=True)
    out.drop_duplicates(inplace=True)
    
    return out.to_dict('records')
 
    

