from fastapi import APIRouter
from app.core.settings import DATA_PATH
from pathlib import Path 
import pandas as pd
import requests 

edge_color = '#f0f0f0'
selected_cpg_color = '#f0027f'
selected_edge_color = '#386cb0'

network = APIRouter()

chromosome_distance = {'1':0, 
'2': 249250621,
'3': 492449994,
'4': 690472424,
'5': 881626700,
'6': 1062541960,
'7': 1233657027,
'8': 1392795690,
'9': 1539159712,
'10': 1680373143,
'11': 1815907890,
'12': 1950914406,
'13': 2084766301,
'14': 2199936179,
'15': 2307285719,
'16': 2409817111,
'17': 2500171864,
'18': 2581367074,
'19': 2659444322,
'20': 2718573305,
'21': 2781598825,
'22': 2829728720,
'X': 2881033286,
'Y': 3036303846}


GODMC_API_URL = 'http://api.godmc.org.uk/v0.1/query'
EWAS_API_URL = 'http://ewascatalog.org/api/?cpg='


@network.get('/process')
# fix filters
def get_data(file:str,minDistance:int,minAssoc:int,chromosomeList:str,minPval:int):
    
    data = pd.read_csv(Path(DATA_PATH)/file,index_col=0) 
    annot = pd.read_csv(Path(DATA_PATH)/'annot'/'cpg_annotation.csv',index_col=0)
    chr_lst = chromosomeList.split('-')
    data['cpg_chr'] = data['cpg_chr'].apply(str)
    data['snp_chr'] = data['snp_chr'].apply(str)

    data = data[data['cistrans']==False]

    data = data.groupby('cpg') # group data by cpgs
    data = data.filter(lambda x: len(x) >= minAssoc) # filter by number of associations per cpg

    data = data[data['cpg_chr'].isin(chr_lst)]
    data = data[data['pval']<=1*10**-minPval]

    data['cpg_pos_abs'] = data['cpg_chr'].apply(lambda cpg_chr: chromosome_distance[str(cpg_chr)]) + data['cpg_pos'] # calculate absolute distance of cpg
    data['snp_pos_abs'] = data['snp_chr'].apply(lambda snp_chr: chromosome_distance[str(snp_chr)]) + data['snp_pos'] # calculate absolute distance of snp
    data['dist'] = abs(data['cpg_pos_abs'] - data['snp_pos_abs']) # calculate distance between pairs
    data = data[data['dist']>=minDistance] # filter by min distance

    data = data.reset_index(drop=True)

    data = data.merge(annot,how='left',on='cpg')

    data.drop(['MAPINFO','CHR'],axis=1,inplace=True)
    
    data['UCSC_RefGene_Name'] = data['UCSC_RefGene_Name'].fillna('no annotation')
    
    data['id'] = data.index 

    response = data.to_dict('records')
 
    return response

@network.get('/ewas')
def ewas(cpg:str,file:str,level:int):
    
    data = pd.read_csv(Path(DATA_PATH)/file,index_col=0)
   
    cpg_interest = [cpg]
    cpgs = []
    snps = []

    for cpg in cpg_interest:
        
        start = 1
        cpgs.append({'id':cpg,
                'level':0})

        while start <= level:
            
            cpg_included = [cpg['id'] for cpg in cpgs if cpg['level'] < start]
            snps_to_include = data[data['cpg'].isin(cpg_included)]['snp'].values
            [snps.append({'id':snp,'level':start}) for snp in snps_to_include 
                         if snp not in [snp['id'] for snp in snps]]
            snps_included = [snp['id'] for snp in snps if snp['level'] == start]
            cpgs_to_include = data[data['snp'].isin(snps_included)]['cpg'].values
            [cpgs.append({'id':cpg,'level':start}) for cpg in cpgs_to_include 
                         if cpg not in [cpg['id'] for cpg in cpgs]] 
            
            start+= 1

    networkDf = data[data['cpg'].isin([cpg['id'] for cpg in cpgs])]
    networkDf = networkDf[networkDf['snp'].isin([snp['id'] for snp in snps])]
    networkDf = networkDf.reset_index(drop=True)   
    networkDf['id'] = networkDf.index 
    edges_data = networkDf.to_dict('records')
    nodes = cpgs + snps
    networkObj = {'nodeAttributes':{'snp':snps,'cpg':cpgs},'edgeAttributes':edges_data}    
    
    try:
        
        ewasResponse = requests.get(EWAS_API_URL+cpg)
        ewasData = ewasResponse.json()
        
        try :
            uniqueCpgList = networkDf['cpg'].unique().tolist()
            godmcQueryJson = {"cpgs":uniqueCpgList,"cistrans": "trans","pval": 1e-10}
            godmcResponse = requests.post(GODMC_API_URL,json=godmcQueryJson)
            godmcData = godmcResponse.json()
            
            return {'subgraph':networkObj,'ewas':ewasData,'godmc':godmcData}
        
        except:
            return {'subgraph':networkObj,'ewas':ewasData,'godmc':'error'}
    
    except:
        
        return {'subgraph':networkObj,'ewas':'error','godmc':'error'}
        