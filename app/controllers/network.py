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
    chr_lst = chromosomeList.split('-')
    data['cpg_chr'] = data['cpg_chr'].apply(str)
    data['snp_chr'] = data['snp_chr'].apply(str)
    
    data = data[data['cistrans']==False]
    
    data = data.groupby('cpg') # group data by cpgs
    data = data.filter(lambda x: len(x) >= minAssoc) # filter by number of associations per cpg
    
    data = data[data['cpg_chr'].isin(chr_lst)]
    data = data[data['pval']<1*10**-minPval]
    

    data['cpg_pos_abs'] = data['cpg_chr'].apply(lambda cpg_chr: chromosome_distance[str(cpg_chr)]) + data['cpg_pos'] # calculate absolute distance of cpg
    data['snp_pos_abs'] = data['snp_chr'].apply(lambda snp_chr: chromosome_distance[str(snp_chr)]) + data['snp_pos'] # calculate absolute distance of snp
    data['dist'] = abs(data['cpg_pos_abs'] - data['snp_pos_abs']) # calculate distance between pairs
    data = data[data['dist']>=minDistance] # filter by min distance

    data = data.reset_index(drop=True)
        
    data['id'] = data.index 

    response = data.to_dict('records')  
 
    return response

@network.get('/ewas')
def ewas(cpg:str,file:str,targetCpg:str):
     
    data = pd.read_csv(Path(DATA_PATH)/file,index_col=0)
    cpg_cons = data[data['cpg'] == targetCpg]
    snps = cpg_cons['snp'].values
    snp_cons = data[data['snp'].isin(snps)]
    assoc_df = pd.concat([snp_cons,cpg_cons],ignore_index=True)
    assoc_df.drop_duplicates(inplace=True)
    assoc_df = assoc_df.reset_index(drop=True)
    assoc_df['id'] = assoc_df.index 
    
    try:
        
        ewasResponse = requests.get(EWAS_API_URL+cpg)
        ewasData = ewasResponse.json()
        
        try :
            uniqueCpgList = assoc_df['cpg'].unique().tolist()
            godmcQueryJson = {"cpgs":uniqueCpgList,"cistrans": "trans"}
            godmcResponse = requests.post(GODMC_API_URL,json=godmcQueryJson)
        
            godmcData = godmcResponse.json()
            
            return {'subgraph':assoc_df.to_dict('records'),'ewas':ewasData,'godmc':godmcData}
        
        except:
            return {'subgraph':assoc_df.to_dict('records'),'ewas':ewasData,'godmc':'error'}
    
    except:
        
        return {'subgraph':assoc_df.to_dict('records'),'ewas':'error','godmc':'error'}
        