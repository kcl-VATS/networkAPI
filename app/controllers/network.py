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
def get_data(file:str,minDistance:int,minAssoc:int,minChrom:int):
    
    data = pd.read_csv(Path(DATA_PATH)/file,index_col=0) 
    
    data['cpg_pos_abs'] = data['cpg_chr'].apply(lambda cpg_chr: chromosome_distance[str(cpg_chr)]) + data['cpg_pos'] # calculate absolute distance of cpg
    data['snp_pos_abs'] = data['snp_chr'].apply(lambda snp_chr: chromosome_distance[str(snp_chr)]) + data['snp_pos'] # calculate absolute distance of snp
    data['dist'] = abs(data['cpg_pos_abs'] - data['snp_pos_abs']) # calculate distance between pairs

    df_g = data.groupby('cpg') # group data by cpgs
    df_g = df_g.filter(lambda x: len(x) >= minAssoc) # filter by number of associations per cpg
    df_g = df_g[df_g['dist']>=minDistance] # filter by min distance
    df_g = df_g.reset_index(drop=True)
    
    # filtering by minimum number of unique choromosomes
    #num_chrom_unique = df_g.groupby('cpg')['snp_chr'].nunique()
    #df_g = pd.merge(df_g,num_chrom_unique,on=['cpg','snp'])
   
    df_g = df_g[df_g['cpg_chr'] >= minChrom]
    df_g.reset_index()
    #df_g['inter'] = (df_g['cpg_chr'] == df_g['snp_chr'])
    #df_g = df_g[df_g['inter']== False]
    #df_g = df_g.reset_index(drop=True)
    df_g['id'] = df_g.index 
        
    return df_g.to_dict('records')


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
        
    



