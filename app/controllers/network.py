from fastapi import APIRouter
from app.core.settings import DATA_PATH
from pathlib import Path 
import pandas as pd

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


@network.get('/trans')
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
    data = data.query(f"pval <= {1*10**-minPval}") 
    data['cpg_pos_abs'] = data['cpg_chr'].apply(lambda cpg_chr: chromosome_distance[str(cpg_chr)]) + data['cpg_pos'] # calculate absolute distance of cpg
    data['snp_pos_abs'] = data['snp_chr'].apply(lambda snp_chr: chromosome_distance[str(snp_chr)]) + data['snp_pos'] # calculate absolute distance of snp
    data['dist'] = abs(data['cpg_pos_abs'] - data['snp_pos_abs']) # calculate distance between pairs
    data = data[data['dist']>=minDistance] # filter by min distance

    data = data.reset_index(drop=True)

    data = data.merge(annot,how='left',on='cpg')

    data.drop(['MAPINFO','CHR'],axis=1,inplace=True)
    
    data['UCSC_RefGene_Name'] = data['UCSC_RefGene_Name'].fillna('-')
    
    data['id'] = data.index 

    response = data.to_dict('records')
  
    return response


@network.get('/cis')
# fix filters
def get_data(file:str,minAssoc:int,targetChr:str,startPos:str,endPos:str,minPval:int):
    
    data = pd.read_csv(Path(DATA_PATH)/file,index_col=0) 
    annot = pd.read_csv(Path(DATA_PATH)/'annot'/'cpg_annotation.csv',index_col=0)
    data['cpg_chr'] = data['cpg_chr'].apply(str)
    data['snp_chr'] = data['snp_chr'].apply(str)
    data = data.query("cistrans == True")
    data = data.query(f'snp_chr == "{targetChr}"')
    data = data.query(f"snp_pos > {startPos} & snp_pos < {endPos}")
    data = data.groupby('cpg') # group data by cpgs
    data = data.filter(lambda x: len(x) >= minAssoc) # filter by number of associations per cpg
    data = data.query(f"pval <= {1*10**-minPval}") 
    data = data.reset_index(drop=True)
    data = data.merge(annot,how='left',on='cpg')
    data.drop(['MAPINFO','CHR'],axis=1,inplace=True) 
    data['UCSC_RefGene_Name'] = data['UCSC_RefGene_Name'].fillna('-')
    data['id'] = data.index 
    response = data.to_dict('records')
    return response




@network.get('/ewas')
def ewas(cpg:str,file:str,level:int):
    
    data = pd.read_csv(Path(DATA_PATH)/file,index_col=0)
    ewas_res = pd.read_csv(Path(DATA_PATH)/'annot'/'ewas_results.csv',dtype={'Chr':str},index_col=0)
    ewas_stud = pd.read_csv(Path(DATA_PATH)/'annot'/'ewas_stud.csv',index_col = 0)
    
    cpg_interest = cpg.split('-')
    cpgs = []
    snps = []
    level = 1    
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

    network_df = data[data['cpg'].isin([cpg['id'] for cpg in cpgs])]
    network_df = network_df[network_df['snp'].isin([snp['id'] for snp in snps])]
    network_df = network_df.reset_index(drop=True)   
    network_df['id'] = network_df.index 
    edges_data = network_df.to_dict('records')
    network_obj = {'nodeAttributes':{'snp':snps,'cpg':cpgs},'edgeAttributes':edges_data}
    ewas_res = ewas_res[ewas_res['CpG'].isin(network_df['cpg'].unique())]
    ewas_res = ewas_res.merge(ewas_stud,on=["StudyID"],how="left")
    ewas_res.reset_index(drop=True,inplace=True)
    ewas_res.fillna('-',inplace=True)
    ewas_fields = ewas_res.columns.tolist()
    ewas_values = ewas_res.values.tolist()
    ewas_obj = {"results":ewas_values,"fields":ewas_fields}

    
    return {'subgraph':network_obj,'ewas':ewas_obj,'godmc':'error'}
        
        

        