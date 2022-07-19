import os 
import glob 

def get_files_list(DATA_PATH:str):
    return [os.path.basename(x) 
            for x in glob.glob(DATA_PATH+"/*.txt")]

