import os

NoneType = type(None)

from astropy.table import Table
from astropy.io import ascii

import pandas as pd

from tqdm import tqdm

from multiprocessing import Pool

import warnings
warnings.filterwarnings('ignore')

from pipe.dataprocess_config import *
from pipe import process_tools

import sqlite3


class SynthesizeLightcurve : 

    def __init__(self, INPUTVAR) :
        field, year, tel, filt,  = INPUTVAR
        self.data_dirs = process_tools.find_directories(INPUTVAR)
        self.filedata = None

        csv_name = '_'.join(INPUTVAR)+'.csv'
        self.csv_path  = os.path.join(OUTPDIRC, csv_name)
        db_name = f'{field}_{year}.db'
        self.db_path = os.path.join(OUTPDIRC, db_name)
        self.db_table_name = f'{tel}_{filt}'

    def _worker(self, mp_var) :
        source = mp_var
        
        if os.path.isfile(self.csv_path) == False :
            ascii.write(source, self.csv_path)
        
        else :
            with open(self.csv_path, mode='a') as f:
                f.seek(0, os.SEEK_END)  # Some platforms don't automatically seek to end when files opened in append mode
                Table(source).write(f, format='ascii.no_header')

        return 0
    
    def generate(self) :

        mp_var = [self.filedata[i] for i in range(len(self.filedata))]

        with Pool(CORENUM) as p :
            for output in p.imap(self._worker, mp_var) :
                if type(output) == NoneType :
                    print('SOME PROBLEM')

    def generate_all(self) :
        delmag_path = self.data_dirs['delmag']

        filenames = [i for i in os.listdir(delmag_path) if i.endswith('.delmag.cat')]
        tqdm_bar_fmt = '{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}'
        for i in tqdm(range(len(filenames)), total=len(filenames), bar_format=tqdm_bar_fmt) :
            filename = filenames[i]
            filepath = os.path.join(delmag_path, filename)
            self.filedata = ascii.read(filepath)
            self.generate()

    def update_db(self,) :
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        df = pd.read_csv(self.csv_path, sep=' ')

        df.to_sql(self.db_table_name,conn,if_exists='replace')

        cursor.close()
        conn.close()

def run(INPUTVAR) :
    GL = SynthesizeLightcurve(INPUTVAR)
    GL.generate_all()
    GL.update_db()


