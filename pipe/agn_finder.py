import os
import numpy as np

from astropy.io import ascii

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import match_coordinates_sky, SkyCoord
from astroquery.vizier import Vizier

import sqlite3
import pandas as pd


from pipe.dataprocess_config import *
from pipe import process_tools



class AgnFinder :

    def __init__(self, INPUTVAR: list) -> None:
        self.field, self.year, self.tel, self.filt = INPUTVAR
        self.kamp_fname = f'matched_targets_{self.field}.csv'
        self.quaia_cat_fname = f'quaia_{self.field}.csv'
        self.kamp_quaia_fname = f'kamp_quaia_{self.field}.csv'
        pass

    def _load_kamp_cat(self) -> None :
        '''
        Reads KMTNet source catalog file as astropy.table.table.Table format.
        KMTNet source catalog should be named as 'matched_targets_{field}.csv'
        '''
        kamp_cat_path = os.path.join(DATFDIRC, self.kamp_fname)
        self.kamp_cat = ascii.read(kamp_cat_path)
        self.cent_ra = np.median(self.kamp_cat['ra'])
        self.cent_dec = np.median(self.kamp_cat['dec'])
    
    def _load_quaia_cat(self) -> None :
        

        if self.quaia_cat_fname not in os.listdir(DATFDIRC) :

            print(f'{self.quaia_cat_fname} not in local directory. Downloading from Vizier ...')
            Vizier.ROW_LIMIT = -1

            gaia = Vizier.query_region(
                coord.SkyCoord(
                    ra=self.cent_ra, dec=self.cent_dec, 
                    unit=(u.deg, u.deg), frame='icrs'
                ),
                radius=2*u.deg,
                catalog=['Gaia'],
            )

            quaia_cat = gaia['I/356/qsocand']

            output_fpath = os.path.join(DATFDIRC, self.quaia_cat_fname)
            ascii.write(quaia_cat, output_fpath, overwrite=True)

        quaia_cat_path = os.path.join(DATFDIRC, self.quaia_cat_fname)
        self.quaia_cat = ascii.read(quaia_cat_path)


    def match_kamp_quaia(self) -> None :

        kamp_quaia_path = os.path.join(DATFDIRC, self.kamp_quaia_fname)

        if not os.path.isfile(kamp_quaia_path) : 

            self._load_kamp_cat()
            self._load_quaia_cat()

            kamp_radec = SkyCoord(
                ra=self.kamp_cat['ra']*u.deg,
                dec=self.kamp_cat['dec']*u.deg,
            )

            gaia_cat_radec = SkyCoord(
                ra=self.quaia_cat['RA_ICRS']*u.deg,
                dec=self.quaia_cat['DE_ICRS']*u.deg,
                frame='icrs'
            )

            mcat = self.kamp_cat.copy()

            idx, d2d, d3d = kamp_radec.match_to_catalog_sky(gaia_cat_radec)

            mcat['sep'] = d2d.to(u.arcsec)

            keys = ["Source", "SolID", "Class", "PQSO", "PGal", "z"]
            key_dict = {key: "quaia_"+key for key in keys}
            for key in key_dict :
                mcat[key_dict[key]] = self.quaia_cat[idx][key]

            kamp_quaia_cat = mcat[np.where(mcat['sep'] < 2*u.arcsec)]
            ascii.write(kamp_quaia_cat, kamp_quaia_path, overwrite=True)

        self.kamp_quaia_cat = ascii.read(kamp_quaia_path)


def gen_lc(agnname, tel, band) -> None :

    lc_dir = os.path.join(OUTPDIRC, 'agn_lightcurves')
    db_fnames = [f for f in os.listdir(OUTPDIRC) if f.endswith('.db')]
    lc_ls = []

    for db_fname in db_fnames :
        db_fpath = os.path.join(OUTPDIRC, db_fname)
        conn = sqlite3.connect(db_fpath)
        cursor = conn.cursor()
        
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables_list = [i[0] for i in cursor.fetchall()]
        
        tb_name = f'{tel}_{band}'
        if tb_name in tables_list :
            sql_prompt = f"SELECT * FROM {tb_name} WHERE Jname = '{agnname}'"
            lc = pd.read_sql(sql_prompt, conn)
            lc_ls.append(lc)
        
        cursor.close()
        conn.close()

    df = pd.concat(lc_ls)
    if len(df) > len(lc_ls)*30 :
        lc_fname = f'{agnname}_{tel}_{band}.csv'
        lc_fpath = os.path.join(lc_dir, lc_fname)
        df.to_csv(lc_fpath)


def run(INPUTVAR) :
    AF = AgnFinder(INPUTVAR)
    AF.match_kamp_quaia()
    for agnname in AF.kamp_quaia_cat['Jname'] :
        for band in ['B', 'V', 'I'] :
            gen_lc(agnname, 'CTIO', band)