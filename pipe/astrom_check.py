import os
import sys

import numpy as np

from astropy.io import ascii

from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord

from astropy.io import ascii

from pipe.dataprocess_config import *
from pipe import process_tools


class AstrometryChecker :
    # Global Path Setting --------------------------------------------------
    def __init__(self, INPUTVAR) : 
        field, year, tel, filter = INPUTVAR
        gaia_fpath = os.path.join(DATFDIRC, f'gaia_{field}.csv')
        self.gaia_cat = ascii.read(gaia_fpath)

        self.data_dirs = process_tools.find_directories(INPUTVAR)
        self.cat_fnames = [i for i in os.listdir(self.data_dirs['initial']) if i.endswith('.cat')]

        self.badastrom_dirs = {} 

        self.keys = ["image", "initial"]
        for key in self.keys :
            self.badastrom_dirs[key] = os.path.join(self.data_dirs[key], 'BADASTROM')
            os.makedirs(self.badastrom_dirs[key], exist_ok=True)


    # Modules --------------------------------------------------------------

    def match_catalog(self, cat) :
        
        data_coord = coord.SkyCoord(
            ra=cat['ALPHA_J2000'],
            dec=cat['DELTA_J2000'], 
            unit=(u.deg, u.deg)
        )
        
        refcat_coord = coord.SkyCoord(
            ra=self.gaia_cat['RAJ2000'], 
            dec=self.gaia_cat['DEJ2000'], 
            unit=(u.deg, u.deg)
        )
        
        idx, d2d, d3d = data_coord.match_to_catalog_sky(refcat_coord)
        
        mcat_gaia = cat.copy()
        mcat_gaia['sep'] = d2d.to(u.arcsec)
        
        mcat_cut = mcat_gaia[np.where(
            (mcat_gaia['sep'] < 3 *u.arcsec) &
            (mcat_gaia['FLAGS'] == 0)
        )]

        return mcat_cut

        
    def densityTest(self, mcat_cut) :
        
        Xmax = np.max(mcat_cut['X_IMAGE'])
        Xmin = np.min(mcat_cut['X_IMAGE'])
        Ymax = np.max(mcat_cut['Y_IMAGE'])
        Ymin = np.min(mcat_cut['Y_IMAGE'])
        
        if (
            (Xmax > 8500) &
            (Xmin < 1000) &
            (Ymax > 8500) &
            (Ymin < 1000)
        ) :
            result = True
        else :
            result = False

        if result == True : 
            dX = (Xmax - Xmin) / 8 + 1
            dY = (Ymax - Ymin) / 8 + 1
            
            N_star_sect = []
            for i in range(8) :
                for j in range(8) :
                
                    cat_sect = mcat_cut[np.where(
                        (Xmin + i*dX <= mcat_cut['X_IMAGE']) & (mcat_cut['X_IMAGE'] <= Xmin + (i+1)*dX) &
                        (Ymin + j*dY <= mcat_cut['Y_IMAGE']) & (mcat_cut['Y_IMAGE'] <= Ymin + (j+1)*dY)
                    )]
            
                    N_star_sect.append(len(cat_sect))
            
            N_expected = np.median(N_star_sect)
            N_bad_sect = len([i for i in N_star_sect if i < N_expected*0.5])
            
            if N_bad_sect > 10 :
                result = False
            else :
                result = True

        return result


    def astrometry_check_worker(self, cat_fname) :

        # Path Setting
        #cat_fname = cat_fnames[i]                             # cbftoxkmtc.20181111.005071.cat
                 # [cat_path]/cbftoxkmtc.20181111.005071.cat
        fpath_dirs = {}
        cat_id = cat_fname.split('.cat')[0]                   # cbftoxkmtc.20181111.005071
        
        #head_fname = f'{cat_id}.new.head'                     # cbftoxkmtc.20181111.005071.new.head
        data_fname = f'{cat_id}.fits.fz'                         # cbftoxkmtc.20181111.005071.fits
        seg_fname = data_fname.replace('resamp', 'segment')
        
        fpath_dirs['initial']= os.path.join(self.data_dirs['initial'], cat_fname)
        fpath_dirs['image'] = os.path.join(self.data_dirs['image'], data_fname)       # [datapath]/cbftoxkmtc.20181111.005071.new.head
        fpath_dirs['seg'] = os.path.join(self.data_dirs['image'], seg_fname)
        #data_fpath = os.path.join(self.data_dir, data_fname)       # [datapath]/cbftoxkmtc.20181111.005071.fits
        
        # Read and match data
        cat = ascii.read(fpath_dirs['initial'])
        mcat_cut = self.match_catalog(cat)
        
        #chip_checker = []
        #for chipnum in range(4) : 
        #    mcat_cut_chip = mcat_cut[np.where((mcat_cut['AMPNUM']-1) // 8 == chipnum)]
        #    if self.densityTest(mcat_cut_chip) == False :
        #        chip_checker.append(1)
        #    else :
        #        chip_checker.append(0)
        
        chip_checker = False
        #mcat_cut_chip = mcat_cut.copy()
        if len(mcat_cut) > 10 : 
            if self.densityTest(mcat_cut) == False :
                chip_checker = False
            else :
                chip_checker = True
        else :
            chip_checker = False
        
        if chip_checker == False :
        #if 1 in chip_checker :
            #ls.append(cat_fpath) #used for debug
            for key in ["initial", "image", "seg"] :
                if key == "seg" : 
                    os.system(f'mv {fpath_dirs[key]} {self.badastrom_dirs["image"]}')
                else :
                    os.system(f'mv {fpath_dirs[key]} {self.badastrom_dirs[key]}')
            #os.system(f'mv {head_fpath} {self.badastrom_data_dir}')
            #os.system(f'mv {data_fpath} {self.badastrom_data_dir}')
        #print('done') #U4D


def run(INPUTVAR) :
    ast_chckr = AstrometryChecker(INPUTVAR)
    #ast_chckr.gen_badastrom_dir()
    for cat_fname in ast_chckr.cat_fnames :
        ast_chckr.astrometry_check_worker(cat_fname)


'''if __name__ == "__main__" :
    #INPUTVAR = ['N55', '2018', 'CTIO', 'B']
    INPUTVAR = sys.argv[1:]
    main(INPUTVAR)'''


