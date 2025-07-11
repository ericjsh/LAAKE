import os

import numpy as np
import pandas as pd

from astropy.io import ascii
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table, vstack
from astroquery.vizier import Vizier

sep_lim = 2*u.arcsec

from pipe.dataprocess_config import *

from pipe import process_tools


#INPUTVAR = ['N59', '2016', 'CTIO', 'B']

def kmtn_refcat_gen(INPUTVAR, vanilla=False) :
    field, *_ = INPUTVAR

    if vanilla :
        datadir = '/Users/kamp/Desktop/KAMP/LAAKE/vanilla_workspace'
    else :
        data_dirs = process_tools.find_directories(INPUTVAR)
        datadir = data_dirs['initial']

    mef_fsizes = [[f, os.path.getsize(os.path.join(datadir, f)) /1024**2] for f in os.listdir(datadir)]

    # File size sigma-clipping
    df = pd.DataFrame(mef_fsizes, columns=['filename', 'filesize'])

    MED = np.median(df['filesize'])
    STD = np.std(df['filesize'])

    df_cut = df[(MED-STD < df['filesize']) & (df['filesize'] < MED + STD)]

    datadir_fsize_cut = list(df_cut['filename'])

    # Checking all 4 chips exist
    mef_names = [f.split('.new')[0] for f in datadir_fsize_cut]
    mef_idx = 0

    while True :
        mef_name = mef_names[mef_idx]
        
        chip_names = [f for f in datadir_fsize_cut if f.startswith(mef_name) and f.endswith('resamp.cat')]
        
        if len(chip_names) == 4 :
            break
        
        else :
            mef_idx += 1

    chip_fpaths = [os.path.join(datadir, chipname) for chipname in chip_names]
    print(chip_fpaths)
    kmtn_refcat_raw = vstack([ascii.read(chipfpath) for chipfpath in chip_fpaths])

    trgt_radec = SkyCoord(
        ra=kmtn_refcat_raw['ALPHA_J2000'], 
        dec=kmtn_refcat_raw['DELTA_J2000']
    )

    # Download Gaia Catalog
    Vizier.ROW_LIMIT = -1

    cent_ra = np.median(trgt_radec.ra)
    cent_dec = np.median(trgt_radec.dec)

    gaia = Vizier.query_region(
        coord.SkyCoord(
            ra=cent_ra, dec=cent_dec, 
            unit=(u.deg, u.deg), frame='icrs'
        ),
        radius=1.5*u.deg,
        catalog=['Gaia'],
    )

    gaia_cat = gaia[0]

    gaia_radec = SkyCoord(
        ra=gaia_cat['RAJ2000'],
        dec=gaia_cat['DEJ2000']
    )

    # Matching KMTNet and Gaia catalogues
    idx, d2d, d3d = trgt_radec.match_to_catalog_sky(gaia_radec)


    # Creating matched calatog
    mcat = kmtn_refcat_raw.copy()

    mcat['ALPHA_J2000'].name = 'ra'
    mcat['DELTA_J2000'].name = 'dec'

    mcat['sep'] = d2d.to(u.arcsec)

    mcat_cut = mcat[np.where(mcat['sep']<sep_lim)]
    mcat_cut['Jname'] = np.array(['a'*21]*len(mcat_cut))

    for idx in range(len(mcat_cut)) :
        jname = 'J' + ''.join(['{:.7f}'.format(mcat_cut['ra'][idx]), '{:.7f}'.format(mcat_cut['dec'][idx])])
        mcat_cut['Jname'][idx] = jname

    matched_targets_file = mcat_cut[['ra', 'dec', 'Jname', 'sep']]
    matched_targets_filename = f'matched_targets_{field}.csv'
    matched_targets_filepath = os.path.join(DATFDIRC, matched_targets_filename)


    ascii.write(matched_targets_file, matched_targets_filepath, overwrite=True) #write to txt