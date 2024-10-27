import os
NoneType = type(None)

from astropy.table import Table, vstack
from astropy.io import fits, ascii
from astropy.stats import SigmaClip, sigma_clip
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, Angle, match_coordinates_sky
import numpy as np

import sys
sys.setrecursionlimit(10**7)

import pickle

#from tqdm.notebook import tqdm
from tqdm import tqdm

from multiprocessing import Pool

import warnings
warnings.filterwarnings('ignore')

import sys
sys.setrecursionlimit(10**7)


from pipe import process_tools
from pipe.dataprocess_config import *



def isfile_all(ls) :
    summ = 0
    for item in ls :
        summ += os.path.isfile(item)

    return summ == len(ls)


def grp_nums(seg_sect) :

    field = seg_sect.copy()
    
    def flood_fill(x ,y, old, new):
        # we need the x and y of the start position, the old value,
        # and the new value
        # the flood fill has 4 parts
        # firstly, make sure the x and y are inbounds
        if x < 0 or x >= len(field[0]) or y < 0 or y >= len(field):
            return
        # secondly, check if the current position equals the old value
        if field[y][x] != old:
            return
    
        # thirdly, set the current position to the new value
        field[y][x] = new
        # fourthly, attempt to fill the neighboring positions
        flood_fill(x+1, y, old, new)
        flood_fill(x-1, y, old, new)
        flood_fill(x, y+1, old, new)
        flood_fill(x, y-1, old, new)

    grp_num = 0
    while (field==1).sum() > 0 :
        grp_num += 1
        ys,xs = np.where(field==1)
        flood_fill(xs[0], ys[0], 1, 0)

    return grp_num


def run(INPUTVAR: list) -> None :
    #INPUTVAR = ['N55', '2018', 'CTIO', 'V']
    files_dict = process_tools.find_files(INPUTVAR)
    data_dirs = process_tools.find_directories(INPUTVAR)

    for cat_id in range(len(files_dict['initial'])) :
        cat_fname = files_dict['initial'][cat_id]
        chip_id = int(cat_fname.split('.resamp.cat')[0].split('.new.')[1])
        print(f'Running for {cat_fname} ...')
        
        fid_chip = cat_fname.split('.resamp.cat')[0]
        fid = fid_chip.split('new')[0]

        seg_fnames = [i for i in files_dict['seg'] if i.startswith(fid_chip)]
        ampmap_fnames = [i for i in files_dict['ampmap'] if i.startswith(fid_chip)]

        if (len(seg_fnames) == 1) and (len(ampmap_fnames) == 1) :

            seg_fname = seg_fnames[0]
            ampmap_fname = ampmap_fnames[0]

            cat_fpath = os.path.join(data_dirs['initial'], cat_fname)
            seg_fpath = os.path.join(data_dirs['image'], seg_fname)
            ampmap_fpath = os.path.join(data_dirs['ampmap'], ampmap_fname)

            #if isfile_all([cat_fpath, seg_fpath, ampmap_fpath]) : 

            cat = ascii.read(cat_fpath)

            if len(cat) > 1 :
                new_stats_keys = [
                    'ampnum',       #amp position of object
                    'POSA',         #Percentage Of (object placed at) Same Amp
                    'DTAE',         #Distance To Amp Edge
                    'grpnum',       #number of semi-groups within the object
                    'AAPG'          #Average Area Per Group
                ]

                for key in new_stats_keys :
                    cat[key] = np.array([0.])*len(cat)


                seg_fitsfile = fits.open(seg_fpath)
                seg_img_data = seg_fitsfile[1].data

                X_max_data = len(seg_img_data[0])
                Y_max_data = len(seg_img_data)

                ampmap_fitsfile = fits.open(ampmap_fpath)
                ampmap_img_data = ampmap_fitsfile[1].data
                ampmap_img_data = np.around(ampmap_img_data).astype(int)

                cat_cut = cat[np.where(
                    (cat['MAG_AUTO'] < 50) & 
                    (cat['FLAGS'] < 10)
                )]


                #cat_obj = cat_cut[np.where(cat_cut['NUMBER'] == 7221)]
                #objnum = cat_obj['NUMBER']

                for n in range(len(cat_cut)) : 
                    #print(n)
                    #n=0
                    cat_obj = cat_cut[n]

                    X = int(cat_obj['X_IMAGE'])
                    Y = int(cat_obj['Y_IMAGE'])

                    ampnum = ampmap_img_data[Y,X]

                    sect_dim = 100

                    Y_min = max(Y-sect_dim, 0)
                    Y_max = min(Y+sect_dim, Y_max_data)
                    X_min = max(X-sect_dim, 0)
                    X_max = min(X+sect_dim, X_max_data)

                    seg_sect = (seg_img_data[Y_min:Y_max, X_min:X_max] == cat_obj['NUMBER'])
                    ampmap_sect = (ampmap_img_data[Y_min:Y_max, X_min:X_max] == ampnum)


                    if Y_min == 0 :
                        pad = np.ones((sect_dim-Y, len(seg_sect[0])))*(-1)
                        seg_sect = np.block([[pad], [seg_sect]])
                        ampmap_sect = np.block([[pad], [ampmap_sect]])

                    if Y_max == len(seg_img_data) :
                        pad = np.ones((Y - Y_max_data + sect_dim, len(seg_sect[0])))*(-1)
                        seg_sect = np.block([[seg_sect], [pad]])
                        ampmap_sect = np.block([[ampmap_sect], [pad]])

                    if X_min == 0 :
                        pad = np.ones((sect_dim*2, sect_dim-X)) * (-1)
                        seg_sect = np.block([pad, seg_sect])
                        ampmap_sect = np.block([pad, ampmap_sect])

                    if X_max == X_max_data :
                        pad = np.ones((sect_dim*2, X - X_max_data + sect_dim)) * (-1)
                        seg_sect = np.block([seg_sect, pad])
                        ampmap_sect = np.block([ampmap_sect, pad])

                    same_amp = ampmap_sect * seg_sect
                    same_amp_per = (same_amp == 1).sum() / (seg_sect == 1).sum() * 100

                    ampmap_dist = (ampmap_sect == 0).sum() / 200**2

                    if (cat_obj['MAG_AUTO'] > 50) | (cat_obj['FLAGS'] > 10) :
                        grpnum = 99
                        avg_area = 0
                    else :
                        #grpnum = grp_nums(seg_sect)
                        grpnum=1
                        avg_area = (seg_sect==1).sum() / grpnum

                    cat_cut[n]['ampnum'] = int(round(ampnum)) + (chip_id - 1) * 8
                    cat_cut[n]['POSA'] = same_amp_per
                    cat_cut[n]['DTAE'] = ampmap_dist
                    cat_cut[n]['grpnum'] = grpnum
                    cat_cut[n]['AAPG'] = avg_area


                cat_dist_cut = cat_cut[np.where(
                    (cat_cut['DTAE'] == 0) &
                    (cat_cut['ISOAREA_IMAGE'] > 8)
                )]

                datadir_new = data_dirs['dist']
                os.makedirs(datadir_new, exist_ok=True)

                cat_fpath_new = os.path.join(datadir_new, cat_fname)
                ascii.write(cat_dist_cut, cat_fpath_new, overwrite=True)
                
                N_src = len(cat_dist_cut)
                per_survived = round(N_src/len(cat_cut) * 100,2)


                print(f'{per_survived}% ({N_src}) sources survived')
                print(f'file saved at {cat_fpath_new}')
