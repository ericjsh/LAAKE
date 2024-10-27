#import os
from pipe.dataprocess_config import *
from . import *

# TODO: rename this file into utils.py for possible misleads in the future.

def find_directories(INPUTVAR: list) -> dict :
    '''
    Finds directories for associated files
    '''
    subdir = find_subdir(INPUTVAR)
    
    dir_dict = {
        'image'  : os.path.join(IMGEDIRC, subdir),
        'ampmap' : os.path.join(AMPMDIRC, subdir),
        'initial' : os.path.join(INCTDIRC, subdir),
        'dist' : os.path.join(CTLGDIRC, subdir),
        'delmag' : os.path.join(DMAGDIRC, subdir),
    }

    return dir_dict

def find_directories_backup(INPUTVAR: list) -> dict :
    '''
    Finds directories for associated files
    '''
    subdir = find_subdir(INPUTVAR)
    
    dir_dict = {
        'image'  : os.path.join(BKUPDIRC, 'laake_backup', 'image', subdir),
        'ampmap' : os.path.join(BKUPDIRC, 'laake_backup', 'ampmap', subdir),
        'initial' : os.path.join(BKUPDIRC, 'laake_backup', 'initial_cat', subdir),
        'dist' : os.path.join(BKUPDIRC, 'laake_backup', 'dist_cat', subdir),
        'delmag' : os.path.join(BKUPDIRC, 'laake_backup', 'delmag_cat', subdir),
    }

    return dir_dict

def find_subdir(INPUTVAR: list) -> str :
    field, year, tel, filt = INPUTVAR
    subdir = f'{field}_{year}/{tel}/{filt}'
    return subdir

def find_files(INPUTVAR: list) -> dict :

    dir_dict = find_directories(INPUTVAR)

    files_dict = {
        'image' : [i for i in os.listdir(dir_dict['image']) if i.endswith('.resamp.fits.fz')],
        'ampmap' : [i for i in os.listdir(dir_dict['ampmap']) if i.endswith('.ampmap.fits.fz')],
        'initial' : [i for i in os.listdir(dir_dict['initial']) if i.endswith('.resamp.cat')],
        'seg' : [i for i in os.listdir(dir_dict['image']) if i.endswith('.segment.fits.fz')]
    }

    return files_dict
