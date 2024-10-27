import os
import numpy as np
from astropy.io import fits


'''total_file = [i for i in os.listdir('.') if i.endswith('.fits.fz')]
for filename in total_file :
    cat_filename = filename.replace('.fits.fz', '.cat')
    seg_filename = filename.replace('resamp.fits.fz', 'segment.fits')
    apertures = np.arange(2.5, 75.1, 2.5)
    apertures_cmd = ','.join([str(ap) for ap in apertures])
    os.system(
        f'sex {filename} -c kmtnet_phot.sex ' + 
        f'-CATALOG_NAME {cat_filename} ' + 
        '-BACK_SIZE 32 ' + 
        f'-PHOT_APERTURES {apertures_cmd} ' +
        '-CHECKIMAGE_TYPE SEGMENTATION ' + 
        f'-CHECKIMAGE_NAME {seg_filename}'
    )
    os.system(f'fpack -D -Y -v {seg_filename}')'''

total_file = [i for i in os.listdir('.') if i.endswith('.fits.fz')]
for filename in total_file :
    cat_filename = filename.replace('.fits.fz', '.cat')
    seg_filename = filename.replace('resamp.fits.fz', 'segment.fits')
    os.system(
        f'sex {filename} -c kmtnet_phot.sex ' + 
        f'-CATALOG_NAME {cat_filename} ' + 
        '-BACK_SIZE 32 ' + 
        '-CHECKIMAGE_TYPE SEGMENTATION ' + 
        f'-CHECKIMAGE_NAME {seg_filename}'
    )
    os.system(f'fpack -D -Y -v {seg_filename}')