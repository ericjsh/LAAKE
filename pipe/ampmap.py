import os
NoneType = type(None)


from astropy.io import fits, ascii
import numpy as np


from pipe.dataprocess_config import *
from pipe import process_tools


import warnings
warnings.filterwarnings('ignore')


def run(INPUTVAR: list) -> None :
    #INPUTVAR = ['N55', '2018', 'CTIO', 'V']

    datadir = process_tools.find_directories(INPUTVAR)['ampmap']

    head_ls = [i for i in os.listdir(datadir) if i.endswith('.head')]

    for head_fname in head_ls : 
        #head_fname = head_ls[0]
        head_fpath = os.path.join(datadir, head_fname)
        fid = head_fname.split('.new.head')[0]

        #os.system(f'ln -s {head_fpath}')
        working_dir = os.path.join(DATFDIRC, 'swarp_files')
        os.system(f'cp {head_fpath} {working_dir}')

        ls = [np.ones((9232, 1152))*i for i in range(8)]
        pseudo_chip_img = np.block(ls)

        primary_hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([primary_hdu])
        image_hdu = fits.ImageHDU(pseudo_chip_img)
        for _ in range(4) :
            hdul.append(image_hdu)

        pchip_fname = f'{fid}.new.fits'
        hdul.writeto(pchip_fname, overwrite=True)

        os.system(f'mv {pchip_fname} {working_dir}')

        os.chdir(f'{working_dir}')

        os.system(f'swarp {pchip_fname} -c kmtnet.swarp')
        os.system('rm *.resamp.weight.fits')

        os.system('fpack -D -Y -v *.resamp.fits')

        names_ls = [(i, i.replace('resamp', 'ampmap')) for i in os.listdir('./') if i.endswith('resamp.fits.fz')]
        for i in range(4):
            prev, new = names_ls[i]
            os.system(f'mv {prev} {new}')

        os.system(f'mv *.ampmap.fits.fz {datadir}')

        os.system(f'rm {pchip_fname}')
        os.system(f'rm {head_fname}')

        os.chdir(ROOTPATH)