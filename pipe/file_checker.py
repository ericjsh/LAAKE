import os
import numpy as np
from astropy.time import Time

from pipe.dataprocess_config import *

from pipe import process_tools
from pipe import starcand
from pipe import kmtn_refcat

from shutil import copy, move


#INPUTVAR = ['N55', '2018', 'CTIO', 'V']

class CheckFiles : 
    '''
    Check whether pipeline is ready to run. If ready to run, finds the number of files to be processed.
    '''

    def __init__(self, INPUTVAR) :
        self.INPUTVAR = INPUTVAR
        self.field, self.year, self.tel, self.filt = INPUTVAR
        self.dir_dict = process_tools.find_directories(INPUTVAR)
        self.ready = False
        self.missing_files = []

    def check_ready2run(self,) -> None :  
        '''
        Check essential files for pipelines to work.
        Download from samp server : 'JD_{year}.list', 'matched_targets_{field}.csv', 'check_list_{year}.txt'
        Download via Vizier (using pipe.starcand module): 'apass_dr10_{field}.csv', 'gaia_{field}.csv'
        '''
        ready2run = True

        self.jdlist_fname = f'JD_{self.field}_{self.year}.list'
        jdlist_fpath = os.path.join(DATFDIRC, self.jdlist_fname)
        if not os.path.isfile(jdlist_fpath) :
            #print(f'{jdlist_fname} is missing!')
            ready2run = False
            self.missing_files.append(self.jdlist_fname)


        self.kmtn_trgt_fname = f'matched_targets_{self.field}.csv'
        kmtn_trgt_fpath = os.path.join(DATFDIRC, self.kmtn_trgt_fname)
        if not os.path.isfile(kmtn_trgt_fpath) :
            #print(f'{kmtn_trgt_fname} is missing!')
            ready2run = False
            self.missing_files.append(self.kmtn_trgt_fname)


        #checklist_fname = f'check_list_{self.year}.txt'
        #checklist_fpath = os.path.join(DATFDIRC, checklist_fname)
        #if not os.path.isfile(checklist_fpath) :
        #    #print(f'{checklist_fname} is missing!')
        #    ready2run = False
        #    self.missing_files.append(checklist_fname)


        apass_fname = f'apass_dr10_{self.field}.csv'
        apass_fpath = os.path.join(DATFDIRC, apass_fname)
        if not os.path.isfile(apass_fpath) :
            #print(f'{apass_fname} is missing!')
            ready2run = False
            self.missing_files.append(apass_fname)

        gaia_cat_fname = f'gaia_{self.field}.csv'
        gaia_cat_fpath = os.path.join(DATFDIRC, gaia_cat_fname)
        if not os.path.isfile(gaia_cat_fpath) :
            #os.system(f'python starcand.py {field} {kmtn_trgt_fname}')
            ready2run = False
            self.missing_files.append(gaia_cat_fname)

        if self.field == 'N7793' :
            self.check_list_fname = f'check_list_{self.field}_{self.year}.txt'
            check_list_fpath = os.path.join(DATFDIRC, self.check_list_fname)
            if not os.path.isfile(check_list_fpath) :
                ready2run = False
                self.missing_files.append(self.check_list_fname)
        else :
            self.check_list_fname = 'placeholder'


        self.ready = ready2run

    
    def download_files(self) -> None :
        '''
        Download missing files (using pipe.starcand module)
        '''
        missing_files_string = ', '.join(self.missing_files)
        print(f'Not ready to run. Followig files are missing: {missing_files_string}.')
        print('Downloading essential files ...')

        if self.kmtn_trgt_fname in self.missing_files :
            kmtn_refcat.kmtn_refcat_gen(self.INPUTVAR)

        if (self.field == 'N7793') & (self.check_list_fname in self.missing_files) :
            check_list_fname_old = 'check_list.txt'
            check_list_fpath_old = os.path.join(BKUPDIRC, f'{self.field}_{self.year}', check_list_fname_old)
            check_list_fpath_old2new = os.path.join(DATFDIRC, check_list_fname_old)
            check_list_fpath_new = os.path.join(DATFDIRC, self.check_list_fname)
            copy(check_list_fpath_old, check_list_fpath_old2new)
            move(check_list_fpath_old2new, check_list_fpath_new)

        if self.jdlist_fname in self.missing_files :
            jdlist_fname_old = 'JD.list'
            jdlist_fpath_old = os.path.join(BKUPDIRC, f'{self.field}_{self.year}', jdlist_fname_old)
            if os.path.isfile(jdlist_fpath_old) :
                jdlist_fpath_old2new = os.path.join(DATFDIRC, jdlist_fname_old)
                jdlist_fpath_new = os.path.join(DATFDIRC, self.jdlist_fname)
                copy(jdlist_fpath_old, jdlist_fpath_old2new)
                move(jdlist_fpath_old2new, jdlist_fpath_new)
            else :
                self.calc_jd()




        #starcand.starcat_generator([self.field, self.kmtn_trgt_fname])
        starcand.starcat_generator(self.INPUTVAR)
        print('Download complete!')
        
    def calc_jd(self) :
        # ------------------ Calc. JD/MJD -------------------------------
        print("### Make JD/MJD list ################################")

        check_list_fpath = os.path.join(DATFDIRC, self.check_list_fname)
        jdlist_fpath = os.path.join(DATFDIRC, self.jdlist_fname)

        dtype2 = np.dtype([('file', 'U26'), ('exptime', 'int'), ('utc', 'U19')])
        file, exptime, utc = np.loadtxt(check_list_fpath, usecols=[0, 3, 6], dtype=dtype2, unpack=True)

        t = Time(utc, scale='utc', format='isot')
        jd = t.jd + (exptime/2.)/3600./24.
        mjd = t.mjd + (exptime/2.)/3600./24.

        dtype3 = np.dtype([('var1', 'U26'), ('var2', float), ('var3', float)])
        outdata = np.zeros(jd.size, dtype=dtype3)
        outdata['var1'] = file
        outdata['var2'] = jd
        outdata['var3'] = mjd

        np.savetxt(jdlist_fpath, outdata, delimiter='  ', fmt=['%26s', '%.9f', '%.11f'], header="file, mid-JD, mid-MJD")

    def count_file_nums(self) -> None :
        '''
        Count number of image files ('*.fitz.fz') and external header files ('*.head'). 
        '''

        self.img_filenames = [i for i in os.listdir(self.dir_dict['image']) if i.endswith('.fits.fz')]
        self.head_filenames = [i for i in os.listdir(self.dir_dict['ampmap']) if i.endswith('.head')]

        N_img = len(self.img_filenames)
        N_head = len(self.head_filenames)

        if N_img == 8*N_head :
            print(f'Total number of .fits images to be processed : {N_img}')
        else :
            print('Some files are missing')
            self.move_files_badastrom()

    
    def move_files_badastrom(self) -> None :
        '''
        Check files that are either missing image for external header to check the final number of data to be processed.
        '''

        image_set = set([fname.split('.new')[0] for fname in self.img_filenames])
        ampmap_set = set([fname.split('.new')[0] for fname in self.head_filenames])
        ampmap_missing_fid = list(image_set - (image_set & ampmap_set))

        ampmap_missing = []
        for idx in range(len(ampmap_missing_fid)) : 
            ampmap_missing += [f for f in self.img_filenames if f.startswith(ampmap_missing_fid[idx])]

        directory = self.dir_dict['image']
        badastrom_path = os.path.join(directory, 'BADASTROM')
        os.makedirs(badastrom_path, exist_ok=True)

        for filename in ampmap_missing :
            filepath = os.path.join(directory, filename)
            os.system(f'mv {filepath} {badastrom_path}')
            

def run(INPUTVAR: list) -> None :
    '''
    Check whether pipeline is ready to run. If ready to run, finds the number of files to be processed.
    '''
    checkfiles = CheckFiles(INPUTVAR)
    checkfiles.check_ready2run()

    for _ in range(2) :
        if checkfiles.ready :
            print('All essential datafiles exist.')
            checkfiles.count_file_nums()
            break
        else :
            checkfiles.download_files()

    print('Ready to go!')

    

    # TODO: test this module
