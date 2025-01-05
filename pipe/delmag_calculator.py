import os

#from types import NoneType
NoneType = type(None)

from astropy.table import Table, vstack
from astropy.io import fits, ascii
from astropy.stats import SigmaClip, sigma_clip
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, Angle, match_coordinates_sky
import numpy as np

from scipy.interpolate import interp1d
from statsmodels.tsa.stattools import adfuller

import pickle

#from tqdm.notebook import tqdm
from tqdm import tqdm

from multiprocessing import Pool

import warnings
warnings.filterwarnings('ignore')

from pipe.dataprocess_config import *

from pipe import process_tools



class DelmagCalculator :

    def __init__(self, data_info) :
        field, year, tel, filt, filename = data_info
        self.tel = tel
        self.filt = filt
        self.field = field
        self.year = year
        self.filename = filename
        #self.chipnum = self.filename.split('.')[4]
        self.datadirc = os.path.join(CTLGDIRC, f'{field}_{year}', tel, filt)
        self.IDcut = 'none'


    def aps_stcat_gen(self, aps_only=False) :
        '''
        Summary: Call APASS reference catalgue for N55 field
        
        INPUT  : none
        OUTPUT : APASS reference catalogue
        '''
        
        #sID_bad_vis_path = os.path.join(DATFDIRC, 'bad_star_N55_2018.csv')
        #sID_bad_vis = [i[0] for i in ascii.read(sID_bad_vis_path)]

        if aps_only == True :
            aps_stcat_path = os.path.join(DATFDIRC, f'apass_dr10_{self.field}_good.txt')
            
            dtype2 = np.dtype([ 
                ('sID', 'U5'), 
                ('ra', float), 
                ('dec', float), 
                ('Bref', float), 
                ('Vref', float), 
                ('Iref', float), 
                ('e_Bref', float), 
                ('e_Vref', float), 
                ('e_Iref', float) 
            ])
            
            aps_stcat_tab = Table(
                np.genfromtxt(
                    aps_stcat_path, 
                    comments='#', 
                    dtype=dtype2, unpack=True
                ),
                dtype=dtype2,
            )



            if self.IDcut == 'none' :
                result = aps_stcat_tab
                
            else :
                result = Table(dtype=dtype2)

                
                #IDcut_final = [i for i in IDcut if i not in sID_bad_vis]
                for sID in self.IDcut :
                    result.add_row(
                        aps_stcat_tab[np.where(aps_stcat_tab['sID'] == sID)][0]
                    )
        else :
            aps_stcat_path = os.path.join(DATFDIRC, f'starcand_{self.field}.csv')
            aps_stcat_tab = ascii.read(aps_stcat_path)
            
            if self.IDcut == 'none' :
                #result = aps_stcat_tab[np.where(np.isin(aps_stcat_tab['sID'], sID_bad_vis) == False)]
                result = aps_stcat_tab
                
            else :
                #IDcut_final = [i for i in self.IDcut if i not in sID_bad_vis]
                IDcut_final = [i for i in self.IDcut]
                result = aps_stcat_tab[np.where(np.isin(aps_stcat_tab['sID'], IDcut_final))]
                with open('stars_check', 'wb') as fp :
                    pickle.dump(IDcut_final, fp)
                #result = aps_stcat_tab[np.where(np.isin(aps_stcat_tab['sID'], self.IDcut))]
            
        return result


    def aps_mag_names(self):
        filter = self.filt
        kmtn_band_ls = ['B', 'V', 'I']
        aps_mag_name_ls  = ['Bref', 'V_gaia', 'I_gaia']
        aps_emag_name_ls = ['e_'+i for i in aps_mag_name_ls]
        
        band_idx = kmtn_band_ls.index(filter)
        
        aps_mag_name = aps_mag_name_ls[band_idx]
        aps_emag_name = aps_emag_name_ls[band_idx]

        return aps_mag_name, aps_emag_name


    def refcat_gen(self) :
        '''
        Summary : Select required band data; apply magnitude cut 13.5 < m < 17

        INPUT   : filter
        OUTPUT  : APASS reference catalog with required band data
        '''
        
        aps_stcat_tab = self.aps_stcat_gen()
        aps_mag_name, aps_emag_name = self.aps_mag_names()
    
        refcat_aps = aps_stcat_tab['sID', 'ra', 'dec', aps_mag_name, aps_emag_name]
        refcat_aps_cut = refcat_aps[np.where(
            (13.5 < refcat_aps[aps_mag_name]) & (refcat_aps[aps_mag_name] < 17)
        )]
    
        return refcat_aps_cut


    def find_mjd(self) :
        '''
        Summary : Find corresponding MJD
        INPUT   : filename
        OUTPUT  : corresponding MJD
        '''
        year = self.year
        filename = self.filename
    
        jdlist_path = os.path.join(DATFDIRC, f'JD_{self.field}_{year}.list')
        
        dtype1 = np.dtype([
            ('org_name', 'U26'), 
            ('mjd_num', 'U17')
        ])
        org_name, mjd_num = np.genfromtxt(
            jdlist_path, 
            usecols=[0, 2], 
            comments='#', 
            dtype=dtype1, unpack=True
        ) 
        
        ### MJD info 
        
        obsid0 = org_name.copy()
        for njd in range(len(org_name)):
            obsdate1 = org_name[njd].split('.')[1]
            obsdate2 = org_name[njd].split('.')[2]
            obsid0[njd] = obsdate1+'.'+obsdate2  # '20160116.044712'
        
        obsid = filename.split('.')[1] + '.' + filename.split('.')[2]
        mjd = mjd_num[np.where(obsid == obsid0)][0]
    
        return float(mjd)


    def gen_mcat(self, refcat, cat_name) :
        '''
        Summary : Match KMTN data and APASS reference catalogue. Then clean the matched catalogue.
        INPUT   : filename - KMTN data name ; refcat - APASS reference catalogue
        OUTPUT  : mcat_cut - Matched catalogue
        '''
        #filename = self.filename
        
        sep_lim = 2. * u.arcsec
    
        filepath = os.path.join(self.datadirc, self.filename)
        aps_mag_name, aps_emag_name = self.aps_mag_names()
        
        #mjd = find_mjd(filename)
        data = ascii.read(filepath)
        
        data_coord = SkyCoord(
            ra=data['ALPHA_J2000'],
            dec=data['DELTA_J2000'], 
            unit=(u.deg, u.deg)
        )
        
        refcat_coord = SkyCoord(
            ra=refcat['ra'], 
            dec=refcat['dec'], 
            unit=(u.deg, u.deg)
        )
        
        idx, d2d, d3d = data_coord.match_to_catalog_sky(refcat_coord)
        
        mcat = data.copy()
        mcat['sep'] = d2d.to(u.arcsec)
    
        if cat_name == 'APASS' :
            mcat['ref_m'] = refcat[idx][aps_mag_name]
            mcat['ref_merr'] = refcat[idx][aps_emag_name]
            mcat['sID']   = refcat[idx]['sID']
            
            mcat_cut = mcat[np.where(mcat['sep'] < sep_lim)]
    
            fwhm_mn = np.median(mcat_cut['FWHM_IMAGE'])
            fwhm_std = np.std(mcat_cut['FWHM_IMAGE'])
            
            result = mcat_cut[np.where(
                (fwhm_mn-3*fwhm_std < mcat_cut['FWHM_IMAGE']) & (mcat_cut['FWHM_IMAGE'] < fwhm_mn+3*fwhm_std) &
                (mcat_cut['FLAGS'] == 0) &
                (mcat_cut['ELLIPTICITY'] < 0.1) &
                (mcat_cut['MAG_AUTO'] < 99)
            )]
        elif cat_name == 'GAIA' :
            mcat_cut = mcat[np.where(mcat['sep'] < sep_lim)]
    
            result = mcat_cut
    
        '''for nf in range(len(result)) :
            result_row = result[nf]
            ampnum_prev = result_row['ampnum']
            ampnum_new = (float((self.filename.split('new.')[1]).split('.resamp')[0])-1)*8 + ampnum_prev
            result_row['ampnum'] = ampnum_new'''
        return result


    def delmag_engine_2(self, mcat) : 
        '''
        Summary : Calculate delmag. Then clean catalogue using delmag.
    
        INPUT   : mcat - matched catalogue.
        OUTPUT  : result - catalogue with delmag and delmagerr.
        '''
        
        mcat['offset'] = mcat['ref_m'] - mcat['MAG_AUTO']
        
        dmag_mn = np.median(mcat['offset'])
        dmag_std = np.std(mcat['offset'])
    
        mcat_cut = mcat[np.where(
            (dmag_mn-dmag_std < mcat['offset']) &
            (mcat['offset'] < dmag_mn+dmag_std)
        )]
    
        if len(mcat_cut) > 3 :
            dmag_mn = np.median(mcat_cut['offset'])
            dmag_std = np.std(mcat_cut['offset'])
    
            mcat_cut['delmag'] = dmag_mn
            mcat_cut['dmag_err'] = dmag_std
    
            result_sub = mcat_cut
    
        else : 
            mcat['delmag'] = dmag_mn
            mcat['dmag_err'] = dmag_std
    
            result_sub = mcat

        aps_mag_name, aps_emag_name = self.aps_mag_names()

        result_sub['mjd'] = self.find_mjd()
        result_sub['mag_inst'] = result_sub['MAG_AUTO'] + result_sub['delmag']
        result_sub['magerr_inst'] = np.sqrt(result_sub['MAGERR_AUTO']**2 + result_sub['dmag_err']**2)
        #result = result_sub['ALPHA_J2000', 'DELTA_J2000', 'X_IMAGE', 'Y_IMAGE', 'mag_inst', 'magerr_inst', 'MAG_AUTO', 'MAGERR_AUTO', 'MAG_ISO', 'MAGERR_ISO', 'sID', 'mjd', 'delmag', 'dmag_err', 'ref_m', 'ref_merr', 'chip']
        result = result_sub.copy()
        #result['ALPHA_J2000'].name = 'ra'
        #result['DELTA_J2000'].name = 'dec'
        result['filename'] = self.filename
    
        return result
    
    def delmag_engine(self, mcat) : 
        '''
        Summary : Calculate delmag. Then clean catalogue using delmag.

        INPUT   : mcat - matched catalogue.
        OUTPUT  : result - catalogue with delmag and delmagerr.
        '''
        ls_result = []
        if self.filt == 'B' : 
            std_lim = 0.06
        else :
            std_lim = 0.03

        for ampnum in range(0,32) : 
            mcat_amp = mcat[np.where(mcat['ampnum'] == ampnum)]
            if len(mcat_amp) >= 5 :
        
                mcat_amp['offset'] = mcat_amp['ref_m'] - mcat_amp['MAG_AUTO']
                #mcat_amp['seeing'] = np.median(mcat_amp['FWHM_IMAGE'])
                
                seeing = np.median(mcat_amp['FWHM_IMAGE'])
                dmag_mn = np.median(mcat_amp['offset'])
                dmag_std = np.std(mcat_amp['offset'])
                
                mcat_cut = mcat_amp[np.where(
                    (dmag_mn-dmag_std < mcat_amp['offset']) &
                    (mcat_amp['offset'] < dmag_mn+dmag_std)
                )]
                
                if len(mcat_cut) >= 3 :
                    dmag_mn_cut = np.median(mcat_cut['offset'])
                    dmag_std_cut = np.std(mcat_cut['offset'])
                    seeing_mn_cut = np.median(mcat_cut['FWHM_IMAGE'])

                    if dmag_std_cut < std_lim :
                
                        mcat_cut['delmag'] = dmag_mn_cut
                        mcat_cut['dmag_err'] = dmag_std_cut
                        mcat_cut['ampnum'] = ampnum
                        mcat_cut['seeing'] = seeing_mn_cut
        
                        ls_result.append(mcat_cut)
        
                else : 
                    if dmag_std < std_lim : 
                        mcat_amp['delmag'] = dmag_mn
                        mcat_amp['dmag_err'] = dmag_std
                        mcat_amp['ampnum'] = ampnum
                        mcat_amp['seeing'] = seeing
        
                        ls_result.append(mcat_amp)

        if len(ls_result) > 0 :

            result_sub = vstack(ls_result)
            
            result_sub['mjd'] = self.find_mjd()
            result_sub['mag_inst'] = result_sub['MAG_AUTO'] + result_sub['delmag']
            result_sub['magerr_inst'] = np.sqrt(result_sub['MAGERR_AUTO']**2 + result_sub['dmag_err']**2)
            #result = result_sub['ALPHA_J2000', 'DELTA_J2000', 'X_IMAGE', 'Y_IMAGE', 'mag_inst', 'magerr_inst', 'MAG_AUTO', 'MAGERR_AUTO', 'sID', 'mjd', 'delmag', 'dmag_err', 'ref_m', 'ref_merr', 'AMPNUM']

            result = result_sub.copy()
            #result['ALPHA_J2000'].name = 'ra'
            #result['DELTA_J2000'].name = 'dec'
        else :
            result = self.result_init()

        return result


    def result_init(self) :

        #dtype_result = np.dtype([('NUMBER', '<i8'), ('ALPHA_J2000', '<f8'), ('DELTA_J2000', '<f8'), ('X_IMAGE', '<f8'), ('Y_IMAGE', '<f8'), ('MAG_AUTO', '<f8'), ('MAGERR_AUTO', '<f8'), ('MAG_ISO', '<f8'), ('MAGERR_ISO', '<f8'), ('MAG_ISOCOR', '<f8'), ('MAGERR_ISOCOR', '<f8'), ('MAG_PETRO', '<f8'), ('MAGERR_PETRO', '<f8'), ('MAG_BEST', '<f8'), ('MAGERR_BEST', '<f8'), ('MAG_APER', '<f8'), ('MAG_APER_1', '<f8'), ('MAG_APER_2', '<f8'), ('MAG_APER_3', '<f8'), ('MAG_APER_4', '<f8'), ('MAG_APER_5', '<f8'), ('MAG_APER_6', '<f8'), ('MAG_APER_7', '<f8'), ('MAG_APER_8', '<f8'), ('MAG_APER_9', '<f8'), ('MAG_APER_10', '<f8'), ('MAG_APER_11', '<f8'), ('MAG_APER_12', '<f8'), ('MAG_APER_13', '<f8'), ('MAG_APER_14', '<f8'), ('MAG_APER_15', '<f8'), ('MAG_APER_16', '<f8'), ('MAG_APER_17', '<f8'), ('MAG_APER_18', '<f8'), ('MAG_APER_19', '<f8'), ('MAG_APER_20', '<f8'), ('MAG_APER_21', '<f8'), ('MAG_APER_22', '<f8'), ('MAG_APER_23', '<f8'), ('MAG_APER_24', '<f8'), ('MAG_APER_25', '<f8'), ('MAG_APER_26', '<f8'), ('MAG_APER_27', '<f8'), ('MAG_APER_28', '<f8'), ('MAG_APER_29', '<f8'), ('MAGERR_APER', '<f8'), ('MAGERR_APER_1', '<f8'), ('MAGERR_APER_2', '<f8'), ('MAGERR_APER_3', '<f8'), ('MAGERR_APER_4', '<f8'), ('MAGERR_APER_5', '<f8'), ('MAGERR_APER_6', '<f8'), ('MAGERR_APER_7', '<f8'), ('MAGERR_APER_8', '<f8'), ('MAGERR_APER_9', '<f8'), ('MAGERR_APER_10', '<f8'), ('MAGERR_APER_11', '<f8'), ('MAGERR_APER_12', '<f8'), ('MAGERR_APER_13', '<f8'), ('MAGERR_APER_14', '<f8'), ('MAGERR_APER_15', '<f8'), ('MAGERR_APER_16', '<f8'), ('MAGERR_APER_17', '<f8'), ('MAGERR_APER_18', '<f8'), ('MAGERR_APER_19', '<f8'), ('MAGERR_APER_20', '<f8'), ('MAGERR_APER_21', '<f8'), ('MAGERR_APER_22', '<f8'), ('MAGERR_APER_23', '<f8'), ('MAGERR_APER_24', '<f8'), ('MAGERR_APER_25', '<f8'), ('MAGERR_APER_26', '<f8'), ('MAGERR_APER_27', '<f8'), ('MAGERR_APER_28', '<f8'), ('MAGERR_APER_29', '<f8'), ('FLAGS', '<i8'), ('CLASS_STAR', '<f8'), ('ELLIPTICITY', '<f8'), ('ISOAREA_IMAGE', '<i8'), ('FWHM_IMAGE', '<f8'), ('KRON_RADIUS', '<f8'), ('ampnum', '<f8'), ('POSA', '<f8'), ('DTAE', '<f8'), ('grpnum', '<f8'), ('AAPG', '<f8'), ('sep', '<f8'), ('ref_m', '<f8'), ('ref_merr', '<f8'), ('sID', '<U5'), ('offset', '<f8'), ('delmag', '<f8'), ('dmag_err', '<f8')])
        
        dtype_result = np.dtype([('NUMBER', '<i8'), ('ALPHA_J2000', '<f8'), ('DELTA_J2000', '<f8'), ('X_IMAGE', '<f8'), ('Y_IMAGE', '<f8'), ('MAG_AUTO', '<f8'), ('MAGERR_AUTO', '<f8'), ('FLAGS', '<i8'), ('CLASS_STAR', '<f8'), ('ELLIPTICITY', '<f8'), ('ISOAREA_IMAGE', '<i8'), ('FWHM_IMAGE', '<f8'), ('KRON_RADIUS', '<f8'), ('ampnum', '<f8'), ('POSA', '<f8'), ('DTAE', '<f8'), ('grpnum', '<f8'), ('AAPG', '<f8'), ('sep', '<f8'), ('ref_m', '<f8'), ('ref_merr', '<f8'), ('sID', '<U5'), ('offset', '<f8'), ('delmag', '<f8'), ('dmag_err', '<f8')])

        result = Table(dtype=dtype_result)
    
        return result


    def worker(self, mp_var) :
        filename, refcat = mp_var
        self.filename = filename
        #self.chipnum = filename.split('.')[4]
        #print(self.filename)
        mcat = self.gen_mcat(refcat, cat_name='APASS')
        if len(mcat) > 50 :
            cat = self.delmag_engine(mcat)
            return cat


    def calc_delmag_all(self) :
    
        refcat = self.refcat_gen()
        filenames = [i for i in os.listdir(self.datadirc) if i.endswith('.cat')]
        mp_var = [(filename, refcat) for filename in filenames]
        
        result = self.result_init()
    
        #CORENUM = 10
    
        tqdm_bar_fmt = '{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}'

        ls_result = []
    
        with Pool(CORENUM) as p :
            for output in tqdm(p.imap(self.worker, mp_var), total=len(filenames), bar_format=tqdm_bar_fmt) :
                if type(output) != NoneType :
                    ls_result.append(output)

        '''for i in tqdm(range(len(mp_var))) :
            output = self.worker(mp_var[i])
            if type(output) != NoneType :
                for output_i in output :
                    result.add_row(output_i)'''
    
        result = vstack(ls_result)

        return result


    def reject_star2(self, result, ensemble=False) :
        sID_ls = list(set(result['sID']))
        
        dtype_lc = np.dtype([
            ('sID', 'U5'),
            ('med_mag', float),
            ('RMS', float),
            ('N_ep', int),
            ('med_err', float),
            ('fvar', float)
        ])
        lc_result = Table(dtype=dtype_lc)
        
        for i in range(len(sID_ls)) :
            sID_i = sID_ls[i]
            sID_lc = result[np.where(
                (result['sID'] == sID_i) &
                (np.isnan(result['delmag']) == False)
            )]
            
            mag_mn = np.median(sID_lc['mag_inst'])
            rms = np.std(sID_lc['mag_inst'])
            mn_err = np.median(sID_lc['magerr_inst'])
            if rms**2 > mn_err**2 :
                fvar = np.sqrt(rms**2 - mn_err**2)
            else : 
                fvar = 0
            N_ep = len(sID_lc)
            
            if N_ep > 10 :
                lc_result.add_row([sID_i, mag_mn, rms, N_ep, mn_err, fvar])

        if ensemble == True : 
        
            avg_rms = np.median(lc_result['RMS'])
            std_rms = np.std(lc_result['RMS'])
            
            lc_result_cut = lc_result[np.where(
                (avg_rms-std_rms<lc_result['RMS']) & (lc_result['RMS']<avg_rms+std_rms)
            )]
            
            sID_good = list(lc_result_cut['sID'])

        else : 
            N_bad = len(lc_result[np.where(lc_result['fvar'] > 0)])

            if N_bad > 20 :
                lc_result_cut = lc_result[np.where(
                    (lc_result['fvar'] == 0 ) &
                    (lc_result['RMS'] < 0.1)
                )]

                sID_good = list(set(lc_result_cut['sID']))
                req_iteration = True
            else :
                sID_good = list(set(lc_result['sID']))
                req_iteration = False


        return sID_good, req_iteration
    
    def reject_star_worker(self, mp_var) :
        result, sID_i = mp_var
        #sID_ls = list(set(result['sID']))

        #sID_i = sID_ls[i]
        if self.filt == 'B' :
            lc = result[np.where(
                (result['sID'] == sID_i)
            )]
        else :
            lc = result[np.where(
                (result['sID'] == sID_i) &
                (result['CLASS_STAR'] > 0.8)
            )]

        if len(lc) > 30 :
            mjd_med = np.median(lc['mjd'])
            mjd_std = np.std(lc['mjd'])
            
            lc_mjd_cut = lc[np.where(
                (mjd_med-3*mjd_std < lc['mjd']) & (lc['mjd'] < mjd_med + 3*mjd_std)
            )]

            if len(lc_mjd_cut) > 20 :
            
                lc_mjd_cut.sort('mjd')
                
                xs = lc_mjd_cut['mjd']
                ys = lc_mjd_cut['mag_inst']
                
                xmin = np.min(lc_mjd_cut['mjd'])
                xmax = np.max(lc_mjd_cut['mjd'])
                
                interp_func = interp1d(xs, ys)
                
                xnew = np.arange(xmin, xmax, 0.1)
                newarr = interp_func(xnew)
                
                result = adfuller(newarr, autolag="AIC")
                pvalue = result[1]
                STD = np.std(ys)
                MED = np.median(ys)

                if self.filt == 'B' :
                    if (STD < MED*0.01) :
                        star_bool = sID_i
                    else :
                        star_bool = False
                else :
                    if (STD < MED*0.003) :
                        star_bool = sID_i
                    else :
                        star_bool = False

            else :
                star_bool = False
        else :
            star_bool = False
        
        return star_bool


    def reject_star(self, result) :

        sID_ls = list(set(result['sID']))
        
        tqdm_bar_fmt = '{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}'
        #CORENUM = 10
        mp_var = [(result, sID_i) for sID_i in sID_ls]

        sID_good = []
        with Pool(CORENUM) as p :
            for output in tqdm(p.imap(self.reject_star_worker, mp_var), total=len(mp_var), bar_format=tqdm_bar_fmt) :
                if type(output) == np.str_ :
                    sID_good.append(output)

        return sID_good, False


    def find_bad_epoch(self, delmag_p2) :
    
        bad_ep_cand = []
        
        sID_ls = list(set(delmag_p2['sID']))
        for i in range(len(sID_ls)) :
            sID_i = sID_ls[i]
            sID_lc = delmag_p2[np.where(
                (delmag_p2['sID'] == sID_i) 
            )]
            
            mn = np.median(sID_lc['mag_inst'])
            std = np.std(sID_lc['mag_inst'])
            
            mn_err = np.median(sID_lc['magerr_inst'])
            std_err = np.std(sID_lc['magerr_inst'])
            
            bad_ep = sID_lc[np.where(
                (mn-2*std > sID_lc['mag_inst']) | (sID_lc['mag_inst'] > mn+2*std) |
                (sID_lc['magerr_inst'] > mn_err + 2*std_err)
            )]
            
            for bad_ep_i in bad_ep :
                bad_ep_cand.append((bad_ep_i['mjd'], bad_ep_i['ampnum']))
        
        bad_ep_ls = []
        
        for i in range(len(list(set(bad_ep_cand)))) : 
            bad_ep_cand_ind = list(set(bad_ep_cand))[i]
            if bad_ep_cand.count(bad_ep_cand_ind) > 1 :
                bad_ep_ls.append(bad_ep_cand_ind)
        
        return bad_ep_ls
    
    def reject_bad_epoch_worker(self, mp_var) :
        #i=0
        #print(i)
        delmag_p2, bad_ep = mp_var
        mjd, amp = bad_ep
        
        bad_ep_idx = list(np.where(
            (delmag_p2['mjd'] == mjd) &
            (delmag_p2['ampnum'] == amp)
        )[0])
        
        #bad_ep_idx_all = bad_ep_idx_all + bad_ep_idx
        return bad_ep_idx


    def reject_bad_epoch(self, delmag_p2) :

        bad_ep_ls = self.find_bad_epoch(delmag_p2)
        print(f'Total of {len(bad_ep_ls)} epochs')
        bad_ep_idx_all = []
        
        '''for i in range(len(bad_ep_ls)) :
            print(i)
            mjd, chip = bad_ep_ls[i]
            
            bad_ep_idx = list(np.where(
                (delmag_p2['mjd'] == mjd) &
                (delmag_p2['chip'] == chip)
            )[0])
            
            bad_ep_idx_all = bad_ep_idx_all + bad_ep_idx'''
        
        tqdm_bar_fmt = '{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}'
        #CORENUM = 10
        mp_var = [(delmag_p2, bad_ep) for bad_ep in bad_ep_ls]

        with Pool(CORENUM) as p :
            for output in tqdm(p.imap(self.reject_bad_epoch_worker, mp_var), total=len(mp_var), bar_format=tqdm_bar_fmt) :
                if type(output) != NoneType :
                    bad_ep_idx_all = bad_ep_idx_all + output
        
        all_idx = np.where(delmag_p2)[0]
        good_ep_idx = all_idx[np.where(np.isin(all_idx, bad_ep_idx_all) == False)]

        return delmag_p2[good_ep_idx]


    def gen_delmag_tab(self,delmag_p3) :

        dtype_result = np.dtype([
            ('mjd', float),
            ('ampnum', int),
            ('delmag', float),
            ('dmag_err', float),
            ('seeing', float)
        ])
        
        result = Table(dtype=dtype_result)
        missing = []
        
        mjd_ls = list(set(delmag_p3['mjd']))
        #chip_ls = ['K', 'M', 'T', 'N']
        amp_ls = list(range(32))
        
        summ = 0
        
        for mjd in mjd_ls :
            for amp in amp_ls : 
                delmag_tab = delmag_p3[np.where(
                    (delmag_p3['mjd'] == mjd) &
                    (delmag_p3['ampnum'] == amp)
                )]
        
                if len(delmag_tab) > 0 :
                    delmag = delmag_tab['delmag'][0]
                    dmag_err = delmag_tab['dmag_err'][0]
                    seeing = delmag_tab['seeing'][0]
                    
                    row = np.array([mjd, amp, delmag, dmag_err, seeing])
                    result.add_row(row)
                else :
                    missing.append((mjd, amp))
        
        return result, missing


    def calculate(self, return_type='delmag') :

        print(f'############# Running DELMAG for {self.tel} {self.field} {self.year} {self.filt} #################')

        print('Calculating initial delmag ...')
        delmag_p1 = self.calc_delmag_all()
        #ascii.write(delmag_p1, 'delmag_p1.cat', overwrite=True)
        

        sID_good, req_itr = self.reject_star(delmag_p1)
        print(f'Number of stars that will be used for calculating delmag: {len(sID_good)}')
        self.IDcut = sID_good
        sID_arr = np.array(sID_good)
        sID_arr = sID_arr[:, np.newaxis]
        sID_good_return = Table(sID_arr)
        #ascii.write(sID_good_return, 'starcand.cat', overwrite=True)

        #itr_num = 1
        #while req_itr == True :
            #itr_num+=1
        #print(f'Calculating delmag for iteration number : {itr_num}')
        #print(type(self.IDcut))
        delmag_p2 = self.calc_delmag_all()
        #sID_good, req_itr = self.reject_star(delmag_p1)
        #print(f'Number of stars that will be used for calculating delmag: {len(sID_good)}')
            #if req_itr == True :
            #    print('Recalculation of delmag required ...')
            #else :
            #    print('Delmag calculation complete.')

            #req_itr = False
        print('Rejecting bad epochs ...')
        delmag_p3 = self.reject_bad_epoch(delmag_p2)
        #ascii.write(delmag_p3, 'delmag_p3.cat', overwrite=True)

        delmag_tab = self.gen_delmag_tab(delmag_p3)

        if return_type == 'delmag' :
            result = delmag_tab
        elif return_type == 'starcat' :
            result = delmag_p3

        print('Finalizing ...')
        return result

    
    def mag_inst_init(self) :
        
        #reference = read_gaia()
        
        #mcat = cdelmag18.gen_mcat(filename, refcat=reference, cat_name='GAIA')
        
        #filename = 'cbftoxkmtc.20180613.034722.new.0002.resamp.cat'
        #filepath = os.path.join(EXTDDIRC, 'N55_2018', 'CTIO', 'B', filename)
        #filepath = '/Users/kamp/Desktop/KAMP/test/data2/N55_2018/CTIO/B/cbftoxkmtc.20180804.048526.new.0003.resamp.cat'
        filepath = os.path.join(self.datadirc, self.filename)
        data = ascii.read(filepath)
        data_cut = data[np.where(
            (data['MAG_AUTO'] < 50) &
            (data['MAGERR_AUTO'] < 50) &
            #(data['MAG_APER'] < 50) &
            #(data['MAGERR_APER'] < 50) &
            (data['FLAGS'] < 5)
        )]
        
        trgt_name_path = os.path.join(DATFDIRC, f'matched_targets_{self.field}.csv')
        trgt_name = ascii.read(trgt_name_path)
        
        data_coord = SkyCoord(
            ra=data_cut['ALPHA_J2000']*u.deg,
            dec=data_cut['DELTA_J2000']*u.deg
        )
        
        trgt_coord = SkyCoord(
            ra=trgt_name['ra']*u.deg,
            dec=trgt_name['dec']*u.deg
        )
        
        idx, d2d, d3d = data_coord.match_to_catalog_sky(trgt_coord)
        
        mcat = data_cut.copy()
        mcat['sep_trgt'] = d2d.to(u.arcsec)
        mcat['Jname'] = trgt_name['Jname'][idx]
        
        mcat_cut = mcat[np.where(mcat['sep_trgt'] < 2.*u.arcsec)]
        
        return mcat_cut
    
    
    def gen_newfilepath(self) : 
        #ls = self.filename.split('.')
        #ls.insert(3,'delmag')
        #filename_new = '.'.join(ls)
        filename_new = self.filename.replace('resamp', 'delmag')
        os.makedirs(os.path.join(DMAGDIRC, f'{self.field}_{self.year}', self.tel, self.filt,), exist_ok=True)
        filepath_new = os.path.join(DMAGDIRC, f'{self.field}_{self.year}', self.tel, self.filt, filename_new)
        return filepath_new

    def mag_inst_engine(self, data, delmag_tab) : 
        #print(self.filename)
        mjd = self.find_mjd()
        
        #delmag_tab = gen_delmag_tab(delmag_p3)
        
        delmag_tab_u = delmag_tab[0]
        delmag_tab_u_chip = delmag_tab_u[np.where(
            (delmag_tab_u['mjd'] == mjd)
        )]
        
        result_ls = []

        if len(delmag_tab_u_chip) > 0 :
        
            for ampnum in range(0, 32) :

                delmag_chip_amp = delmag_tab_u_chip[np.where(delmag_tab_u_chip['ampnum'] == ampnum)]
                data_amp = data[np.where(data['ampnum'] == ampnum)]

                if (len(delmag_chip_amp) > 0) & (len(data_amp) > 0) : 
                    
                    delmag_chip_amp = delmag_chip_amp[0]
                    data_amp['mag_inst'] = data_amp['MAG_AUTO'] + delmag_chip_amp['delmag']
                    data_amp['magerr_inst'] = np.sqrt(data_amp['MAGERR_AUTO']**2 + delmag_chip_amp['dmag_err']**2)
                    data_amp['delmag'] = delmag_chip_amp['delmag']
                    data_amp['dmag_err'] = delmag_chip_amp['dmag_err']
                    data_amp['seeing'] = delmag_chip_amp['seeing']

                    #result = data_amp['']
                    #result_i = data_amp['Jname', 'ALPHA_J2000', 'DELTA_J2000', 'X_IMAGE', 'Y_IMAGE', 'mag_inst', 'magerr_inst', 'FLAGS', 'ELLIPTICITY'].copy()


                    result_i = data_amp.copy()
                    
                    result_ls.append(result_i)

            if len(result_ls) > 0 :
        
                result = vstack(result_ls)
                result['mjd'] = np.array([mjd])
                new_fpath = self.gen_newfilepath()
                ascii.write(result, new_fpath, overwrite=True)

    def apply_delmag_worker(self, mp_var) :
        filename, delmag_tab = mp_var
        self.filename = filename
        #self.chipnum = filename.split('.')[4]
        data = self.mag_inst_init()
        self.mag_inst_engine(data, delmag_tab)

    def apply_delmag(self) :
        delmag_tab = self.calculate()
        filenames = [i for i in os.listdir(self.datadirc) if i.endswith('.cat')]
        mp_var = [(filename, delmag_tab) for filename in filenames]
        
        #self.apply_delmag_worker(mp_var[0])

        tqdm_bar_fmt = '{desc:<5.5}{percentage:3.0f}%|{bar:10}{r_bar}'
        #CORENUM = 10

        with Pool(CORENUM) as p :
            for output in tqdm(p.imap(self.apply_delmag_worker, mp_var), total=len(mp_var), bar_format=tqdm_bar_fmt) :
                if type(output) != NoneType :
                    print('error!!!!')
        #for filename in filenames
        #filename = filenames[0]


def run(INPUTVAR) :
    #data_info = ['CTIO', 'V', 'N55', '2018', '']
    data_info = INPUTVAR + ['']
    cdelmag = DelmagCalculator(data_info)
    cdelmag.apply_delmag()
    #cdelmag.calculate()


'''if __name__ == "__main__" :
    main()'''