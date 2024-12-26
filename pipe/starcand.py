import os

import numpy as np

from astropy.io import ascii

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import match_coordinates_sky, SkyCoord
from astroquery.vizier import Vizier

from astropy.table.table import Table

from pipe.dataprocess_config import *


class ReadFiles :

    '''
    Reading KMTN, GAIA, and APASS catalogues.
    DATFDIRC is where all the catalogues are saved at.
    '''

    def __init__(self, INPUTVAR: list) :
        self.field, self.year, self.tel, self.filt = INPUTVAR
        self.kmtn_fname = f'matched_targets_{self.field}.csv'

    def _kmtn_cat(self) -> Table :
        '''
        Reads KMTNet source catalog file as astropy.table.table.Table format.
        KMTNet source catalog should be named as 'matched_targets_{field}.csv'
        '''
        #kmtn_cat_path = os.path.join(DATFDIRC, 'matched_targets_N55.csv')
        kmtn_cat_path = os.path.join(DATFDIRC, self.kmtn_fname)
        kmtn_cat = ascii.read(kmtn_cat_path)
        self.cent_ra = np.median(kmtn_cat['ra'])
        self.cent_dec = np.median(kmtn_cat['dec'])
    
        return kmtn_cat
    
    
    def _gaia_cat(self) -> Table :
        '''
        Calling GAIA DR3 via Vizier as astropy.table.table.Table format.
        Catalogue cut at gaia_cat_cut illustrated at:
        https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
        Catalouge should be named as 'gaia_{field}.csv'
        '''
        gaia_cat_fname = f'gaia_{self.field}.csv'
        
        if gaia_cat_fname not in os.listdir(DATFDIRC) :

            print(f'{gaia_cat_fname} not in local directory. Downloading from Vizier ...')
            Vizier.ROW_LIMIT = -1
            
            gaia = Vizier.query_region(
                coord.SkyCoord(
                    ra=self.cent_ra, dec=self.cent_dec, 
                    unit=(u.deg, u.deg), frame='icrs'
                ),
                radius=2*u.deg,
                catalog=['Gaia'],
            )
            
            gaia_cat = gaia['I/355/gaiadr3']

            output_fpath = os.path.join(DATFDIRC, gaia_cat_fname)
            ascii.write(gaia_cat, output_fpath, overwrite=True)

        gaia_cat_path = os.path.join(DATFDIRC, gaia_cat_fname)
        gaia_cat = ascii.read(gaia_cat_path)
    
        gaia_cat_cut = gaia_cat[np.where(
            (-0.5 < gaia_cat['BP-RP']) & (gaia_cat['BP-RP'] < 2.75)
        )]
    
        return gaia_cat_cut
    
    
    def _apass_cat(self) -> Table :
        '''
        Reads KMTNet source catalog file as astropy.table.table.Table format.
        KMTNet source catalog should be named as 'apass_dr_10_{field}.csv'
        '''
        aps_cat_fname = f'apass_dr10_{self.field}.csv'

        if aps_cat_fname not in os.listdir(DATFDIRC) : 

            print(f'{aps_cat_fname} not in local directory. Downloading from Vizier ...')
            Vizier.ROW_LIMIT = -1
            
            aps = Vizier.query_region(
                coord.SkyCoord(
                    ra=self.cent_ra, dec=self.cent_dec, 
                    unit=(u.deg, u.deg), frame='icrs'
                ),
                radius=2*u.deg,
                catalog=['APASS'],
            )
            
            aps_cat = aps[0]

            output_fpath = os.path.join(DATFDIRC, aps_cat_fname)
            ascii.write(aps_cat, output_fpath, overwrite=True)

        aps_cat_path = os.path.join(DATFDIRC, f'apass_dr10_{self.field}.csv')
        aps_cat = ascii.read(aps_cat_path)
    
        aps_cat_cln = aps_cat[(~aps_cat['Bmag'].mask) & (~aps_cat['Vmag'].mask)]
    
        return aps_cat_cln
    

def _match_kmtn_gaia(kmtn_cat: Table, gaia_cat: Table) -> Table :
    '''
    Cross matching sources from preKAMP NGC55 and GAIA DR3
    '''
    kmtn_radec = SkyCoord(
        ra=kmtn_cat['ra']*u.deg,
        dec=kmtn_cat['dec']*u.deg
    )
    
    gaia_cat_radec = SkyCoord(
        ra=gaia_cat['RAJ2000']*u.deg,
        dec=gaia_cat['DEJ2000']*u.deg
    )
    
    mcat = kmtn_cat.copy()
    
    idx, d2d, d3d = kmtn_radec.match_to_catalog_sky(gaia_cat_radec)
    
    sigmaG_0 = 0.0027553202
    sigmaGBP_0 = 0.0027901700
    sigmaGRP_0 = 0.0037793818
    
    mcat['sep_KG'] = d2d.to(u.arcsec)
    mcat['Plx'] = gaia_cat['Plx'][idx]
    mcat['BP'] = gaia_cat['BPmag'][idx]
    mcat['e_BP'] = np.sqrt((-2.5/np.log(10)*gaia_cat['e_FBP'][idx]/gaia_cat['FBP'][idx])**2 + sigmaGBP_0**2)
    mcat['RP'] = gaia_cat['RPmag'][idx]
    mcat['e_RP'] = np.sqrt((-2.5/np.log(10)*gaia_cat['e_FRP'][idx]/gaia_cat['FRP'][idx])**2 + sigmaGRP_0**2)
    mcat['Gmag'] = gaia_cat['Gmag'][idx]
    mcat['e_Gmag'] = np.sqrt((-2.5/np.log(10)*gaia_cat['e_FG'][idx]/gaia_cat['FG'][idx])**2 + sigmaG_0**2)
    
    mcat_cut = mcat[np.where(
        (mcat['sep_KG'] < 2.*u.arcsec)
    )]

    return mcat_cut


def match_all(INPUTVAR: list) -> Table :
    '''
    Cross matches KMTNet source catalog, GAIA DR3, and APASS DR10
    '''
    catalogs = ReadFiles(INPUTVAR)
    kmtn_cat = catalogs._kmtn_cat()
    gaia_cat = catalogs._gaia_cat()
    aps_cat = catalogs._apass_cat()

    mcat_kmtn_gaia = _match_kmtn_gaia(kmtn_cat, gaia_cat)

    mcat_kmtn_gaia_radec = SkyCoord(
        ra=mcat_kmtn_gaia['ra']*u.deg,
        dec=mcat_kmtn_gaia['dec']*u.deg
    )
    
    aps_cat_radec = SkyCoord(
        ra=aps_cat['RAJ2000']*u.deg,
        dec=aps_cat['DEJ2000']*u.deg
    )
    
    mcat = mcat_kmtn_gaia.copy()
    
    idx, d2d, d3d = mcat_kmtn_gaia_radec.match_to_catalog_sky(aps_cat_radec)
    
    mcat['sep_MA'] = d2d.to(u.arcsec)
    mcat['Bref'] = aps_cat['Bmag'][idx]
    mcat['e_Bref'] = aps_cat['e_Bmag'][idx]
    mcat['Vref'] = aps_cat['Vmag'][idx]
    mcat['e_Vref'] = aps_cat['e_Vmag'][idx]
    
    mcat_cut = mcat[np.where(mcat['sep_MA'] < 2.*u.arcsec)]

    return mcat_cut


def _gaia2JC_function(x: float, filter: str) -> float :
    '''
    GAIA polynomial, f(x) = a0 + a1*x + a2*x**2, as defined in 
    https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html

    G - V = f(GBP - GRP) --> V = G - f(GBP - GRP)
    G - I = f(GBP - GRP) --> I = G - f(GBP - GRP)
    
    '''
    if filter == 'V' :
        a0 = -0.01760
        a1 = -0.006860
        a2 = -0.1732
    elif filter == 'I' :
        a0 = 0.02085
        a1 = 0.7419
        a2 = -0.09631
        
    return a0 + a1*x + a2*x*x


def gaia2JC(mcat_all: Table) -> Table :
    '''
    Run _gaia2JC_function for all rows
    '''

    mcat_all['V_gaia'] = mcat_all['Gmag'] - _gaia2JC_function(mcat_all['BP'] - mcat_all['RP'], 'V')
    mcat_all['e_V_gaia'] = np.sqrt(mcat_all['e_Gmag']**2 + (mcat_all['BP']/mcat_all['V_gaia']*mcat_all['e_BP'])**2 + (mcat_all['RP']/mcat_all['V_gaia']*mcat_all['e_RP'])**2)
    
    mcat_all['I_gaia'] = mcat_all['Gmag'] - _gaia2JC_function(mcat_all['BP'] - mcat_all['RP'], 'I')
    mcat_all['e_I_gaia'] = np.sqrt(mcat_all['e_Gmag']**2 + (mcat_all['BP']/mcat_all['I_gaia']*mcat_all['e_BP'])**2 + (mcat_all['RP']/mcat_all['I_gaia']*mcat_all['e_RP'])**2)

    return mcat_all


def parallax_cut(mcat_all: Table, filt: str) -> Table :
    '''
    Select star candidates by removing sources that exist outside of our Milky Way
    (i.e. remove sources with parallex > 0.038)
    '''
    mcat_all['sID'] = np.array([f's{i}' for i,_ in enumerate(mcat_all)])

    result = mcat_all[np.where(
        (mcat_all['e_Gmag'] < 0.01) &
        (mcat_all['e_BP'] < 0.01) &
        (mcat_all['e_RP'] < 0.01)
    )]

    if filt == 'B' : 
    
        result['Vref'] = result['Vref'].astype(float)
        result['e_Vref'] = result['e_Vref'].astype(float)
        result['Bref'] = result['Bref'].astype(float)
        result['e_Bref'] = result['e_Bref'].astype(float)
        
        result_sub = result[np.where(
            np.abs(result['Vref'].astype(float) - result['V_gaia']) < 0.086
        )]
    else :
        result_sub = result.copy()
    
    result_final = result_sub[np.where(
        (result_sub['Plx'] > 0.05)
    )]


    return result_final


def starcat_generator(INPUTVAR: list) -> None :
    '''
    Generates star candidate catalog named as 'starcand_{field}.csv'.
    INPUTVAR should be in the format of [{field}, {KMTNet source catalog}]
    '''
    #INPUTVAR = ['N55', '2018', 'CTIO', 'I']
    field, _, _, filt = INPUTVAR

    print("##### Running STARCAND #####")
    print("Loading KMTN, GAIA, APASS catalogs ...")
    if filt == 'B' :
        mcat_all = match_all(INPUTVAR)

    else : 
        catalogs = ReadFiles(INPUTVAR)
        kmtn_cat = catalogs._kmtn_cat()
        gaia_cat = catalogs._gaia_cat()

        mcat_all = _match_kmtn_gaia(kmtn_cat, gaia_cat)

    print("Converting GAIA filter system to Johnson-Cousin filter system ...")
    mcat = gaia2JC(mcat_all)
    
    print("Finializing ...")
    starcand_cat = parallax_cut(mcat, filt)


    output_fname = f'starcand_{field}.csv'
    output_fpath = os.path.join(DATFDIRC, output_fname)
    ascii.write(starcand_cat, output_fpath, overwrite=True)

    print(f"Done! Total of {len(starcand_cat)} star candidates saved at {output_fname}")


'''if __name__ == "__main__" :
    INPUTVAR = sys.argv[1:]
    main(INPUTVAR)'''


