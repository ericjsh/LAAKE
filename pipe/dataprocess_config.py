import os
'''
ROOTPATH : path where pipeline is operated.
DATAPATH : path for data (configuration sepcified for data in data_config.json)

WORKDIRC : .ipynb workspace (will be depreciated in the future versions).
PIPEDIRC : path for all pipelines.
DATFDIRC : Extetnal files (apass cat., gaia cat, star cat. etc.) location.

IMGEDIRC : path for running SourceExtractor and storing science image files.
INCTDIRC : path for storing SExtracted catalog files.
AMPMDRIC : path for storing files used for generating ampmaps.
CTLGDIRC : path for storing amp distance measured catalog files that will be used for delmag calculation.
DMAGDIRC : path for storing delmag calculated catalog files.
'''

ROOTPATH = os.getcwd()
DATAPATH = os.path.join(ROOTPATH, 'data')

WORKDIRC = os.path.join(ROOTPATH, 'workspace')
PIPEDIRC = os.path.join(ROOTPATH, 'pipe')
DATFDIRC = os.path.join(ROOTPATH, 'datafile')
OUTPDIRC = os.path.join(ROOTPATH, 'output')

IMGEDIRC = os.path.join(DATAPATH, 'image')
INCTDIRC = os.path.join(DATAPATH, 'initial_cat')
AMPMDIRC = os.path.join(DATAPATH, 'ampmap')
CTLGDIRC = os.path.join(DATAPATH, 'dist_cat')
DMAGDIRC = os.path.join(DATAPATH, 'delmag_cat')

BKUPDIRC = '/Volumes/pKAMP2'
CORENUM = 10

import json

CNFGPATH = os.path.join(PIPEDIRC, 'data_config.json')
with open(CNFGPATH) as f :
    data_config = json.load(f)
