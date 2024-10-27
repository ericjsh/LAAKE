import os
from shutil import copy

from pipe.dataprocess_config import *

from pipe import process_tools
from pipe import movefile



class PhotometryExecutor :

    def __init__(self, INPUTVAR: list) -> None :
        self.data_dirs = process_tools.find_directories(INPUTVAR)
    
    def run_source_extractor(self) -> None :

        sex_files_path = os.path.join(DATFDIRC, 'sex_files')
        sex_files_name = os.listdir(sex_files_path)

        movefile.locate_files(
            execute='copy',
            filenames=sex_files_name,
            src=sex_files_path,
            dest=self.data_dirs['image']
        )

        os.chdir(self.data_dirs['image'])
        os.system('python kmtnet_sex.py')
        os.chdir(ROOTPATH)


    def path_finalize(self) -> None :

        movefile.run(
            execute='move',
            src=self.data_dirs['image'],
            dest=self.data_dirs['initial'],
            extension='.cat'
        )


def run(INPUTVAR: list) -> None :
    #INPUTVAR = ['N55', '2018', 'CTIO', 'V']
    PE = PhotometryExecutor(INPUTVAR)
    PE.run_source_extractor()
    PE.path_finalize()

