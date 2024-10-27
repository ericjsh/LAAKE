import os
import sys
from shutil import copy, move

from pipe.dataprocess_config import *
from pipe import *


def initialize_laake(INPUTVAR) :

    # Move external headers and science images to local machine
    data_dirs = process_tools.find_directories(INPUTVAR)
    subdir = process_tools.find_subdir(INPUTVAR)

    host_dir = os.path.join(BKUPDIRC, subdir)

    movefile.run(
        execute='copy',
        src=host_dir, 
        dest=data_dirs['image'], 
        extension=data_config['image']['bp_extensions'][0],
        head='c'
    )

    movefile.run(
        execute='copy',
        src=host_dir, 
        dest=data_dirs['ampmap'], 
        extension=data_config['ampmap']['bp_extensions'][0]
    )

    # Move JD.list to local machine
    field, year, *_ = INPUTVAR
    jdlist_fname_old = 'JD.list'
    fy_id = f'{field}_{year}'
    jdlist_fname_new = f'JD_{fy_id}.list'
    jdlist_fpath_new = os.path.join(DATFDIRC, jdlist_fname_new)


    if not os.path.isfile(jdlist_fpath_new) :
        jdlist_fpath_old = os.path.join(BKUPDIRC, fy_id, jdlist_fname_old)
        jdlist_fpath_old2new = os.path.join(DATFDIRC, jdlist_fname_old)
        copy(jdlist_fpath_old, jdlist_fpath_old2new)
        move(jdlist_fpath_old2new, jdlist_fpath_new)



def execute_laake(INPUTVAR) :

    photometry.run(INPUTVAR)

    file_checker.run(INPUTVAR)

    starcand.starcat_generator(INPUTVAR)

    astrom_check.run(INPUTVAR)

    ampmap.run(INPUTVAR)

    dist_measure.run(INPUTVAR)

    delmag_calculator.run(INPUTVAR)

    synthesize.run(INPUTVAR)


def finalize_laake(INPUTVAR) :
    data_dirs_local = process_tools.find_directories(INPUTVAR)
    data_dirs_backup = process_tools.find_directories_backup(INPUTVAR)

    backup_keys = ['image', 'ampmap', 'initial', 'delmag', 'dist']

    for backup_key in backup_keys :

        print(f'Backing up for {backup_key}...')
        backup_dir = data_dirs_backup[backup_key]
        print(backup_dir)

        local_dir = data_dirs_local[backup_key]

        movefile.run_all(
            execute='move',
            src=local_dir, 
            dest=backup_dir
        )

def main(INPUTVAR) :
    field, year, tel, filt = INPUTVAR

    print(
        f'''
        =============================================================
             __           ____       ____     __   ___   ________
            |  |         /    \     /    \   |  | /  /  |   _____|
            |  |        /  /\  \   /  /\  \  |  |/  /   |  |_____
            |  |       |  /__\  | |  /__\  | |      \   |   _____|
            |  |_____  |   __   | |   __   | |  |\   \  |  |_____
            |________| |__|  |__| |__|  |__| |__| \___\ |________|
            Lightcurve  Assembly Architecture for KAMP   Extracts
        =============================================================
        LAAKE ver 1.0.0
        Running for:
            - field: {field}
            - year: {year}
            - telescope: {tel}
            - filter: {filt}
        '''
        )
    
    initialize_laake(INPUTVAR)
    execute_laake(INPUTVAR)
    finalize_laake(INPUTVAR)

if __name__ == "__main__" :
    INPUTVAR = sys.argv[1:]
    #INPUTVAR = ['N55', '2017', 'CTIO', 'V']
    main(INPUTVAR)



