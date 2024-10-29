import os
import sys
import argparse
from shutil import copy, move

from pipe.dataprocess_config import *
from pipe import *
import pipe

def initialize_laake(INPUTVAR) :

    # Move external headers and science images to local machine
    data_dirs = pipe.process_tools.find_directories(INPUTVAR)
    subdir = pipe.process_tools.find_subdir(INPUTVAR)

    host_dir = os.path.join(BKUPDIRC, subdir)

    pipe.movefile.run(
        execute='copy',
        src=host_dir, 
        dest=data_dirs['image'], 
        extension=data_config['image']['bp_extensions'][0],
        head='c'
    )

    pipe.movefile.run(
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

    pipe.photometry.run(INPUTVAR)

    pipe.file_checker.run(INPUTVAR)

    pipe.starcand.starcat_generator(INPUTVAR)

    pipe.astrom_check.run(INPUTVAR)

    pipe.ampmap.run(INPUTVAR)

    pipe.dist_measure.run(INPUTVAR)

    pipe.delmag_calculator.run(INPUTVAR)

    pipe.synthesize.run(INPUTVAR)


def finalize_laake(INPUTVAR) :
    data_dirs_local = pipe.process_tools.find_directories(INPUTVAR)
    data_dirs_backup = pipe.process_tools.find_directories_backup(INPUTVAR)

    backup_keys = ['image', 'ampmap', 'initial', 'delmag', 'dist']

    for backup_key in backup_keys :

        print(f'Backing up for {backup_key}...')
        backup_dir = data_dirs_backup[backup_key]
        print(backup_dir)

        local_dir = data_dirs_local[backup_key]

        pipe.movefile.run_all(
            execute='move',
            src=local_dir, 
            dest=backup_dir
        )


def laake_emblem() :
    banner = f'''
        =============================================================
             __           ____       ____     __   ___   ________
            |  |         /    \     /    \   |  | /  /  |   _____|
            |  |        /  /\  \   /  /\  \  |  |/  /   |  |_____
            |  |       |  /__\  | |  /__\  | |      \   |   _____|
            |  |_____  |   __   | |   __   | |  |\   \  |  |_____
            |________| |__|  |__| |__|  |__| |__| \___\ |________|
            Lightcurve  Assembly Architecture for KAMP   Extracts
        =============================================================
        LAAKE ver {pipe.__version__}
    '''
    return banner

def laake_pipeline(INPUTVAR) :
    field, year, tel, filt = INPUTVAR

    print(laake_emblem() + 
        f'''
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

    parser = argparse.ArgumentParser()

    command_group = parser.add_mutually_exclusive_group()
    command_group.add_argument('-v', '--version', action='store_true', help='check LAAKE version')
    command_group.add_argument('-e', '--emblem', action='store_true', help='print LAAKE banner')
    command_group.add_argument('-g', '--gui', action='store_true', help='open light curve GUI')
    command_group.add_argument('-r', '--run', action='store_true', help='run data processing pipeline')

    required = parser.add_argument_group('required arguments for -r, --run')
    required.add_argument('-f', '--field', action='append', type=str, help='input field name')
    required.add_argument('-y', '--year', action='append', type=str, help='input year')
    required.add_argument('-b', '--band', action='append', type=str, help='input band (filter) name')
    required.add_argument('-t', '--tel', action='append', type=str, help='input telescope name')

    args = parser.parse_args()
    #print(args.__dict__)
    #print(len(args.__dict__))

    if args.emblem :
        print(laake_emblem())
        sys.exit(0)

    if args.gui :
        print("Developer's Note: GUI Will be updated in later version")
        sys.exit(0)

    if args.version :
        print(pipe.__version__)
        sys.exit(0)

    if args.run :

        if (args.year==None) | (args.band==None) | (args.tel==None) | (args.field==None) :
            raise ValueError("Some parameters are missing")

        if not args.year[0].isnumeric() : 
            raise ValueError(f"{args.year[0]} is not a valid year")

        if args.band[0] not in ['B', 'V', 'I'] :
            raise ValueError(f"{args.band[0]} is not a valid band (filter) name")

        if args.tel[0] not in ['CTIO', 'SSO', 'SAAO'] :
            raise ValueError(f"{args.tel[0]} is not a valid telescope name")

        inputvar = [args.field[0], args.year[0], args.tel[0], args.band[0]]
        #INPUTVAR = sys.argv[1:]
        #INPUTVAR = ['N55', '2017', 'CTIO', 'V']
        laake_pipeline(inputvar)

    else :
        #parser.print_help()
        print(laake_emblem() + '''
        Written by Sungho JUNG <eric2912@snu.ac.kr>
            
        visit https://congruous-value-af8.notion.site/LAAKE-12c48de1964b80e09902ce791bc1d833?pvs=4
            
        > SYNTAX python laake.py [-h] [-v | -e | -g | -r] [-f FIELD] [-y YEAR] [-b BAND] [-t TEL]
        > to show help: python laake.py -h
        > example command: python laake.py -r 
        ''')



