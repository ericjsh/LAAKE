import os
import sys
import argparse
from shutil import copy, move
import time

import pipe.agn_finder
from pipe.dataprocess_config import *
from pipe import *
import pipe

def laake_logger(INPUTVAR, process) :
    field, year, tel, band = INPUTVAR

    logger_fname = "laake_logger.txt"

    if not os.path.isfile(logger_fname) : 
        with open(logger_fname, "w") as f :
            f.write("Timestamp            Field     Year      Band      Telescope LAAKE ver Process         ")
            f.write('\n')

    with open("laake_logger.txt", "a") as f :
        timestamp = time.strftime("%Y-%m-%dT%H:%M:%S")
        f.write(f'{timestamp}  {field:<9} {year}      {band}         {tel:<9} {pipe.__version__:<9} {process}')
        f.write('\n')


def initialize_laake(INPUTVAR) :

    laake_logger(INPUTVAR, 'Initializing LAAKE')

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


def execute_laake(INPUTVAR) :

    laake_logger(INPUTVAR, 'laake.photometry')
    pipe.photometry.run(INPUTVAR)

    laake_logger(INPUTVAR, 'laake.file_checker')
    pipe.file_checker.run(INPUTVAR)

    laake_logger(INPUTVAR, 'laake.starcand.starcat_generator')
    pipe.starcand.starcat_generator(INPUTVAR)

    laake_logger(INPUTVAR, 'laake.astrom_check')
    pipe.astrom_check.run(INPUTVAR)

    laake_logger(INPUTVAR, 'laake.ampmap')
    pipe.ampmap.run(INPUTVAR)

    laake_logger(INPUTVAR, 'laake.dist_measure')
    pipe.dist_measure.run(INPUTVAR)

    laake_logger(INPUTVAR, 'laake.delmag_calculator')
    pipe.delmag_calculator.run(INPUTVAR)

    laake_logger(INPUTVAR, 'laake.synthesize')
    pipe.synthesize.run(INPUTVAR)


def finalize_laake(INPUTVAR) :

    laake_logger(INPUTVAR, 'Finalizing LAAKE')

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
    laake_logger(INPUTVAR, 'Complete')

if __name__ == "__main__" :

    parser = argparse.ArgumentParser()

    command_group = parser.add_mutually_exclusive_group()
    command_group.add_argument('-v', '--version', action='store_true', help='check LAAKE version')
    command_group.add_argument('-e', '--emblem', action='store_true', help='print LAAKE banner')
    command_group.add_argument('-g', '--gui', action='store_true', help='open light curve GUI')
    command_group.add_argument('-u', '--update', action='store_true', help='update light curves')
    command_group.add_argument('-q', '--queue', action='store_true', help='run data processing pipeline from queue')
    command_group.add_argument('-r', '--run', action='store_true', help='run data processing pipeline')

    required = parser.add_argument_group('required arguments for -r, --run')
    required.add_argument('-f', '--field', action='append', type=str, help='input field name')
    required.add_argument('-y', '--year', action='append', type=str, help='input year')
    required.add_argument('-b', '--band', action='append', type=str, help='input band (filter) name')
    required.add_argument('-t', '--tel', action='append', type=str, help='input telescope name')

    args = parser.parse_args()

    if args.emblem :
        print(laake_emblem())
        sys.exit(0)

    if args.gui :
        print("Developer's Note: GUI Will be updated in later version")
        sys.exit(0)

    if args.version :
        print(pipe.__version__)
        sys.exit(0)

    if args.update :
        if (args.tel==None) | (args.field==None) :
            raise ValueError("Some parameters are missing")

        if args.tel[0] not in ['CTIO', 'SSO', 'SAAO'] :
            raise ValueError(f"{args.tel[0]} is not a valid telescope name")
        
        inputvar = [args.field[0], '2016', args.tel[0], 'B']
        print(f'Generating AGN light curves for {args.field[0]} {args.tel[0]}')
        pipe.agn_finder.run(inputvar)
        print(f'Complete')
        sys.exit(0)

    if args.queue :

        queue_fname = "laake_pipe_queue.txt"

        if not os.path.isfile(queue_fname) :
            with open(queue_fname, 'w') as f :
                f.write("Field Year Telescope Band")
                f.write("\n")

        while True :

            with open(queue_fname, "r+") as f :
                lines = f.readlines()
                if len(lines) > 2 :
                    queue = lines[1]
                    f.seek(0)
                    for line in lines :
                        if line != queue :
                            f.write(line)
                    f.truncate()
                else :
                    print('end')
                    break

            inputvar = queue.split('\n')[:-1][0].split(' ')
            laake_pipeline(inputvar)

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
        #INPUTVAR = ['N55', '2017', 'CTIO', 'V']
        laake_pipeline(inputvar)

    else :
        print(laake_emblem() + '''
        Written by Sungho JUNG <eric2912@snu.ac.kr>
            
        visit https://congruous-value-af8.notion.site/LAAKE-12c48de1964b80e09902ce791bc1d833?pvs=4
            
        > SYNTAX python laake.py [-h] [-v | -e | -g | -r] [-f FIELD] [-y YEAR] [-b BAND] [-t TEL]
        > to show help: python laake.py -h
        > example command: python laake.py -r 
        ''')



