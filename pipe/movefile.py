from pipe.dataprocess_config import *

import os
from shutil import move, copy
from concurrent.futures import ThreadPoolExecutor

from pipe import process_tools

#src_inp = '/Users/eric/Desktop/Home/College/Intern/KAMP/test/data2/./N55_2018/CTIO/V'
#dest_inp = '/Volumes/KAMP/data1/./N55_2018/CTIO/V'


def locate_files(execute , filenames, src, dest) -> None:

    func_dict = {
        'move' : move,
        'copy' : copy
    }

    for filename in filenames:
 
        src_path = os.path.join(src, filename)
        dest_path = os.path.join(dest, filename)
        func_dict[execute](src_path, dest_path)
        #print(f'.moved {src_path} to {dest_path}', flush=True)


def run(execute, src: str, dest: str, extension: str, head=None) -> None:
    os.makedirs(dest, exist_ok=True)
    if type(head) == str :
        files = [name for name in os.listdir(src) if name.startswith(head) & name.endswith(extension)]
        #files.sort()
        #files = files[:100]
    else :
        files = [name for name in os.listdir(src) if name.endswith(extension)]
        #files.sort()
        #files = files[:100]
    print(len(files))
    n_workers = 100 + 5*(len(files) < 10)
    chunksize = round(len(files) / n_workers)
    if chunksize == 0 :
        chunksize += 1
    with ThreadPoolExecutor(n_workers) as exe:
        for i in range(0, len(files), chunksize):
            filenames = files[i:(i + chunksize)]
            _ = exe.submit(locate_files, execute, filenames, src, dest)
    print('Done')


def run_all(execute, src: str, dest: str) -> None:
    os.makedirs(dest, exist_ok=True)
    files = [name for name in os.listdir(src)]

    if len(files) > 1:
        n_workers = 100 + 5*(len(files) < 100)
        chunksize = round(len(files) / n_workers)
        with ThreadPoolExecutor(n_workers) as exe:
            for i in range(0, len(files), chunksize):
                filenames = files[i:(i + chunksize)]
                _ = exe.submit(locate_files, execute, filenames, src, dest)
        print('Done')


'''INPUTVAR = ['N55', '2017', 'CTIO', 'V']

data_dirs = process_tools.find_directories(INPUTVAR)
subdir = process_tools.find_subdir(INPUTVAR)

backup_rootpath = '/Volumes/DHSON_HDD_Backup/pkamp-dataprocess_backup/data'


backup_keys = ['image', 'ampmap', 'initial', 'delmag', 'dist']

for backup_key in backup_keys :

    print(f'Backing up for {backup_key}...')
    backup_dir = os.path.join(backup_rootpath, backup_key, subdir)

    local_dir = data_dirs[backup_key]

    run_all(
        execute='move',
        src=local_dir, 
        dest=backup_dir
    )

data_dirs'''