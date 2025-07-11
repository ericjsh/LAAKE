# dev. ver. 250227


# STEP1: change the variables
import os
from pipe import movefile, kmtn_refcat

extdir = '/Volumes/preKAMP2017'
fieldname = 'N3585'
year = '2017'

vanilla_dir = '/Users/kamp/Desktop/KAMP/LAAKE/vanilla_workspace'
tel = 'CTIO'
band = 'V'
INPUTVAR = [fieldname, year, tel, band]
extpath = os.path.join(extdir, f'{fieldname}_{year}', tel, band)


# STEP2: check if correct ones came in
filenames = [f for f in os.listdir(extpath) if f.endswith('.resamp.cat')]
print(filenames[:10])
print(len(filenames))


# STEP3: Move files and go see if files correctly came in
movefile.run(execute='copy', src=extpath, dest=vanilla_dir, extension='.resamp.cat')

# STEP4: Generate reference epoch
kmtn_refcat.kmtn_refcat_gen(INPUTVAR, vanilla=True)

# STEP5: don't forget to erase vanilla!!!!