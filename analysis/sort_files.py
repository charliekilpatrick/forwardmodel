import glob
from astropy.io import fits
from numpy import *
import os
import shutil

def sort_files(indir, channel='ch1', objname=None):

    drs = glob.glob(os.path.join(indir, 'r????????'))
    drs.sort()
    objs = []
    mjds = []

    for dr in drs:
        print(f'Checking {dr}')

        fls = glob.glob(os.path.join(dr, f'{channel}/bcd/*_cbcd.fits'))
        fls.sort()

        print(f'First file {fls[0]}')

        f = fits.open(fls[0], mode='readonly')
        obj = f[0].header['OBJECT']
        if objname:
            obj = objname

        mjd = f[0].header['MJD_OBS']

        mjds.append(mjd)
        objs.append(obj)

        for subdir in glob.glob(os.path.join(dr, '*')):
            if channel not in subdir:
                print(f'Deleting {subdir}')
                if os.path.isdir(subdir):
                    # Remove directory
                    shutil.rmtree(subdir)
                else:
                    # Remove file
                    os.remove(subdir)


    drs = array(drs)
    objs = array(objs)
    mjds = array(mjds)

    inds = argsort(mjds)

    drs = drs[inds]
    objs = objs[inds]
    mjds = mjds[inds]

    assert all(mjds[1:] - mjds[:-1] > 0)

    epochs = {}
    for obj in unique(objs):
        epochs[obj] = 1

    for i in range(len(drs)):
        outname = f'{objs[i]}:{str(epochs[objs[i]]).zfill(3)}'
        outname = os.path.basename(outname)
        outname = os.path.join(indir, outname)
        print(f'Moving {drs[i]}->{outname}')
        shutil.move(drs[i], outname)
        epochs[objs[i]] += 1
