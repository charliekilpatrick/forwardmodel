import glob
from astropy.io import fits
from numpy import *
import os
import shutil
from astropy.time import Time

import warnings
warnings.filterwarnings('ignore')

def sort_files(indir, channel='ch1', objname=None, one_epoch=False):

    drs = glob.glob(os.path.join(indir, 'r????????'))
    drs.sort()
    objs = []
    mjds = []

    for dr in drs:
        print(f'Checking {dr}')

        fls = glob.glob(os.path.join(dr, f'{channel}/bcd/*_cbcd.fits'))
        fls.sort()

        if len(fls)==0: continue

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

    groups = []
    for i in range(len(drs)):
        if any([i in g['idx'] for g in groups]):
            continue

        group = {}
        group['idx']=[i]
        group['drs']=[drs[i]]
        group['mjds']=[mjds[i]]
        expr = os.path.join(drs[i], channel, 'bcd', '*cbcd.fits')
        group['files']=list(glob.glob(expr))
        group['epoch']=None

        for j in range(len(drs)):
            if j != i:
                if ((mjds[j] > mjds[i] - 20. and mjds[j] < mjds[i] + 20.) or
                    one_epoch):
                    group['idx'].append(j)
                    group['drs'].append(drs[j])
                    group['mjds'].append(mjds[j])
                    expr = os.path.join(drs[j], channel, 'bcd', '*cbcd.fits')
                    group['files'].extend(list(glob.glob(expr)))

        # Sort files by name
        group['files']=sorted(group['files'])

        groups.append(group)

    print('There are {0} groups'.format(len(groups)))

    # Sort in order of ascending MJD
    groups = sorted(groups, key=lambda x: min(x['mjds']))

    for i,group in enumerate(groups):
        t0 = Time(min(group['mjds']), format='mjd')
        newbasedir = t0.datetime.strftime('ut%y%m%d')
        newbasedir = os.path.join(indir, newbasedir, channel)
        if not os.path.exists(newbasedir):
            os.makedirs(newbasedir)

        print(f'Copying files for group {i}')

        for file in group['files']:
            maskfile = file.replace('cbcd.fits','bimsk.fits')
            uncrfile = file.replace('cbcd.fits','cbunc.fits')

            # Update original file with specific header values
            hdu = fits.open(file, mode='update')
            h = hdu[0].header
            hdu[0].header['CHNLNUM']='Ch'+str(h['CHNLNUM']).strip()
            hdu[0].header['MJD']=h['MJD_OBS']
            if h['BUNIT']=='MJy/sr':
                convert_factor = 1.0e6 / (180 / pi * 3600)**2
                pixscale = sqrt(h['CD1_1']**2+h['CD1_2']**2) * 3600
                zpt = -2.5 * log10(convert_factor * pixscale**2 / 3631)

                hdu[0].header['ZPTPSCAL']=pixscale
                hdu[0].header['ZPTMAG']=zpt

            hdu.close()

            shutil.copyfile(file, os.path.join(newbasedir,
                os.path.basename(file)))
            shutil.copyfile(maskfile, os.path.join(newbasedir,
                os.path.basename(maskfile)))
            shutil.copyfile(uncrfile, os.path.join(newbasedir,
                os.path.basename(uncrfile)))

    # Clean up original directories
    for dr in drs:
        if os.path.exists(dr) and os.path.isdir(dr):
            print(f'Deleting {dr}')
            shutil.rmtree(dr)

