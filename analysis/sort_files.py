import glob
from astropy.io import fits
import numpy as np
import os
import shutil
from astropy.time import Time

import warnings
warnings.filterwarnings('ignore')

def sort_files(indir, channel='ch1', objname=None, one_epoch=False):

    dr = os.path.join(indir, 'rawdata')
    objs = []
    mjds = []
    files = []

    groups = []

    fls = glob.glob(os.path.join(dr, '*_cbcd.fits'))
    fls.sort()

    i=0
    for fl in fls:
        f = fits.open(fl, mode='update')
        h = f[0].header
        
        obj = h['OBJECT']
        if objname:
            obj = objname

        mjd = h['MJD_OBS']

        f[0].header['CHNLNUM']='Ch'+str(h['CHNLNUM']).strip()
        f[0].header['MJD']=h['MJD_OBS']

        if h['BUNIT']=='MJy/sr':
            convert_factor = 1.0e6 / (180 / np.pi * 3600)**2
            pixscale = np.sqrt(h['CD1_1']**2+h['CD1_2']**2) * 3600
            zpt = -2.5 * np.log10(convert_factor * pixscale**2 / 3631)

            f[0].header['ZPTPSCAL']=pixscale
            f[0].header['ZPTMAG']=zpt

        mjds = [g['mjds'] for g in groups]
        new_group=True
        for j,m in enumerate(mjds):
            if ((np.min(m) > mjd - 20. and np.min(m) < mjd + 20.) or one_epoch):
                groups[j]['idx'].append(i)
                groups[j]['mjds'].append(mjd)
                groups[j]['files'].append(fl)
                new_group=False
                break

        if new_group:
            group = {}
            group['idx']=[i]
            group['mjds']=[mjd]
            group['files']=[fl]
            groups.append(group)

    print('There are {0} groups'.format(len(groups)))

    # Sort in order of ascending MJD
    groups = sorted(groups, key=lambda x: np.min(x['mjds']))

    for i,group in enumerate(groups):
        t0 = Time(np.min(group['mjds']), format='mjd')
        newbasedir = t0.datetime.strftime('ut%y%m%d')
        newbasedir = os.path.join(indir, newbasedir, channel)
        if not os.path.exists(newbasedir):
            os.makedirs(newbasedir)

        print(f'Copying files for group {i}')

        for file in group['files']:
            maskfile = file.replace('cbcd.fits','bimsk.fits')
            uncrfile = file.replace('cbcd.fits','cbunc.fits')

            shutil.copyfile(file, os.path.join(newbasedir,
                os.path.basename(file)))
            shutil.copyfile(maskfile, os.path.join(newbasedir,
                os.path.basename(maskfile)))
            shutil.copyfile(uncrfile, os.path.join(newbasedir,
                os.path.basename(uncrfile)))

