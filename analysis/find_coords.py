from numpy import *
import glob
from astropy import wcs
from astropy.io import fits
import os

def find_coords(basedir):
    mosaic = os.path.join(basedir, '*:001/working_dir/Combine/mosaic.fits')
    fls = glob.glob(mosaic)

    print('There are {0} total files to analyze in {1}'.format(len(fls),
        basedir))

    outdata = {}

    for fl in fls:
        print(f'Analyzing {fl}')

        f = fits.open(fl)
        w = wcs.WCS(f[0].header)
        dat = f[0].data

        newdat = dat*0.

        for i in range(int(len(dat)/2 - 60), int(len(dat)/2 + 60)):
            if i % 100 == 0:
                print("*")
            for j in range(int(len(dat[0])/2 - 60), int(len(dat[0])/2 + 60)):
                newdat[i,j] = median(dat[i-2:i+3, j-2:j+3])

        inds = where(isnan(newdat))
        newdat[inds] = 0.

        print(newdat)
        print(newdat.max())

        inds = where(newdat == newdat.max())

        print(inds)

        x = inds[1][0] + 1
        y = inds[0][0] + 1

        name = fl.replace(basedir, '')
        if name.startswith('/'): name=name[1:]
        name = name.split('/')[0]
        name = name.split(':')[0]

        world = w.all_pix2world([[x, y]], 1)

        if name not in outdata.keys():
            outdata[name]=list(world[0])

        f.close()

    return(outdata)

