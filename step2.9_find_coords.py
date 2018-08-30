from numpy import *
import glob
from astropy import wcs
from astropy.io import fits

fls = glob.glob("*:1/working_dir/Combine/mosaic.fits")

print len(fls)

for fl in fls:
    print fl

    f = fits.open(fl)
    w = wcs.WCS(f[0].header)
    dat = f[0].data
    
    newdat = dat*0.

    for i in range(len(dat)/2 - 60, len(dat)/2 + 60):
        if i % 100 == 0:
            print "*"
        for j in range(len(dat[0])/2 - 60, len(dat[0])/2 + 60):
            newdat[i,j] = median(dat[i-2:i+3, j-2:j+3])

    inds = where(isnan(newdat))
    newdat[inds] = 0.

    print newdat
    print newdat.max()

    inds = where(newdat == newdat.max())
    
    print inds

    x = inds[1][0] + 1
    y = inds[0][0] + 1


    
    world = w.all_pix2world([[x, y]], 1)
    print '"' + fl.split(":")[0] + '": ' + str(list(world[0])) + ","
    f.close()

