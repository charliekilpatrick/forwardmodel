import glob
from astropy.io import fits
from numpy import *
import commands

drs = glob.glob("r*")

assert len(drs) == 320

objs = []
mjds = []

for dr in drs:
    print dr

    fls = glob.glob(dr + "/ch1/bcd/*_cbcd.fits")

    print fls[0]
    
    f = fits.open(fls[0])
    obj = f[0].header["OBJECT"]
    mjd = f[0].header["MJD_OBS"]
    f.close()
    
    mjds.append(mjd)
    objs.append(obj)

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
    cmd = "mv " + drs[i] + " " + objs[i].replace(" ", "_") + ":" + str(epochs[objs[i]])
    print cmd
    print commands.getoutput(cmd)

    epochs[objs[i]] += 1
