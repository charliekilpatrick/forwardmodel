import os
import glob
import multiprocessing as mp
import subprocess
import shutil
import copy
from astropy.io import fits
from astropy import wcs

def filter0(list_of_fls):
    return [item for item in list_of_fls if item.count('_0000_0000_') == 0]

def write_out_inlist(fls, outfile):
    with open(outfile, 'w') as f:
        f.write('\n'.join(fls)+'\n')

def run_dir(var):
    dr, mopex, ra, dec, channel, min_exptime = var

    ra = float(ra) ; dec = float(dec) ; min_exptime = float(min_exptime)

    wd = os.path.join(dr, 'working_dir')

    if os.path.exists(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)

    # Delete unnecessary files
    for file in glob.glob(os.path.join(dr, channel, 'bcd/*')):
        if (('cbcd.fits' not in file) and ('cbunc.fits' not in file) and
            ('bimsk.fits' not in file)):
            os.remove(file)

    fls = filter0(glob.glob(os.path.join(dr, channel, 'bcd/*cbcd.fits')))
    fls.sort()

    new_fls = []
    for fl in fls:
        hdu = fits.open(fl, mode='readonly')
        if min_exptime:
            exptime=hdu[0].header['EXPTIME']
            if exptime<min_exptime:
                print(f'Skipping {fl}: EXPTIME={exptime}<{min_exptime}')
                continue
        w = wcs.WCS(hdu[0].header)
        x,y = w.all_world2pix([[ra, dec]], 1)[0]

        if (x>hdu[0].header['NAXIS1'] or x<0 or
            y>hdu[0].header['NAXIS2'] or y<0):
            print(f'Skipping {fl}: ra={ra}, dec={dec} not in image')
            continue

        new_fls.append(fl)

    fls = copy.copy(new_fls)
    nfls = len(fls)

    print(f'We need to reduce {nfls} files for {dr}')

    ufls = [fl.replace("cbcd.fits", "cbunc.fits") for fl in fls]
    mfls = [fl.replace("cbcd.fits", "bimsk.fits") for fl in fls]

    if not os.path.exists(os.path.join(wd, 'input')):
        os.makedirs(os.path.join(wd, 'input'))

    write_out_inlist(fls, os.path.join(wd, 'input/images.list'))
    write_out_inlist(ufls, os.path.join(wd, 'input/sigma.list'))
    write_out_inlist(mfls, os.path.join(wd, 'input/mask.list'))

    os.system(f"cp -r {mopex}/cal {wd}")
    os.system(f"cp -r {mopex}/cdf {wd}")
    os.system(f"cp {mopex}/mopex-script-env.csh {wd}")

    inp = '-I input/images.list -S input/sigma.list -d input/mask.list'

    with open(os.path.join(wd, 'run.sh'), 'w') as f:
        ch = channel.replace('ch','I')
        f.write(f'cd {wd} \n')
        f.write('source mopex-script-env.csh \n')
        f.write(f'overlap.pl -n overlap_{ch}.nl {inp} > overlap.log \n')
        f.write(f'mosaic.pl -n mosaic_{ch}.nl {inp} > mosaic.log \n')

    cmd = f"/bin/csh {wd}/run.sh"
    print(cmd)
    print(subprocess.call(cmd, shell=True))

    for subdr in ["BoxOutlier", "ReInterp", "Overlap_Corr", "Detect",
        "DualOutlier", "Outlier", "Interp", "Medfilter", "Coadd",
        "Combine", "cal", "cdf", "addkeyword.txt", "*.nl",
        "*.fits", "FIF.tbl", "header_list.tbl"]:
        cmd = f"rm -fvr {os.path.join(wd, subdr)}"
        print(cmd)
        os.system(cmd)

def initial_process(basedir, mopex, ra, dec, channel='ch1', min_exptime=None,
    nprocesses=8):

    pool = mp.Pool(processes=nprocesses)

    drs = glob.glob(os.path.join(basedir, "*:*"))
    var = [(dr, mopex, str(ra), str(dec), channel, str(min_exptime))
        for dr in drs]
    pool.map(run_dir, var)

