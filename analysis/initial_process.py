import os
import glob
import multiprocessing as mp
import subprocess
import shutil
import copy
import time
import astropy
from astropy.io import fits
from astropy import wcs
from astropy.time import Time

def create_fif(ra, dec, imsize=750, pixscale=0.000166670,
    filename='mosaic_fif.tbl'):

    crpix = float(imsize)/2.0 + 0.5

    with open(filename, 'w') as outfile:
        outfile.write('\\char comment = Custom FIF by CDK \n')
        outfile.write('\\char Date-Time = '+time.asctime()+'\n')
        outfile.write('\\char PROJTYPE = \'TAN\' \n')
        outfile.write(f'\\real  CRVAL1 = {ra} \n')
        outfile.write(f'\\real  CRVAL2 = {dec} \n')
        outfile.write(f'\\real  CRPIX1 = {crpix} \n')
        outfile.write(f'\\real  CRPIX2 = {crpix} \n')
        outfile.write(f'\\real  CROTA2 = 0.0 \n')
        outfile.write(f'\\int   NAXIS1 = {imsize} \n')
        outfile.write(f'\\int   NAXIS2 = {imsize} \n')
        outfile.write(f'\\real  CDELT1 = -{pixscale} \n')
        outfile.write(f'\\real  CDELT2 =  {pixscale} \n')
        outfile.write(f'\\char  CTYPE1 =  \"RA---TAN\" \n')
        outfile.write(f'\\char  CTYPE2 =  \"DEC--TAN\" \n')

def filter0(list_of_fls):
    return [item for item in list_of_fls if item.count('_0000_0000_') == 0]

def write_out_inlist(fls, outfile):
    with open(outfile, 'w') as f:
        f.write('\n'.join(fls)+'\n')

def run_dir(var):
    dr, mopex, ra, dec, channel, min_exptime,\
        min_mjd, max_mjd, objname, use_fif = var

    ra = float(ra) ; dec = float(dec) ; min_exptime = float(min_exptime)
    min_mjd = float(min_mjd) ; max_mjd = float(max_mjd)

    use_fif = bool(use_fif)

    wd = os.path.abspath(os.path.join(dr, channel))

    if not os.path.exists(wd):
        os.makedirs(wd)

    # Delete unnecessary files
    for file in glob.glob(os.path.join(dr, channel, '*')):
        if (('cbcd.fits' not in file) and ('cbunc.fits' not in file) and
            ('bimsk.fits' not in file)):
            os.remove(file)

    fls = filter0(glob.glob(os.path.join(dr, channel, '*cbcd.fits')))
    fls.sort()

    new_fls = []
    for fl in fls:
        hdu = fits.open(fl, mode='readonly')
        mjd = hdu[0].header['MJD_OBS']

        if mjd < min_mjd or mjd > max_mjd:
            print(f'Skipping {fl}: outside MJD {min_mjd}->{max_mjd}')
            continue

        if min_exptime:
            exptime=hdu[0].header['EXPTIME']
            if exptime<min_exptime:
                print(f'Skipping {fl}: EXPTIME={exptime}<{min_exptime}')
                continue

        w = wcs.WCS(hdu[0].header)
        try:
            x,y = w.all_world2pix([[ra, dec]], 1)[0]
        except astropy.wcs.wcs.NoConvergence:
            print(f'Skipping {fl}: ra={ra}, dec={dec} not in image')
            continue

        if (x>hdu[0].header['NAXIS1'] or x<0 or
            y>hdu[0].header['NAXIS2'] or y<0):
            print(f'Skipping {fl}: ra={ra}, dec={dec} not in image')
            continue

        new_fls.append(fl)

    fls = copy.copy(new_fls)
    nfls = len(fls)

    print(f'We need to reduce {nfls} files for {dr}')

    first_file = fls[0]
    mjd = fits.open(first_file)[0].header['MJD_OBS']
    aorkey = str(fits.open(first_file)[0].header['AORKEY'])
    t = Time(mjd, format='mjd')
    datestr = t.datetime.strftime('ut%y%m%d')

    baseoutname = f'{objname}.{channel}.{datestr}.{aorkey}_stk.fits'
    baseoutcov = baseoutname.replace('_stk.fits','_cov.fits')
    baseoutunc = baseoutname.replace('_stk.fits','_unc.fits')

    ufls = [fl.replace("cbcd.fits", "cbunc.fits") for fl in fls]
    mfls = [fl.replace("cbcd.fits", "bimsk.fits") for fl in fls]

    if not os.path.exists(os.path.join(wd, 'input')):
        os.makedirs(os.path.join(wd, 'input'))

    write_out_inlist(fls, os.path.join(wd,  'input/images.list'))
    write_out_inlist(ufls, os.path.join(wd, 'input/sigma.list'))
    write_out_inlist(mfls, os.path.join(wd, 'input/mask.list'))

    os.system(f"cp -r {mopex}/cal {wd}")
    os.system(f"cp -r {mopex}/cdf {wd}")
    os.system(f"cp {mopex}/mopex-script-env.csh {wd}")

    inp = f'-I input/images.list -S input/sigma.list -d input/mask.list -O {wd}'

    if use_fif:
        fif_file = os.path.join(wd, 'mosaic_fif.tbl')
        create_fif(ra, dec, imsize=750, pixscale=0.000166670,
            filename=fif_file)
        inp += f' -F {fif_file}'

    with open(os.path.join(wd, 'run.sh'), 'w') as f:
        ch = channel.replace('ch','I')
        f.write(f'cd {wd} \n')
        f.write('source mopex-script-env.csh \n')
        if use_fif:
            f.write(f'overlap.pl -n overlap_{ch}_nofid.nl {inp} > overlap.log \n')
            f.write(f'mosaic.pl -n mosaic_{ch}_nofid.nl {inp} > mosaic.log \n')
        else:
            f.write(f'overlap.pl -n overlap_{ch}.nl {inp} > overlap.log \n')
            f.write(f'mosaic.pl -n mosaic_{ch}.nl {inp} > mosaic.log \n')

    cmd = f"/bin/csh {wd}/run.sh"
    print(cmd)
    print(subprocess.call(cmd, shell=True))

    for subdr in ["BoxOutlier", "ReInterp", "Overlap_Corr", "Detect",
        "DualOutlier", "Outlier", "Interp", "Medfilter", "Coadd",
        "cal", "cdf", "addkeyword.txt", "*.nl",
        "FIF.tbl", "header_list.tbl"]:
        cmd = f"rm -fvr {os.path.join(wd, subdr)}"
        print(cmd)
        os.system(cmd)

    cmd1=f"mv -v {wd}/Combine/mosaic.fits {wd}/{baseoutname}"
    cmd2=f"mv -v {wd}/Combine/mosaic_cov.fits {wd}/{baseoutcov}"
    cmd3=f"mv -v {wd}/Combine/mosaic_unc.fits {wd}/{baseoutunc}"
    cmd4=f"rm -fvr {wd}/Combine"

    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)


def initial_process(basedir, mopex, ra, dec, objname, channel='ch1',
    min_exptime=None, date_range=[], nprocesses=8, use_fif=False):

    pool = mp.Pool(processes=nprocesses)

    drs = glob.glob(os.path.join(basedir, "ut*"))
    if len(date_range)==2:
        min_mjd = date_range[0]
        max_mjd = date_range[1]
    else:
        min_mjd = -999999
        max_mjd = 999999
    var = [(dr, mopex, str(ra), str(dec), channel, str(min_exptime),
        str(min_mjd), str(max_mjd), str(objname), str(use_fif)) for dr in drs]
    pool.map(run_dir, var)

