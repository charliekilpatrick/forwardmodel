import os
import glob
import multiprocessing as mp
import subprocess
import shutil
import copy
import time
import astropy
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.time import Time


def format_header(img):

    hdu = fits.open(img)
    hdr = hdu[0].header

    hdu[0].header['MJD-OBS'] = Time(hdr['DATE_OBS']).mjd 
    hdu[0].header['GAIN'] = 1.0
    hdu[0].header['RDNOISE'] = 0.0

    if 'stk.fits' in img:
        hdu[0].header['OBSTYPE']='OBJECT'
    elif 'stk.noise.fits' in img:
        hdu[0].header['OBSTYPE']='NOISE'
    elif 'stk.mask.fits' in img:
        hdu[0].header['OBSTYPE']='MASK'
    elif 'stk.cov.fits' in img:
        hdu[0].header['OBSTYPE']='COVERAGE'

    if 'ch1' in img:
        #https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/9/
        # Assume 2.8 mJy saturation level
        hdu[0].header['SATURATE']=330.905902179
    elif 'ch2' in img:
        hdu[0].header['SATURATE']=342.723970114
    elif 'ch3' in img:
        hdu[0].header['SATURATE']=3190.87834245
    elif 'ch4' in img:
        hdu[0].header['SATURATE']=3427.23970114

    hdu.writeto(img, overwrite=True, output_verify='silentfix')

def rescale_img(img, target_zpt=27.5, zpt_err=0.01):

    hdu = fits.open(img)
    hdr = hdu[0].header

    if 'stk.fits' in img:
        hdu[0].data = hdu[0].data - np.min(hdu[0].data)


    w = wcs.WCS(hdr)

    pscale = np.mean(wcs.utils.proj_plane_pixel_scales(w)) * 3600.0

    zpt=-2.5*np.log10(pscale*pscale*23.5045*1.0e-6/3631.0)

    rescale = 10**(0.4*(target_zpt - zpt))

    hdu[0].header['ZPTMAG']=target_zpt

    # Conservatively estimate 1% uncertainty:
    # https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/52/
    hdu[0].header['ZPTMUNC']=zpt_err

    hdu[0].data = hdu[0].data * rescale
    # Rescale saturate
    hdu[0].header['SATURATE']=hdr['SATURATE']*rescale

    hdu.writeto(img, overwrite=True, output_verify='silentfix')


def sanitize_and_create_mask(img):

    hdu = fits.open(img)
    data = hdu[0].data

    # Get parts of image that are nan
    mask = np.isnan(data)

    # Set pixels to 129 in mask for nan values
    newdata = np.zeros(data.shape)
    newdata[mask]=129

    # Reset image to 0 where it was previously masked to nan
    hdu[0].data[mask] = 0

    newdata = newdata.astype('uint16')
    newhdu = fits.PrimaryHDU()

    newhdu.data = newdata
    newhdu.header = hdu[0].header
    newhdu.header['BITPIX']=16

    maskfile = img.replace('.fits','.mask.fits')

    newhdu.writeto(maskfile, overwrite=True, output_verify='silentfix')
    hdu.writeto(img, overwrite=True, output_verify='silentfix')


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
        min_mjd, max_mjd, objname, use_fif, one_epoch, basedir, photpipe = var

    ra = float(ra) ; dec = float(dec) ; min_exptime = float(min_exptime)
    min_mjd = float(min_mjd) ; max_mjd = float(max_mjd)

    use_fif = str(use_fif)=='True'
    one_epoch = str(one_epoch)=='True'

    if photpipe=='None':
        photpipe = None

    wd = os.path.abspath(os.path.join(dr, channel))

    if not os.path.exists(wd):
        os.makedirs(wd)

    # Delete unnecessary files
    for file in glob.glob(os.path.join(dr, channel, '*')):
        if (('cbcd.fits' not in file) and ('cbunc.fits' not in file) and
            ('bimsk.fits' not in file)):
            if os.path.isdir(file):
                shutil.rmtree(file)
            else:
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

    if nfls==0:
        print(f'WARNING: no files to reduce.  Skipping reduction of {dr}')
        print(f'Deleting {dr}')
        shutil.rmtree(dr)
        return

    first_file = fls[0]
    mjd = fits.open(first_file)[0].header['MJD_OBS']
    aorkey = str(fits.open(first_file)[0].header['AORKEY'])
    t = Time(mjd, format='mjd')

    if one_epoch:
        datestr = 'all'
    else:
        datestr = t.datetime.strftime('ut%y%m%d')

    baseoutname = f'{objname}.{channel}.{datestr}.{aorkey}_stk.fits'
    baseoutcov = baseoutname.replace('.fits','.cov.fits')
    baseoutunc = baseoutname.replace('.fits','.noise.fits')
    baseoutmsk = baseoutname.replace('.fits','.mask.fits')

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

    fmt_files = [f'{basedir}/{baseoutname}',f'{basedir}/{baseoutunc}',
        f'{basedir}/{baseoutcov}',f'{basedir}/{baseoutmsk}']

    cmd1=f"mv -v {wd}/Combine/mosaic.fits {basedir}/{baseoutname}"
    cmd2=f"mv -v {wd}/Combine/mosaic_cov.fits {basedir}/{baseoutcov}"
    cmd3=f"mv -v {wd}/Combine/mosaic_unc.fits {basedir}/{baseoutunc}"
    cmd4=f"rm -fvr {wd}/Combine"

    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)

    outimg = f'{basedir}/{baseoutname}'
    if os.path.exists(outimg):
        sanitize_and_create_mask(outimg)

    for file in [f'{basedir}/{baseoutname}',f'{basedir}/{baseoutunc}',
        f'{basedir}/{baseoutcov}',f'{basedir}/{baseoutmsk}']:
        if os.path.exists(file):
            format_header(file)

    for file in [f'{basedir}/{baseoutname}',f'{basedir}/{baseoutunc}']:
        if os.path.exists(file):
            rescale_img(file)

    if photpipe:
        rawdir = os.path.join(photpipe, 'rawdata', objname, channel)
        if not os.path.exists(rawdir):
            os.makedirs(rawdir)

        outbaseimg = os.path.join(rawdir, baseoutname)
        shutil.copyfile(f'{basedir}/{baseoutname}', outbaseimg)

        workdir = os.path.join(photpipe, 'workspace', objname, channel)
        if not os.path.exists(workdir):
            os.makedirs(workdir)

        outbaseunc = os.path.join(workdir, baseoutunc)
        if os.path.exists(f'{basedir}/{baseoutunc}'):
            shutil.copyfile(f'{basedir}/{baseoutunc}', outbaseunc)

        outbasemsk = os.path.join(workdir, baseoutmsk)
        if os.path.exists(f'{basedir}/{baseoutmsk}'):
            shutil.copyfile(f'{basedir}/{baseoutmsk}', outbasemsk)

def initial_process(basedir, mopex, ra, dec, objname, channel='ch1',
    min_exptime=None, date_range=[], nprocesses=8, use_fif=False,
    one_epoch=False, photpipe=None):

    pool = mp.Pool(processes=nprocesses)

    drs = glob.glob(os.path.join(basedir, "ut*"))
    if len(date_range)==2:
        min_mjd = date_range[0]
        max_mjd = date_range[1]
    else:
        min_mjd = -999999
        max_mjd = 999999
    var = [(dr, mopex, str(ra), str(dec), channel, str(min_exptime),
        str(min_mjd), str(max_mjd), str(objname), str(use_fif),
        str(one_epoch), str(basedir), str(photpipe)) for dr in drs]
    pool.map(run_dir, var)

