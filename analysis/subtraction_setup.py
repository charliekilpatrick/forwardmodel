import glob
from astropy.io import fits
from numpy import *
import sys
import tqdm
import os
import shutil
from astropy.time import Time
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils import CircularAperture

from . import find_coords
from . import param_data
from . import inject_fake_stars

def do_it(cmd):
    print(cmd)
    print(os.system(cmd))

def filter0(fls):
    return [item for item in fls if item.count('_0000_0000_') == 0]

max_ims = 1000 # 40 # Limit this to limit RAM or CPU time

def setup_subtractions(basedir, channel, ra, dec,
    email='ckilpatrick@northwestern.edu', clobber=True, interactive=True,
    date_range=[], offset=[0.0, 0.0], stamp_size=29, nprocesses=32,
    prf_version=4, sci_err_scale=20.0, niter=4, masking=False,
    mask_radius=2.0, fake_stars=False, fake_radius=4.0, fake_min_mag=18.0,
    fake_max_mag=24.0):

    if not interactive and not date_range:
        print('ERROR: need to be in interactive mode or provide a date range')
        raise Exception('STOP!')

    drs = glob.glob(os.path.join(basedir, "ut*"))
    drs.sort()

    mjds = []
    for dr in drs:
            print(f'Directory: {dr}')

            expr = os.path.join(dr, f'{channel}/*cbcd*fits')

            first_file = glob.glob(expr)[0]
            print(f'First file: {first_file}')
            f = fits.open(first_file)
            mjds.append(f[0].header['MJD_OBS'])
            f.close()

    for dr in drs:
            subtraction_dir = os.path.join(dr, 'subtraction')
            if os.path.exists(subtraction_dir) and clobber:
                shutil.rmtree(subtraction_dir)
            if not os.path.exists(subtraction_dir):
                os.makedirs(subtraction_dir)

            expr = os.path.join(dr, f'{channel}/*cbcd*.fits')
            cbcds = filter0(glob.glob(expr))
            cbcds.sort()

            for cbcd in cbcds:
                flname = os.path.basename(cbcd)

                f2name = cbcd.replace("cbcd.fits", "cbunc.fits")
                f3name = os.path.join(dr, channel, "Rmask",
                    flname.replace("cbcd.fits", "cbcd_rmask.fits"))

                if not os.path.exists(f2name):
                    print(f'No uncertainty file {f2name}, continuing...')
                    continue
                if not os.path.exists(f3name):
                    print(f'No mask file {f3name}, continuing...')
                    continue

                f1merge = os.path.join(dr, 'subtraction',
                    flname.replace("cbcd.fits", "cbcd_merged.fits"))
                if os.path.exists(f1merge):
                    print(f'{f1merge} exists, continuing...')
                    continue

                f1 = fits.open(cbcd)
                f1.append(f1[0])
                f1[0].data = zeros([0,0], dtype=float32)
                f1[1].name = "SCI"

                f2 = fits.open(cbcd.replace("cbcd.fits", "cbunc.fits"))
                f1.append(f2[0])
                f1[2].name = "ERR"
                f2.close()

                f3 = fits.open(f3name)
                f1.append(f3[0])
                f1[3].name = "DQ"
                f3.close()

                f1.writeto(f1merge, overwrite = True)
                f1.close()



    groups = []
    for i in range(len(drs)):
        if any([i in g['idx'] for g in groups]):
            continue

        group = {}
        group['idx']=[i]
        group['drs']=[drs[i]]
        group['mjds']=[mjds[i]]
        expr = os.path.join(basedir, drs[i], 'subtraction/*merged.fits')
        group['files']=list(glob.glob(expr))
        group['epoch']=None

        for j in range(len(drs)):
            if j != i:
                if mjds[j] > mjds[i] - 20. and mjds[j] < mjds[i] + 20.:
                    group['idx'].append(j)
                    group['drs'].append(drs[j])
                    group['mjds'].append(mjds[j])
                    expr = os.path.join(basedir, drs[j],
                        'subtraction/*merged.fits')
                    group['files'].extend(list(glob.glob(expr)))

        # Sort files by name
        group['files']=sorted(group['files'])

        groups.append(group)

    print('There are {0} groups'.format(len(groups)))

    # Sort in order of ascending MJD
    groups = sorted(groups, key=lambda x: min(x['mjds']))

    # Decide which epoch the groups belong to
    current_epoch = 1
    if interactive:
        for i,group in enumerate(groups):
            min_mjd = min(group['mjds'])
            max_mjd = max(group['mjds'])
            t0 = Time(min_mjd, format='mjd')
            t1 = Time(max_mjd, format='mjd')
            print('Group start time is {0}'.format(
                t0.datetime.strftime('%Y-%m-%dT%H:%M:%S')))
            print('Group end time is {0}'.format(
                t1.datetime.strftime('%Y-%m-%dT%H:%M:%S')))
            x=input('Transient is detectable in this epoch [y/n]? ')
            if 'y' in x:
                # Group is a science epoch
                groups[i]['epoch']=current_epoch
                current_epoch += 1
            else:
                # Assume group is a template epoch
                groups[i]['epoch']=0
    # Assume date range is in MJD
    elif date_range:
        print(f'Date range is {date_range}')
        for i,group in enumerate(groups):
            min_mjd = min(group['mjds'])
            max_mjd = max(group['mjds'])
            t0 = Time(min_mjd, format='mjd')
            t1 = Time(max_mjd, format='mjd')
            print('Group start time is MJD={0}, {1}'.format(min_mjd,
                t0.datetime.strftime('%Y-%m-%dT%H:%M:%S')))
            print('Group end time is MJD={0}, {1}'.format(max_mjd,
                t1.datetime.strftime('%Y-%m-%dT%H:%M:%S')))
            if min_mjd > date_range[0] and max_mjd < date_range[1]:
                # Group is a science epoch
                groups[i]['epoch']=current_epoch
                current_epoch += 1
            else:
                groups[i]['epoch']=0
            epoch = groups[i]['epoch']
            print(f'Epoch is {epoch}')
    else:
        print('ERROR: subtraction_setup assumes interactive=True or '+\
            'date_range is passed.')
        raise Exception('STOP!')

    wd = os.path.join(basedir, 'subtraction')

    lines = param_data.data

    epochs_str = ''
    err_str = ''
    fls_group = []
    for i,group in enumerate(groups):
        fls_group.extend(group['files'])
        err_val = 1.
        epoch_val = group['epoch']
        if epoch_val>0 and not masking:
            err_val = sci_err_scale

        if masking and epoch_val!=0:
            for file in group['files']:
                print(f'Adding mask to {file}')
                hdu = fits.open(file)
                w = wcs.WCS(hdu["DQ"].header)
                snra = float(ra)+float(offset[0])
                sndec = float(dec)+float(offset[1])
                x,y = w.all_world2pix([[float(snra), float(sndec)]], 1)[0]

                # pixel scale in arcsec/pix
                pixscale = sqrt(float(hdu["DQ"].header['CD1_1'])**2 +\
                                   float(hdu["DQ"].header['CD1_2'])**2)
                pixscale = pixscale * 3600.0

                backmask = zeros(hdu["DQ"].data.shape, dtype=bool)
                backmask[int(y-mask_radius/pixscale):int(y+mask_radius/pixscale),
                    int(x-mask_radius/pixscale):int(x+mask_radius/pixscale)]=True

                # Bit #5 is not used for other purposes, so flip this for object
                # Should be 0 if no hit, 32 if hit
                radmask = hdu["DQ"].data & 32
                radmask = radmask!=32

                mask = backmask & radmask

                # Add 8 to all pixels in the mask thus flipping radmask
                hdu["DQ"].data[mask] += 32
                hdu.writeto(file, overwrite=True, output_verify='silentfix')

        # Need to format with + if this is not the first entry
        if not err_str:
            err_str += '[{0}]*{1}'.format(err_val, len(group['files']))
        else:
            err_str += ' + [{0}]*{1}'.format(err_val, len(group['files']))
        if not epochs_str:
            # Update - making this a unique value for every group
            epochs_str += '[{0}]*{1}'.format(i, len(group['files']))
        else:
            epochs_str += ' + [{0}]*{1}'.format(i, len(group['files']))


    psf_names = setup_psf_files(basedir, fls_group, channel, ra, dec,
        prf_version)

    # Inject fake stars into all non-template images
    if fake_stars:
        for i,group in enumerate(groups):
            epoch_val = group['epoch']
            if epoch_val!=0:
                # Create range of fake stars at input radius
                fwhm = 2.02 # This is actually FWHM of Ch2
                num_stars = 2 * pi * fake_radius / (1.5 * 2 * fwhm)
                num_stars = int(floor(num_stars))

                ras = [] ; decs = []
                coord = SkyCoord(ra, dec, unit='deg')

                if num_stars>0:
                    print(f'Injecting {num_stars} into epoch {epoch_val}')
                    mags = random.uniform(fake_min_mag, fake_max_mag, num_stars)

                    # Now get RA/Decs of stars
                    pas = 360.0 / num_stars * arange(num_stars)
                    for pa in pas:
                        newcoord = coord.directional_offset_by(pa*u.deg,
                            fake_radius * u.arcsec)
                        ras.append(newcoord.ra.degree)
                        decs.append(newcoord.dec.degree)

                for i,file in enumerate(group['files']):
                    inject_fake_stars.inject_stars(file, psf_names[i], ras,
                        decs, mags)

    lines = lines.replace("IIIII", str(fls_group))
    lines = lines.replace("PPPPP", epochs_str)
    lines = lines.replace("EEEEE", err_str)

    lines = lines.replace("TTTTT", str(niter))

    lines = lines.replace("RRRRR", str(ra))
    lines = lines.replace("DDDDD", str(dec))

    lines = lines.replace("RRRAOFFSET", str(offset[0]))
    lines = lines.replace("DDECOFFSET", str(offset[1]))

    lines = lines.replace("QQQQQ", str(stamp_size))
    lines = lines.replace("SSSSS", "31")
    lines = lines.replace("AAAAA", "1")

    lines = lines.replace("FFFFF", channel)
    lines = lines.replace("BBBBB", '"'+basedir+'"')

    lines = lines.replace("NNNNN", str(nprocesses))

    lines = lines.replace("VVVVV", str(psf_names))

    # OVRSAMPL should be in fits header
    prf = param_data.get_prf(int(channel.replace('ch','')), prf_version)
    hdu = fits.open(prf)
    if 'OVRSAMPL' not in hdu[0].header.keys():
        raise Exception(f'ERROR: OVRSAMPL not in {prf} header.  Stop!')
    else:
        oversample = hdu[0].header['OVRSAMPL']
        lines = lines.replace("OOOOO", str(oversample))


    param_file = os.path.join(basedir, 'paramfile.txt')
    f = open(param_file, 'w')
    f.write(lines)
    f.close()

    return(param_file)


def setup_psf_files(datadir, images, channel, ra, dec, prf_version):

    prf = param_data.get_prf(int(channel.replace('ch','')), prf_version)
    spatial = param_data.get_spatially_varying(int(channel.replace('ch','')),
        prf_version)

    # Going to make a bunch of new psf files - first make PRF dir
    psf_dir = os.path.join(datadir, 'psfs')
    if not os.path.exists(psf_dir):
        os.makedirs(psf_dir)

    if spatial:

        # Need to do this image by image for given x/y
        dx_file = prf.replace('base','dx')
        dy_file = prf.replace('base','dy')

        psf_hdu = fits.open(prf)
        dx_hdu = fits.open(dx_file)
        dy_hdu = fits.open(dy_file)

        # Make one PSF file for each image
        psf_names = []
        for i,img in enumerate(images):

            hdu = fits.open(img)
            w = wcs.WCS(hdu[0].header)
            pix_xy = w.all_world2pix([[float(ra), float(dec)]], 1)[0]
            pix_xy = array(around(pix_xy), dtype=int32)

            psf_data = psf_hdu[0].data + (pix_xy[0]-128)*dx_hdu[0].data +\
                (pix_xy[1]-128)*dy_hdu[0].data

            img_base = os.path.basename(img)
            psf_name = img_base.replace('.fits','_psf.fits')
            psf_name = os.path.join(psf_dir, psf_name)

            newhdu = fits.PrimaryHDU()
            newhdu.header = psf_hdu[0].header
            newhdu.data = psf_data

            newhdu.writeto(psf_name, overwrite=True,
                output_verify='silentfix')
            psf_names.append(psf_name)

    else:
        psf_names = [prf]*len(images)

    return(psf_names)
