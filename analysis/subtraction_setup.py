import glob
from astropy.io import fits
from numpy import *
import sys
import tqdm
import os
import shutil
from astropy.time import Time

from . import find_coords
from . import param_data

def do_it(cmd):
    print(cmd)
    print(os.system(cmd))

def filter0(fls):
    return [item for item in fls if item.count('_0000_0000_') == 0]

max_ims = 1000 # 40 # Limit this to limit RAM or CPU time

def setup_subtractions(basedir, channel, ra, dec,
    email='ckilpatrick@northwestern.edu', clobber=True, interactive=True,
    date_range=[], offset=[0.0, 0.0], stamp_size=29):

    if not interactive and not date_range:
        print('ERROR: need to be in interactive mode or provide a date range')
        raise Exception('STOP!')

    drs = glob.glob(os.path.join(basedir, "*:*"))
    drs.sort()

    mjds = []

    for dr in drs:
            print(f'Directory: {dr}')

            expr = os.path.join(dr, f'{channel}/bcd/*cbcd*fits')

            print(expr)

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

            expr = os.path.join(dr, f'{channel}/bcd/*cbcd*.fits')
            cbcds = filter0(glob.glob(expr))
            cbcds.sort()

            for cbcd in cbcds:
                flname = cbcd.split("/")[-1]

                f2name = cbcd.replace("cbcd.fits", "cbunc.fits")
                f3name = os.path.join(dr, 'working_dir/Rmask',
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

                print(f1.info())

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
    groups = sorted(groups, key=lambda x: argmin(x['mjds']))

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
        for i,group in enumerate(groups):
            min_mjd = min(group['mjds'])
            max_mjd = max(group['mjds'])
            if min_mjd > date_range[0] and max_mjd < date_range[1]:
                # Group is a science epoch
                groups[i]['epoch']=current_epoch
                current_epoch += 1
            else:
                groups[i]['epoch']=0
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
        if epoch_val>0:
            err_val = 20.

        # Need to format with + if this is not the first entry
        if not err_str:
            err_str += '[{0}]*{1}'.format(err_val, len(group['files']))
        else:
            err_str += ' + [{0}]*{1}'.format(err_val, len(group['files']))
        if not epochs_str:
            epochs_str += '[{0}]*{1}'.format(epoch_val, len(group['files']))
        else:
            epochs_str += ' + [{0}]*{1}'.format(epoch_val, len(group['files']))

    lines = lines.replace("IIIII", str(fls_group))
    lines = lines.replace("PPPPP", epochs_str)
    lines = lines.replace("EEEEE", err_str)

    lines = lines.replace("RRRRR", str(ra))
    lines = lines.replace("DDDDD", str(dec))

    lines = lines.replace("RRRAOFFSET", str(offset[0]))
    lines = lines.replace("DDECOFFSET", str(offset[1]))

    lines = lines.replace("QQQQQ", str(stamp_size))
    lines = lines.replace("SSSSS", "31")
    lines = lines.replace("AAAAA", "1")

    lines = lines.replace("FFFFF", channel)

    param_file = os.path.join(basedir, 'paramfile.txt')
    f = open(param_file, 'w')
    f.write(lines)
    f.close()

    return(param_file)

