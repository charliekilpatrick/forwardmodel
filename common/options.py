import argparse
import sys
from astropy.coordinates import SkyCoord
import numpy as np
from astropy import units as u
from astropy.time import Time

def message(msg):
    print('\n\n'+msg+'\n'+'#'*80+'\n'+'#'*80+'\n\n')

def parse_two_floats(value):
    values = value.split()
    if len(values) != 2:
        raise argparse.ArgumentError
    values = map(float, values)
    return(values)

def is_number(num):
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)

def parse_coord(ra, dec):
    if (not (is_number(ra) and is_number(dec)) and
        (':' not in ra and ':' not in dec)):
        error = 'ERROR: cannot interpret: {ra} {dec}'
        print(error.format(ra=ra, dec=dec))
        return(None)

    if (':' in ra and ':' in dec):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        error = 'ERROR: Cannot parse coordinates: {ra} {dec}'
        print(error.format(ra=ra,dec=dec))
        return(None)

# Channel name needs to be ch1, ch2, ch3, ch4
def parse_channel_name(channel):
    if not channel:
        print(f'ERROR: channel={channel} cannot be parsed')
        sys.exit()
    if '1' in str(channel):
        return('ch1')
    elif '2' in str(channel):
        return('ch2')
    elif '3' in str(channel):
        return('ch3')
    elif '4' in str(channel):
        return('ch4')

    print(f'ERROR: channel={channel} cannot be parsed')
    sys.exit()

def parse_date_range(date_range_data):
    t0 = None
    date_range = []
    if is_number(date_range_data[0]):
        # Assume init date is MJD
        try:
            t0 = Time(date_range_data[0], format='mjd')
        except:
            pass
    else:
        try:
            t0 = Time(date_range_data[0])
        except:
            pass
    reldate = None ; t1 = None
    if is_number(date_range_data[1]):
        if float(date_range_data[1])<10000:
            reldate = float(date_range_data[1])
        else:
            try:
                t1 = Time(date_range_data[1], format='mjd')
            except:
                pass
    else:
        try:
            t1 = Time(date_range_data[1])
        except:
            pass

    if t0 and reldate:
        if reldate < 0:
            date_range = [t0.mjd + reldate, t0.mjd]
            print(f'Interactive=False, date range={args.date_range}')
        else:
            date_range = [t0.mjd, t0.mjd + reldate]
    elif t0 and t1:
        if t0 < t1:
            date_range = [t0.mjd, t1.mjd]
        else:
            date_range = [t1.mjd, t0.mjd]

    return(date_range)

def parse_arguments(usage=''):

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('datadir', type=str, default='.',
                        help='Input directory with r* subdirs from Spitzer')
    parser.add_argument('ra', type=str, default='0.0',
                        help='Right ascension of region to analyze in images')
    parser.add_argument('dec', type=str, default='0.0',
                        help='Declination of region to analyze in images')
    parser.add_argument('--band','--channel', type=str, default='ch1',
                        help='Spitzer/IRAC band to reduce')
    parser.add_argument('--download', type=str, default=None,
                        help='Download files for input RA/Dec to input dir.')
    parser.add_argument('--object', type=str, default=None,
                        help='Overwrite object name for all files with object')
    parser.add_argument('--mopex-dir', type=str, default='/data/software/mopex',
                        help='Base directory for mopex')
    parser.add_argument('--nprocesses', type=int, default=4,
                        help='Number of threads to use with multiprocessing')
    parser.add_argument('--min-exptime', type=float, default=10.0,
                        help='Minimum exposure time images to analyze')
    parser.add_argument('--no-clobber', default=False, action='store_true',
                        help='Do not clobber previous runs of pipeline')
    parser.add_argument('--science-date-range', default=None, type=str, nargs=2,
                        help='UT date (MJD, ISO, etc.) for date range of '+\
                        'range of science images.')
    parser.add_argument('--data-date-range', default=None, type=str, nargs=2,
                        help='UT date (MJD, ISO, etc.) for date range of '+\
                        'all data.')
    parser.add_argument('--stamp-size', default=29, type=float,
                        help='Stamp size for analyzing Spitzer data '+\
                        '(must be odd integer).')
    parser.add_argument('--prf-version', default=2, type=int,
                        help='Version of the PRF to use from data directory.')
    parser.add_argument('--sn-offset', action='store', type=float,
                  help='Offsets (in RA/Dec) between input RA/Dec and target',
                  default=[0.0, 0.0], nargs=2)
    parser.add_argument('--sci-err-scale', default=20.0, type=float,
                        help='Scale factor to use for error scaling of model '+\
                        'in science images (i.e., with transient emission).')
    parser.add_argument('--niterations', default=4, type=int,
                        help='number of iterations for forward model code')
    parser.add_argument('--cluster', default=False, action='store_true',
                        help='Run forward_model in cluster mode')
    parser.add_argument('--skip-sort', default=False, action='store_true',
                        help='Skip initial file sorting')
    parser.add_argument('--skip-initial-process', default=False,
        action='store_true', help='Skip initial file processing')
    parser.add_argument('--skip-subtraction', default=False,
        action='store_true', help='Skip subtraction setup and forward model')
    parser.add_argument('--skip-insert-subtractions', default=False,
        action='store_true', help='Skip final mosaicing with subtractions')
    parser.add_argument('--use-fif', default=False,
        action='store_true', help='Use a FIF table generated from ra/dec')
    parser.add_argument('--masking', default=False,
        action='store_true', help='Use masking instead of error scaling for '+\
        'science epochs.')
    parser.add_argument('--mask-radius', default=2.0, type=float,
        help='Radius (in arcsec) of masking.')
    parser.add_argument('--fake-stars', default=False,
        action='store_true', help='Inject fake stars into science images.')
    parser.add_argument('--fake-radius', default=4.0, type=float,
        help='Radius (in arcsec) to inject fake stars around RA/Dec.')
    parser.add_argument('--fake-min-mag', default=18.0, type=float,
        help='Minimum magnitude of injected fake stars.')
    parser.add_argument('--fake-max-mag', default=21.0, type=float,
        help='Maximum magnitude of injected fake stars.')
    parser.add_argument('--redo-mosaic', default=False,
        action='store_true', help='Redo the mosaic during subtraction setup.')
    parser.add_argument('--email', type=str,
        default='ckilpatrick@northwestern.edu',
        help='Email to message when subtraction process is finished')
    parser.add_argument('--elliptical', default=False,
        action='store_true', help='Using 1D elliptical model forward model.')
    parser.add_argument('--one-epoch', default=False,
        action='store_true', help='Reduce all data as a single epoch.')

    if len(sys.argv) < 4: print(usage) ; sys.exit(1)
    else: coord = parse_coord(sys.argv[2], sys.argv[3])

    if not coord: print(usage) ; sys.exit(1)

    sys.argv[2] = str(coord.ra.degree) ; sys.argv[3] = str(coord.dec.degree)

    args = parser.parse_args()
    args.band = parse_channel_name(args.band)

    args.date_range = []
    args.all_date_range = []
    args.interactive = True

    args.stamp_size = int(args.stamp_size)
    if args.stamp_size % 2 == 0:
        print('WARNING: stamp size must be an odd integer')
        stamp_size = args.stamp_size + 1
        print(f'Setting stamp size = {stamp_size}')
        args.stamp_size = int(stamp_size)

    if args.science_date_range:
        date_range = parse_date_range(args.science_date_range)
        if len(date_range)==2:
            args.interactive = False
            args.date_range = date_range

    if args.data_date_range:
        date_range = parse_date_range(args.science_date_range)
        if len(date_range)==2:
            args.all_date_range = date_range

    return(args)
