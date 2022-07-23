import argparse
import sys
from astropy.coordinates import SkyCoord
import numpy as np
from astropy import units as u

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
    parser.add_argument('--object', type=str, default=None,
                        help='Overwrite object name for all files with object')
    parser.add_argument('--mopex-dir', type=str, default='/data/software/mopex',
                        help='Base directory for mopex')
    parser.add_argument('--nprocesses', type=int, default=8,
                        help='Number of threads to use with multiprocessing')
    parser.add_argument('--min-exptime', type=float, default=10.0,
                        help='Minimum exposure time images to analyze')
    parser.add_argument('--no-clobber', default=False, action='store_true',
                        help='Do not clobber previous runs of pipeline')
    parser.add_argument('--email', type=str,
        default='ckilpatrick@northwestern.edu',
        help='Email to message when subtraction process is finished')

    if len(sys.argv) < 4: print(usage) ; sys.exit(1)
    else: coord = parse_coord(sys.argv[2], sys.argv[3])

    if not coord: print(usage) ; sys.exit(1)

    sys.argv[2] = str(coord.ra.degree) ; sys.argv[3] = str(coord.dec.degree)

    args = parser.parse_args()
    args.band = parse_channel_name(args.band)

    return(args)
