import argparse
import sys

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
    parser.add_argument('--band','--channel', type=str, default='ch1',
                        help='Spitzer/IRAC band to reduce')
    parser.add_argument('--object', type=str, default=None,
                        help='Overwrite object name for all files with object')
    parser.add_argument('--mopexdir', type=str, default='/data/software/mopex',
                        help='Base directory for mopex')
    parser.add_argument('--nprocesses', type=int, default=8,
                        help='Number of threads to use with multiprocessing')

    args = parser.parse_args()
    args.band = parse_channel_name(args.band)

    return(args)
