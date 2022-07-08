#!/usr/bin/env python
import argparse
import sys

# Import from repository
import sort_files
import initial_process

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


parser = argparse.ArgumentParser(description='Parse command line options')
parser.add_argument('datadir', type=str, default='.',
                    help='Input directory with r* subdirs from Spitzer')
parser.add_argument('--band', type=str, default='ch1',
                    help='Spitzer/IRAC band to reduce')
parser.add_argument('--object', type=str, default=None,
                    help='Overwrite object name for all files with object')
parser.add_argument('--mopexdir', type=str, default='/data/software/mopex',
                    help='Base directory for mopex')
parser.add_argument('--nprocesses', type=int, default=8,
                    help='Number of threads to use with multiprocessing')

args = parser.parse_args()
args.channel = parse_channel_name(args.channel)

# Do file sorting
sort_files.sort_files(args.datadir, channel=args.band, objname=args.object)

if not os.path.exists(args.mopexdir):
    print(f'mopex directory {args.mopexdir} does not exist')
    print('mopex is required to continue')
    print('Install mopex and pass the directory as --mopexdir to continue')
    sys.exit()

initial_process.initial_process(args.datadir, args.mopexdir,
    channel=args.band, nprocesses=args.nprocesses)
