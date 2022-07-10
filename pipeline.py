#!/usr/bin/env python
# Main wrapper for spitzer_sn pipeline
# v0.0 by CDK - 7/10/2022
import sys, os

# See options file to check available command line options
import options

# Import other modules from repository
import sort_files
import initial_process

if __name__ == '__main__':
    args = options.parse_arguments(usage='pipeline.py datadir [options]')

    # Do file sorting
    sort_files.sort_files(args.datadir, channel=args.band, objname=args.object)

    if not os.path.exists(args.mopexdir):
        print(f'mopex directory {args.mopexdir} does not exist')
        print('mopex is required to continue')
        print('Install mopex and pass the directory as --mopexdir to continue')
        sys.exit()

    initial_process.initial_process(args.datadir, args.mopexdir,
        channel=args.band, nprocesses=args.nprocesses)
