#!/usr/bin/env python
# Main wrapper for spitzer_sn pipeline
# v0.0 by CDK - 7/10/2022
import sys, os

# See options file to check available command line options
import options

# Import other modules from repository
from analysis import sort_files
from analysis import initial_process
from analysis import subtraction_setup
from analysis import insert_into_subtractions

if __name__ == '__main__':
    args = options.parse_arguments(usage='pipeline.py datadir ra dec [options]')

    # Do file sorting
    sort_files.sort_files(args.datadir,
        channel=args.band, objname=args.object)

    if not os.path.exists(args.mopex_dir):
        print(f'mopex directory {args.mopex_dir} does not exist')
        print('mopex is required to continue')
        print('Install mopex and pass the directory as --mopexdir to continue')
        sys.exit()

    initial_process.initial_process(args.datadir, args.mopex_dir, args.ra,
        args.dec, channel=args.band, nprocesses=args.nprocesses,
        min_exptime=args.min_exptime)

    clobber = not args.no_clobber

    run_file = subtraction_setup.setup_subtractions(args.datadir, args.band,
        args.ra, args.dec, email=args.email, clobber=clobber,
        interactive=args.interactive, date_range=args.date_range,
        offset=args.sn_offset, stamp_size=args.stamp_size)

    # Run new_phot.py script
    basedir=args.datadir
    cmd = f'new_phot.py {run_file} > {basedir}/new_phot_initial.out'
    print(f'Running: {cmd}')
    os.system(cmd)

    run_files = insert_into_subtractions.insert_into_subtractions(args.datadir,
        args.mopex_dir, args.band, args.object, email=args.email)

    # Run insert subtraction files
    for file in sorted(run_files):
        basename = os.path.basename(file)
        logbase = basename.replace('.sh','.log')
        logfile = os.path.join(basedir, logbase)
        cmd = f'/bin/csh {file} > {logfile}'
        print(f'Running: {cmd}')
        os.system(cmd)
