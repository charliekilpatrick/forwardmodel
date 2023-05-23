#!/usr/bin/env python
# Main wrapper for spitzer_sn pipeline
# v0.0 by CDK - 7/10/2022
import sys, os, time, glob, shutil

# See options file to check available command line options
# download_spitzer has methods for grabbing all data files from SHA
from common import options
from common import download_spitzer

# Import other modules from repository
from analysis import sort_files
from analysis import initial_process
from analysis import subtraction_setup
from analysis import insert_into_subtractions
from analysis import param_data

if __name__ == '__main__':
    start = time.time()
    args = options.parse_arguments(
        usage='pipeline.py datadir telescope ra dec [options]')

    command = ' '.join(sys.argv)
    options.message(f'Starting: {command}')

    if args.download:
        options.message(f'Downloading to: {args.download}')
        coord = download_spitzer.parse_coord(args.ra, args.dec)
        download_spitzer.download_from_coord(coord, outdir=args.download)

        # Copy directories from download directory into base directory
        for subdir in sorted(glob.glob(os.path.join(args.download, 'r*'))):
            basedirname = os.path.basename(subdir)
            targdir = os.path.join(args.datadir, basedirname)
            if not os.path.exists(targdir):
                print(f'Copying {subdir}->{targdir}')
                shutil.copytree(subdir, targdir)
            else:
                print(f'{targdir} already exists')

    # Do file sorting
    if not args.skip_sort:
        options.message('Sorting files in '+args.datadir)
        sort_files.sort_files(args.datadir, channel=args.band,
            objname=args.object, one_epoch=args.one_epoch)

    if not args.skip_initial_process:
        if not os.path.exists(args.mopex_dir):
            print(f'mopex directory {args.mopex_dir} does not exist')
            print('mopex is required to continue')
            print('Install mopex and use --mopex-dir to continue')
            sys.exit()

        options.message('Initial processing')
        initial_process.initial_process(args.datadir, args.mopex_dir, args.ra,
            args.dec, args.object, channel=args.band,
            nprocesses=args.nprocesses, min_exptime=args.min_exptime,
            date_range=args.all_date_range, use_fif=args.use_fif)

    clobber = not args.no_clobber

    if not args.skip_subtraction:
        options.message('Subtraction setup')
        run_file = subtraction_setup.setup_subtractions(args.datadir, args.band,
            args.ra, args.dec, email=args.email, clobber=clobber,
            interactive=args.interactive, date_range=args.date_range,
            offset=args.sn_offset, stamp_size=args.stamp_size,
            nprocesses=args.nprocesses, prf_version=args.prf_version,
            sci_err_scale=args.sci_err_scale, niter=args.niterations,
            masking=args.masking, mask_radius=args.mask_radius,
            fake_stars=args.fake_stars, fake_radius=args.fake_radius,
            fake_min_mag=args.fake_min_mag, fake_max_mag=args.fake_max_mag)

        # Run new_phot.py script
        options.message('Running photometry script')
        basedir=args.datadir
        #if args.elliptical:
        if True:
            script = os.path.join(param_data.analysis_dir,
                'new_phot_elliptical.py')
        #else:
        #    script = os.path.join(param_data.analysis_dir, 'new_phot.py')
        cmd = f'{script} {run_file}'
        os.chdir(basedir)
        print(cmd)
        os.system(cmd)

        # TODO - finish implementing parallelization for new fm code
        #from analysis import forward_model
        #fm = forward_model.forward_model(run_file)
        #fm.do_main_reduction(fm.parsed, fm.settings)

    if not args.skip_insert_subtractions:
        options.message('Subtracting models from data')
        run_files = insert_into_subtractions.insert_into_subtractions(
            args.datadir, args.mopex_dir, args.band, args.object,
            args.fake_stars, email=args.email)

        # Run insert subtraction files
        for file in sorted(run_files):
            basename = os.path.basename(file)
            logbase = basename.replace('.sh','.log')
            logfile = os.path.join(args.datadir, logbase)
            cmd = f'/bin/csh {file} > {logfile}'
            print(f'Running: {cmd}')
            os.system(cmd)

    total_time = time.time()-start
    message = f'Finished with: {command}\n'
    message += f'It took {total_time} seconds to complete this script.'
    options.message(message)
