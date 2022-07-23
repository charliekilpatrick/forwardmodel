#!/usr/bin/env python
from astropy.io import fits
from numpy import *
import gzip
import pickle
import sys
import tqdm
import os
import glob
import shutil

def create_mopex_cmd(basecmd, idx, channel):

    cmd = ''
    if '.pl' not in basecmd:
        cmd = basecmd
        cmd += '.pl'
    else:
        cmd = basecmd
        basecmd = cmd.replace('.pl','')

    ch = channel.replace('ch','I')

    if idx==0:
        cmd += f' -n {basecmd}_{ch}.nl'
    else:
        cmd += f' -n {basecmd}_{ch}_nofid.nl'

    cmd += ' -I images.list -S sigma.list -d mask.list'

    if idx!=0:
        cmd += ' -F mosaic_fif.tbl'

    cmd += f' > log_{basecmd}.txt'

    return(cmd)

def create_nofid_params(mopex, basecmd, channel):

    ch = channel.replace('ch','I')

    outparamfile = os.path.join(mopex, 'cdf', f'{basecmd}_{ch}_nofid.pl')
    if os.path.exists(outparamfile):
        return(None)

    paramfile = outparamfile.replace('_nofid','')

    with open(paramfile, 'r') as params:
        data = params.readlines()

    newdata = data.replace('run_fiducial_image_frame = 1',
        'run_fiducial_image_frame = 0')

    with open(outparamfile, 'w') as params:
        params.write(newdata)


def insert_into_subtractions(basedir, mopex, channel, objname,
    email='ckilpatrick@northwestern.edu'):

    out_runfiles=[]

    match = os.path.join(basedir, 'run_stacks_*.sh')
    for file in glob.glob(match):
        if os.path.exists(file):
            print(f'Deleting {file}')
            os.remove(file)

    pkl_file = os.path.join(basedir, 'fit_results.pickle')
    pkl_data = gzip.open(pkl_file, 'rb')
    [all_data, parsed, settings, SNCmat, Cmat] = pickle.load(pkl_data)

    for key in settings:
        print(f"settings: {key}")

    for key in all_data:
        print(f"all_data: {key}")

    print('epoch_names:',settings['epoch_names'])
    print('epochs:',settings['epochs'])

    first_epoch_dir = ''
    for ee in sorted(settings['epoch_names']):

        wd = os.path.join(basedir, "sub_stack_epoch{0}".format(
            str(ee).zfill(3)))
        wd_im = os.path.join(wd, "ims")

        if os.path.exists(wd):
            print(f'Deleting: {wd}')
            shutil.rmtree(wd)

        os.makedirs(wd)
        os.makedirs(wd_im)

        filename = os.path.join(basedir, "run_stacks_{0}.sh".format(
            str(ee).zfill(3)))
        f_mopex = open(filename, 'w')
        f_mopex.write("sleep 2 \n")

        print("Figuring out which images to look at...")

        images_to_work_with = []
        bad_images = []
        for idx in where(settings['epochs']==ee)[0]:
            if (settings["images"][idx].count(basedir) and
                os.path.exists(settings["images"][idx])):
                images_to_work_with.append(idx)
            else:
                bad_images.append(idx)

        assert len(images_to_work_with) > 0
        n=len(images_to_work_with)
        print(f"Found {n} images to work with")
        if len(bad_images)>0:
            print(f"Bad images: {bad_images}")
        else:
            print("Found no bad images")

        resid_file = os.path.join(basedir, 'residuals.fits')

        f = fits.open(resid_file)
        subtractions = f[0].data
        print(f"subtractions shape: {subtractions.shape}")
        f.close()

        f_ilist = open(os.path.join(wd,"images.list"), 'w')
        f_slist = open(os.path.join(wd,"sigma.list"), 'w')
        f_mlist = open(os.path.join(wd,"mask.list"), 'w')

        total_im = len(images_to_work_with)
        print(f'Writing out {total_im} images')
        for imind in images_to_work_with:
            # Construct full path to new image
            newim_base = os.path.basename(settings["images"][imind])
            newim_base = newim_base.replace("cbcd_merged.fits", "cbcd.fits")
            newim = os.path.join(wd_im, newim_base)

            assert newim != settings["images"][imind]

            origim_base = settings["images"][imind].replace("cbcd_merged.fits",
                "cbcd.fits")
            origim_base = os.path.basename(origim_base)
            origim_path = os.path.split(settings["images"][imind])[0]

            # Original image path should be like...
            origim_path = origim_path.replace('subtraction', channel+'/bcd')

            # Full path to original image
            origim = os.path.join(origim_path, origim_base)

            assert os.path.exists(origim)

            f = fits.open(origim)

            pixels_not_modified_by_subtraction = f[0].data*0. + 1

            for i, ii in enumerate(range(all_data["pixelranges"][imind][0],
                all_data["pixelranges"][imind][1])):
                for j, jj in enumerate(range(all_data["pixelranges"][imind][2],
                    all_data["pixelranges"][imind][3])):
                    if subtractions[imind,i,j] != 0:
                        f[0].data[ii, jj] = subtractions[imind,i,j]
                        pixels_not_modified_by_subtraction[ii,jj] = 0

            sky_inds = where((pixels_not_modified_by_subtraction == 1)*\
                (1 - isnan(f[0].data))*(1 - isinf(f[0].data)))
            f[0].data -= median(f[0].data[sky_inds])*\
                pixels_not_modified_by_subtraction

            f.writeto(newim, overwrite = True)
            f.close()

            f_ilist.write(newim + '\n')
            f_slist.write(origim.replace("cbcd.fits", "cbunc.fits") + '\n')
            f_mlist.write(origim.replace("cbcd.fits", "bimsk.fits") + '\n')

        f_ilist.close()
        f_slist.close()
        f_mlist.close()

        print("This needs to run with csh")

        if int(ee)==0:
            # Should be first so first_epoch_dir will be populated
            not_first_epoch = 0
            first_epoch_dir = wd
        else:
            not_first_epoch = 1
            f_mopex.write(f"\ncp {first_epoch_dir}/mosaic_fif.tbl {wd} \n")
            create_nofid_params(mopex, 'overlap', channel)
            create_nofid_params(mopex, 'mosaic', channel)

        ch = channel.replace('ch','I')

        f_mopex.write(f"cp -r {mopex}/cal {wd} \n")
        f_mopex.write(f"cp -r {mopex}/cdf {wd} \n")
        f_mopex.write(f"cp {mopex}/mopex-script-env.csh {wd} \n")

        f_mopex.write(f"cd {wd} \n")
        f_mopex.write("source mopex-script-env.csh \n")
        f_mopex.write(create_mopex_cmd('overlap', ee, channel)+' \n')
        f_mopex.write(create_mopex_cmd('mosaic', ee, channel)+' \n')

        cmd=f"mv -v Combine/mosaic.fits Combine/mosaic_{objname}_sub.fits"
        f_mopex.write(cmd+' \n')

        for subdir_to_remove in ["BoxOutlier", "ReInterp", "Overlap_Corr",
            "DualOutlier", "Outlier", "Interp"]:
            f_mopex.write(f"rm -fr {wd}/{subdir_to_remove} \n")

        try:
            f_mopex.close()
        except:
            print("Couldn't close file! I guess there's not one open.")

        # add to list of output run files
        print(f'Need to run {filename}')
        out_runfiles.append(filename)

    return(out_runfiles)

if __name__ == '__main__':
    # Testing the method on local machine
    insert_into_subtractions('/data/ckilpatrick/spitzer/ch1',
        '/data/software/mopex', 'ch1', 'AT2017gfo',
        email='ckilpatrick@northwestern.edu')
