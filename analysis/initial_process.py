import os
import glob
import multiprocessing as mp
import subprocess
import shutil

def filter0(list_of_fls):
        return [item for item in list_of_fls if item.count("_0000_0000_") == 0]

def write_out_inlist(fls, outfile):
        with open(outfile, 'w') as f:
            f.write('\n'.join(fls)+'\n')

def run_dir(var):
        dr, mopex, channel = var
        wd = os.path.join(dr, "working_dir")

        if os.path.exists(wd):
            shutil.rmtree(wd)
        os.makedirs(wd)

        fls = filter0(glob.glob(os.path.join(dr, channel, 'bcd/*cbcd.fits')))
        fls.sort()

        ufls = [fl.replace("cbcd.fits", "cbunc.fits") for fl in fls]
        mfls = [fl.replace("cbcd.fits", "bimsk.fits") for fl in fls]

        write_out_inlist(fls, os.path.join(wd, 'images.list'))
        write_out_inlist(ufls, os.path.join(wd, 'sigma.list'))
        write_out_inlist(mfls, os.path.join(wd, 'mask.list'))

        os.system(f"cp -r {mopex}/cal {wd}")
        os.system(f"cp -r {mopex}/cdf {wd}")
        os.system(f"cp {mopex}/mopex-script-env.csh {wd}")

        with open(os.path.join(wd, 'run.sh'), 'w') as f:
            bandname = channel.replace('ch','I')
            f.write(f'cd {wd} \n')
            f.write('source mopex-script-env.csh \n')
            f.write(f'overlap.pl -n overlap_{bandname}.nl -I images.list -S sigma.list -d mask.list > log1.txt 2>&1 \n')
            f.write(f'mosaic.pl -n mosaic_{bandname}.nl -I images.list -S sigma.list -d mask.list > log2.txt 2>&1 \n')

        cmd = f"csh {wd}/run.sh"
        print(cmd)
        print(subprocess.call(cmd, shell = True))

        for subdr in "BoxOutlier Detect DualOutlier Interp Medfilter Outlier Overlap_Corr ReInterp cal".split(None):
            cmd = f"rm -fr {os.path.join(wd, subdr)}"
            print(cmd)
            os.system(cmd)

def initial_process(basedir, mopex, channel='ch1', nprocesses=8):

    pool = mp.Pool(processes=nprocesses)

    drs = glob.glob(os.path.join(basedir, "*:*"))
    var = [(dr, mopex, channel) for dr in drs]
    pool.map(run_dir, var)

