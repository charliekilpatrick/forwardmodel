import glob
from astropy.io import fits
from numpy import *
import sys
import tqdm
import find_coords
import os
import shutil
import param_data

def do_it(cmd):
    print(cmd)
    print(os.system(cmd))

def filter0(fls):
    return [item for item in fls if item.count("_0000_0000_") == 0]

max_ims = 1000 # 40 # Limit this to limit RAM or CPU time

def setup_subtractions(basedir, channel, email='ckilpatrick@northwestern.edu'):
    coords = find_coords.find_coords(basedir)
    larger_patch = list(coords.keys())

    print("larger_patch {0}".format(" ".join(larger_patch)))

    run_file = os.path.join(basedir, "run.sh")
    f_all = open(run_file, 'w')

    for gal_name in tqdm.tqdm(larger_patch):

        drs = glob.glob(os.path.join(basedir, gal_name + ":*"))
        drs.sort()

        mjds = []

        for dr in drs:
            print(f'Directory: {dr}')

            expr = os.path.join(dr, f"{channel}/bcd/*cbcd*fits")

            print(expr)

            first_file = glob.glob(expr)[0]
            print(f'First file: {first_file}')
            f = fits.open(first_file)
            mjds.append(f[0].header["MJD_OBS"])
            f.close()

        print(mjds)
        print(drs)

        for dr in drs:
            subtraction_dir = os.path.join(dr, 'subtraction')
            if os.path.exists(subtraction_dir):
                shutil.rmtree(subtraction_dir)
            os.makedirs(subtraction_dir)

            expr = os.path.join(dr, f'{channel}/bcd/*cbcd*.fits')
            cbcds = filter0(glob.glob(expr))
            cbcds.sort()

            for cbcd in cbcds:
                flname = cbcd.split("/")[-1]

                f1 = fits.open(cbcd)
                f1.append(f1[0])
                f1[0].data = zeros([0,0], dtype=float32)
                f1[1].name = "SCI"

                f2 = fits.open(cbcd.replace("cbcd.fits", "cbunc.fits"))
                f1.append(f2[0])
                f1[2].name = "ERR"
                f2.close()

                f3name = os.path.join(dr, 'working_dir/Rmask', flname.replace("cbcd.fits", "cbcd_rmask.fits"))
                f3 = fits.open(f3name)
                f1.append(f3[0])
                f1[3].name = "DQ"
                f3.close()

                print(f1.info())

                f1merge = os.path.join(dr, 'subtraction', flname.replace("cbcd.fits", "cbcd_merged.fits"))
                f1.writeto(f1merge, overwrite = True)
                f1.close()



        for i in range(len(drs)):

            refs = []
            for j in range(len(drs)):
                if j != i:
                    if mjds[j] < mjds[i] - 20. or mjds[j] > mjds[i] + 330.:
                        refs.append(j)

            print(i, drs[i], refs)

            wd = os.path.join(drs[i], "subtraction")

            lines = param_data.data

            fls = []
            first_epoch_count = None
            for ind in [i] + refs:
                expr = os.path.join(basedir, drs[ind], 'subtraction/*merged.fits')
                glob_results = glob.glob(expr)[:max_ims]
                glob_results.sort()
                fls += glob_results
                if first_epoch_count == None:
                    first_epoch_count = len(glob_results)

            lines = lines.replace("IIIII", str(fls))
            lines = lines.replace("EEEEE", "[20.]*%i + [1.]*%i" % (first_epoch_count, len(fls) - first_epoch_count))
            lines = lines.replace("PPPPP", "[0]*%i" % len(fls))

            coord_key = os.path.split(drs[i])[1].split(':')[0]
            lines = lines.replace("RRRRR", str(   coords[coord_key][0]   ))
            lines = lines.replace("DDDDD", str(   coords[coord_key][1]   ))

            if larger_patch.count(gal_name):
                lines = lines.replace("QQQQQ", "29")
                lines = lines.replace("SSSSS", "31")
                lines = lines.replace("AAAAA", "1")
            else:
                lines = lines.replace("QQQQQ", "21")
                lines = lines.replace("SSSSS", "20")
                lines = lines.replace("AAAAA", "0")

            param_file = os.path.join(wd, "paramfile.txt")
            log_file = os.path.join(wd, "log.txt")
            f = open(param_file, 'w')
            f.write(lines)
            f.close()

            f_all.write(f"cd {wd} \n")
            f_all.write(f"new_phot.py {param_file} > {log_file} 2>&1 \n")

    f_all.write(f'echo "Done with forward models" | mailx -s "Done" {email}\n')
    f_all.close()

    return(run_file)

