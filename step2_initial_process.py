import commands
import glob
import multiprocessing as mp
import subprocess


def filter0(list_of_fls):
    return [item for item in list_of_fls if item.count("_0000_0000_") == 0]



def run_dir(dr):
    wd = dr + "/working_dir/"
    commands.getoutput("rm -fr " + wd)
    commands.getoutput("mkdir " + wd)

    fls = glob.glob(commands.getoutput("pwd") + "/" + dr + "/ch1/bcd/*cbcd.fits")
    fls.sort()
    fls = filter0(fls)

    f = open(wd + "/images.list", 'w')
    f.write('\n'.join(fls) + '\n')
    f.close()

    fls = [item.replace("cbcd.fits", "cbunc.fits") for item in fls]
    f = open(wd + "/sigma.list", 'w')
    f.write('\n'.join(fls) + '\n')
    f.close()

    fls = [item.replace("cbunc.fits", "bimsk.fits") for item in fls]
    f = open(wd + "/mask.list", 'w')
    f.write('\n'.join(fls) + '\n')
    f.close()


    commands.getoutput("cp -r /home/drubin/mopex/cal " + wd)

    f = open(wd + "/run.sh", 'w')
    f.write("cd " + pwd + "/" + wd  + """
source ../../mopex-script-env.csh
overlap.pl -n overlap_I1.nl -I images.list -S sigma.list -d mask.list > log1.txt
mosaic.pl -n mosaic_I1.nl -I images.list -S sigma.list -d mask.list > log2.txt
""")
    f.close()

    cmd = "cd " + pwd + "/" + wd + "; csh run.sh"
    print cmd
    print subprocess.call(cmd, shell = True)

    for subdr in "BoxOutlier Detect DualOutlier Interp Medfilter Outlier Overlap_Corr ReInterp cal".split(None):
        cmd = "rm -fr " + wd + "/" + subdr
        print cmd
        print commands.getoutput(cmd)


pwd = commands.getoutput("pwd")

pool = mp.Pool(processes = 16)

drs = glob.glob("*:*")
pool.map(run_dir, drs)

