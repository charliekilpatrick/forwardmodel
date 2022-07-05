import glob
import numpy as np
import sys
import subprocess

print("python step5_forward_model.py NGC_3256  NGC_3256_cand1  120.0 60.0  1 2 3 4 0 0 0 0")

basedir = sys.argv[1]
basedir = basedir.split(":")[0] # In case NGC_3256:8 or something like that is given

cand_name = sys.argv[2]

RA = sys.argv[3]
Dec = sys.argv[4]

epochs = sys.argv[5:]

assert len(epochs) == 8

all_ims = []
all_epochs = []
pwd = subprocess.getoutput("pwd")

for epoch in range(1, 9):
    this_glob = np.sort(glob.glob(pwd + "/" + basedir + ":" + str(epoch) + "/subtraction/*merged.fits"))
    all_ims.extend(this_glob)
    all_epochs.extend([epochs[epoch - 1]]*len(this_glob))

f = open("orig_paramfile_snmod.txt", 'r')
lines = f.read()
f.close()

lines = lines.replace("QQQQQ", "15")
lines = lines.replace("RRRRR", RA)
lines = lines.replace("DDDDD", Dec)
lines = lines.replace("EEEEE", str([1.]*len(all_ims)).replace(" ", ""))
lines = lines.replace("PPPPP", str(all_epochs).replace(" ", "").replace("'", ""))
lines = lines.replace("SSSSS", "16")
lines = lines.replace("AAAAA", "1")
lines = lines.replace("IIIII", str(all_ims).replace(" ", ""))



lines = lines.replace("/home/drubin/spitzer_prf", pwd)

lines = lines.split('\n')

for i in range(len(lines)):
    if lines[i].count("iterative_centroid") == 1:
        lines[i] = "iterative_centroid    0"
    if lines[i].count("bad_pixel_list"):
        lines[i] = "bad_pixel_list   '" + pwd + "/blank_bad_pixel_list.txt'"
    if lines[i].count("pixel_area_map"):
        lines[i] = "pixel_area_map   '" + pwd + "/blank_pam.fits'"



lines = '\n'.join(lines)

subprocess.getoutput("rm -fr " + cand_name)
subprocess.getoutput("mkdir " + cand_name)

f = open(cand_name + "/paramfile.txt", 'w')
f.write(lines)
f.close()



f = open("mod.sh", 'w')

f.write("""#!/bin/bash
#SBATCH --job-name=model
#SBATCH --partition=shared
#SBATCH --time=00-23:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile
cd """ + subprocess.getoutput("pwd") + "/" + cand_name + """
python /home/drubin/new_psf_phot/new_phot.py paramfile.txt > log.txt
""")
