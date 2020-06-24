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

for epoch in range(1, 9):
    this_glob = np.sort(glob.glob(basedir + ":" + str(epoch) + "/subtraction/*merged.fits"))
    all_ims.extend(this_glob)
    all_epochs.extend([epochs[epoch - 1]]*len(this_glob))

f = open("orig_paramfile.txt", 'r')
lines = f.read()
f.close()

lines = lines.replace("QQQQQ", "15")
lines = lines.replace("RRRRR", RA)
lines = lines.replace("DDDDD", Dec)
lines = lines.replace("EEEEE", str([1].*len(all_ims)).replace(" ", ""))
lines = lines.replace("PPPPP", str(all_epochs).replace(" ", ""))
lines = lines.replace("SSSSS", "16")
lines = lines.replace("AAAAA", "1")

lines = lines.split('\n')

for i in range(len(lines)):
    if lines[i].count("iterative_centroid") == 1:
        lines[i] = "iterative_centroid    0"
lines = '\n'.join(lines)

subprocess.getoutput("rm -fr " + cand_name)
subprocess.getoutput("mkdir " + cand_name)

f = open(cand_name + "/paramfile.txt", 'w')
f.write(lines)
f.close()

