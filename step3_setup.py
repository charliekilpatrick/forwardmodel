import glob
from astropy.io import fits
import commands


drs = glob.glob("NGC_232:*")
drs.sort()

mjds = []

for dr in drs:
    print dr

    first_file = glob.glob(dr + "/ch1/bcd/*bcd*fits")[0]
    print first_file
    f = fits.open(first_file)
    mjds.append(f[0].header["MJD_OBS"])
    f.close()

print mjds
print drs



print "Ready to erase any old results? Control-c if not!"
raw_input("")

for dr in drs:
    commands.getoutput("rm -f " + dr + "/subtraction")
    commands.getoutput("mkdir " + dr + "/subtraction")


f_all = open("run.sh", 'w')

for i in range(len(drs)):

    refs = []
    for j in range(len(drs)):
        if j != i:
            if mjds[j] < mjds[i] - 20. or mjds[j] > mjds[i] + 330.:
                refs.append(j)

    print i, refs
    
    wd = "subtraction_" + drs[i].replace("../", "")
    commands.getoutput("mkdir " + wd)

    f = open("orig_paramfile.txt", 'r')
    lines = f.read()
    f.close()

    fls = []
    for ind in [i] + refs:
        fls += glob.glob(drs[ind] + "/*merged.fits")
    fls = ["../" + item for item in fls]

    lines = lines.replace("IIIII", str(fls))
    lines = lines.replace("EEEEE", "[20.]*90 + [1.]*%i" % (90*len(refs)))

    f = open(wd + "/paramfile.txt", 'w')
    f.write(lines)
    f.close()

    assert len(fls) == 90*(len(refs) + 1)
    
    f_all.write("cd " + commands.getoutput("pwd") + "/" + wd + "\n")
    f_all.write("python ../new_phot.py paramfile.txt > log.txt\n")
f_all.close()

