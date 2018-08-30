from astropy.io import fits
from numpy import *
import gzip
import cPickle as pickle
import sys
import commands

def do_it(cmd):
    print cmd
    print commands.getoutput(cmd)


f_mopex = open("run_stacks.sh", 'w')
f_mopex.write("source mopex-script-env.csh\n")


for basedir in sys.argv[1:]:
    [all_data, parsed, settings, SNCmat, Cmat] = pickle.load(gzip.open(basedir + "/subtraction/fit_results.pickle", 'rb'))
    
    for key in settings:
        print "settings:", key
        
    for key in all_data:
        print "all_data:", key

    print "Figuring out which images to look at..."
        
    images_to_work_with = []
    for i in range(settings["n_img"]):
        if settings["images"][i].count(basedir):
            images_to_work_with.append(i)

    assert len(images_to_work_with) > 0
    print "Found these ", images_to_work_with, len(images_to_work_with)
    
    wd = basedir + "/sub_stack/"
    wd_im = wd + "ims/"
    
    do_it("rm -fr " + wd)
    do_it("mkdir -p " + wd_im)
    
    f = fits.open(basedir + "/subtraction/residuals.fits")
    subtractions = f[0].data
    print "subtractions ", subtractions.shape
    f.close()
    
    
    f_ilist = open(basedir + "/sub_stack/images.list", 'w')
    f_slist = open(basedir + "/sub_stack/sigma.list", 'w')
    f_mlist = open(basedir + "/sub_stack/mask.list", 'w')


    for imind in images_to_work_with:
        newim = settings["images"][imind].replace("cbcd_merged.fits", "cbcd.fits").replace("subtraction", "sub_stack/ims")
        assert newim != settings["images"][imind]
        
        origim = settings["images"][imind].replace("cbcd_merged.fits", "cbcd.fits").replace("subtraction", "ch1/bcd")
    
        
        f = fits.open(origim)
        
        
        for i, ii in enumerate(range(all_data["pixelranges"][imind][0], all_data["pixelranges"][imind][1])):
            for j, jj in enumerate(range(all_data["pixelranges"][imind][2], all_data["pixelranges"][imind][3])):
                if subtractions[imind,i,j] != 0:
                    f[0].data[ii, jj] = subtractions[imind,i,j]
                
        f.writeto(newim, overwrite = True)
        f.close()
        
        f_ilist.write(newim + '\n')
        f_slist.write(origim.replace("cbcd.fits", "cbunc.fits") + '\n')
        f_mlist.write(origim.replace("cbcd.fits", "bimsk.fits") + '\n')
    f_ilist.close()
    f_slist.close()
    f_mlist.close()
    
    print "This needs to run with csh"
    
    pwd = commands.getoutput("pwd")
    sub_stack_dir_abs = pwd + "/" + basedir + "/sub_stack"

    f_mopex.write("\ncp -r /home/drubin/mopex/cal " + sub_stack_dir_abs + '\n')
    f_mopex.write("cd " + sub_stack_dir_abs + '\n')
    f_mopex.write("""overlap.pl -n overlap_I1.nl -I images.list -S sigma.list -d mask.list > log1.txt
mosaic.pl -n mosaic_I1.nl -I images.list -S sigma.list -d mask.list > log2.txt\n""")
    f_mopex.write("mv -v Combine/mosaic.fits " + " Combine/mosaic_" + basedir.replace(":", "-")  + "_sub.fits" + '\n')

    f_mopex.write("rm -fr " + sub_stack_dir_abs + "/cal\n")
    

f_mopex.close()
