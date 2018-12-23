from astropy.io import fits
from numpy import *
import gzip
import cPickle as pickle
import sys
import commands
import tqdm

def do_it(cmd):
    print cmd
    print commands.getoutput(cmd)


f_mopex = open("run_stacks.sh", 'w')
f_mopex.write("source mopex-script-env.csh\n")

basedirs = sys.argv[1:]
basedirs.sort()

for basedir in tqdm.tqdm(basedirs):
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
        

        pixels_not_modified_by_subtraction = f[0].data*0. + 1

        for i, ii in enumerate(range(all_data["pixelranges"][imind][0], all_data["pixelranges"][imind][1])):
            for j, jj in enumerate(range(all_data["pixelranges"][imind][2], all_data["pixelranges"][imind][3])):
                if subtractions[imind,i,j] != 0:
                    f[0].data[ii, jj] = subtractions[imind,i,j]
                    pixels_not_modified_by_subtraction[ii,jj] = 0

        sky_inds = where((pixels_not_modified_by_subtraction == 1)*(1 - isnan(f[0].data))*(1 - isinf(f[0].data)))
        f[0].data -= median(f[0].data[sky_inds])*pixels_not_modified_by_subtraction
        
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


    if int(basedir.split(":")[-1]) == 1:
        not_first_epoch = 0
    else:
        not_first_epoch = 1
        f_mopex.write("\ncp " + pwd + "/" + basedir[:-1] + "1/sub_stack/mosaic_fif.tbl" " " + pwd + "/" + basedir + "/sub_stack\n")

    f_mopex.write("\ncp -r /home/drubin/mopex/cal " + sub_stack_dir_abs + '\n')
    f_mopex.write("cd " + sub_stack_dir_abs + '\n')
    f_mopex.write("""overlap.pl -n overlap_I1""" + "_nofid"*not_first_epoch + """.nl -I images.list -S sigma.list -d mask.list """ + "-F mosaic_fif.tbl"*not_first_epoch + """ > log1.txt
mosaic.pl -n mosaic_I1""" + "_nofid"*not_first_epoch + """.nl -I images.list -S sigma.list -d mask.list """ + "-F mosaic_fif.tbl"*not_first_epoch + """ > log2.txt\n""")
    f_mopex.write("mv -v Combine/mosaic.fits " + " Combine/mosaic_" + basedir.replace(":", "-")  + "_sub.fits" + '\n')

    f_mopex.write("rm -fr " + sub_stack_dir_abs + "/cal\n")
    
    


f_mopex.write('\necho "Done with mopex runs" | mailx -s "Done" drubin@stsci.edu\n')
f_mopex.close()
