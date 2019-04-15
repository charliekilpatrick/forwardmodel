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


def open_file(mopex_count, run_stacks_list):
    flnme = "run_stacks_%i.sh" % mopex_count
    
    f_mopex = open(flnme, 'w')
    f_mopex.write("sleep " + str(random.randint(120)) + '\n')
    f_mopex.write("source mopex-script-env.csh\n")
    dir_count = 0
    run_stacks_list.append(flnme)

    return f_mopex, dir_count, run_stacks_list

    

def close_file(f_mopex, mopex_count, num_ims, run_stacks_num_ims):
    f_mopex.close()
    run_stacks_num_ims.append(num_ims)
    
    dir_count = 0
    num_ims = 0

    return dir_count, run_stacks_num_ims



print commands.getoutput("rm -fv run_stacks_*.sh")

run_stacks_list = []
run_stacks_num_ims = []

sne_per_file = int(sys.argv[1])
basedirs = sys.argv[2:]

basedirs.sort()

mopex_count = 0
num_ims = 0

f_mopex, dir_count, run_stacks_list = open_file(mopex_count, run_stacks_list)


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

    num_ims += settings["n_img"]
            
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

    for subdir_to_remove in ["BoxOutlier", "ReInterp", "Overlap_Corr", "DualOutlier", "Outlier", "Interp"]:
        f_mopex.write("rm -fr " + sub_stack_dir_abs + "/" + subdir_to_remove + "\n")

    
    if int(basedir.split(":")[-1]) == 8:
        dir_count += 1
        
        if dir_count >= sne_per_file:
            dir_count, run_stacks_num_ims = close_file(f_mopex, mopex_count, num_ims, run_stacks_num_ims)
            
            if basedir != basedirs[-1]:
                mopex_count += 1

                f_mopex, dir_count, run_stacks_list = open_file(mopex_count, run_stacks_list)
                


try:
    dir_count, run_stacks_num_ims = close_file(f_mopex, mopex_count, num_ims, run_stacks_num_ims)
except:
    print "Couldn't close file! I guess there's not one open."

run_stacks_num_ims = array(run_stacks_num_ims)

print "Final image count for each file:"
for i in range(len(run_stacks_list)):
    print run_stacks_list[i], run_stacks_num_ims[i]

ind = argmax(run_stacks_num_ims)

f_mopex = open(run_stacks_list[ind], 'a')
f_mopex.write('\necho "Done with mopex runs" | mailx -s "Done" drubin@stsci.edu\n')
f_mopex.close()

print "Run these:"
print "/bin/csh"
print "source " + " & ; source ".join(run_stacks_list) + " &"

