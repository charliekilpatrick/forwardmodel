import glob
from astropy.io import fits
import commands
from numpy import *
import sys
import tqdm

def do_it(cmd):
    print cmd
    print commands.getoutput(cmd)


def filt0(fls):
    return [item for item in fls if item.count("_0000_0000_") == 0]

max_ims = 1000 # 40 # Limit this to limit RAM or CPU time

coords = {"NGC_232": [10.690756, -23.561488],
          "NGC_3256": [156.96334, -43.903763],
          "NGC_2623": [129.60029, 25.754693]}

coords = {"NGC_3256": [156.96350816118547, -43.90380553110768],
          "IRAS_23128-5919": [348.94484627544676, -59.0545035495522],
          "NGC_034": [2.777168302672514, -12.107777807931035],
          "NGC_5331_TARG1": [208.067485826228, 2.1015660804773133],
          "NGC_5331_TARG2":  [208.0684141, 2.108639844],
          "Arp_220": [233.73873553513712, 23.503220967430188],
          "MCG+12-02-001": [13.516223942691987, 73.08512446722798],
          "UGC_2369_TARG1": [43.50722187340597, 14.970613677652484],
          "UGC_2369_TARG2": [43.50700477, 14.97645362],
          "MCG+00-29-023": [170.3010610688463, -2.98405254442051],
          "UGC_8782": [208.0742960648129, 31.44641396833412],
          "UGC_4881_TARG1": [138.98143240285972, 44.3328842403561],
          "UGC_4881_TARG2": [138.9778276, 44.33103462],
          "IC_4687": [273.4155250235469, -57.72523290081323],
          "NGC_1614": [68.49983495194611, -8.579331308801555],
          "NGC_6926": [308.2753840314846, -2.0274892570084013],
          "NGC_5256": [204.57208277632142, 48.2757863068072],
          "NGC_3110": [151.00883156228585, -6.47480543268672],
          "IC_2810": [171.43787365478804, 14.676581918097405],
          "IRAS_10565+2448": [164.82573540588095, 24.542994805493556],
          "NGC_6090": [242.92026125149007, 52.4576809439665],
          "UGC_8387": [200.1472907029957, 34.139731801202686],
          "UGC_8335_TARG1": [198.8958362173512, 62.12462478522953],
          "UGC_8335_TARG2": [198.8780454, 62.12922522],
          "NGC_5257": [204.9705853159915, 0.8402193643363522],
          "Mk_848": [229.5256334812985, 42.745903088632886],
          "UGC_5101": [143.96533672946484, 61.35324885146309],
          "NGC_7130": [327.0812647337232, -34.951384389086144],
          "IRAS_18293-3413": [278.1712938706666, -34.191065819336536],
          "Mk_273": [206.1754647105241, 55.88687275927587],
          "MCG+08-18-012": [144.12879892102796, 48.4695610749106],
          "IC_1623_TARG1": [16.947925945241185, -17.507300364289787],
          "IC_1623_TARG2": [16.94389719, -17.50632799],
          "ESO_507-70": [195.71823291198265, -23.921664883700547],
          "NGC_232": [10.690775420819906, -23.561516125017977],
          "Arp_148": [165.9747272411593, 40.85005103483688],
          "NGC_1572": [65.67848216858545, -40.60095444513549],
          "NGC_6240": [253.24555688764707, 2.401172073586949],
          "MCG-03-12-002_TARG1": [65.33336172130593, -18.815999657235167],
          "MCG-03-12-002_TARG2": [65.33310252, -18.81090372],
          "Arp_302": [224.2530504177939, 24.617537193398924],
          "IRAS_03359+1523W": [54.695995883450365, 15.54802903116914],
          "IRAS_17208-0014": [260.84161431982784, -0.2833919830450202],
          "NGC_2623": [129.60046625437928, 25.7546002753708],
          "IC_563": [146.58473358464468, 3.0457444142970593],
          "NGC_3690_TARG1": [172.12959959670823, 58.56115247872474],
          "NGC_3690_TARG2": [172.1403054, 58.5628389]}


f_all = open("run.sh", 'w')


for gal_name in tqdm.tqdm(sys.argv[1:]):
    
    drs = glob.glob(gal_name + ":*")
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
    
    
    
    
    #print "Ready to erase any old results? Control-c if not!"
    #raw_input("")
    
    for dr in drs:
        do_it("rm -fr " + dr + "/subtraction")
        do_it("mkdir " + dr + "/subtraction")
        
        cbcds = filt0(glob.glob(dr + "/ch1/bcd/*cbcd*fits"))
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
            
            # NGC_232:1/working_dir/Rmask/SPITZER_I1_47507968_0084_0000_2_cbcd_rmask.fits
            f3 = fits.open(dr + "/working_dir/Rmask/" + flname.replace("cbcd.fits", "cbcd_rmask.fits"))
            
            f1.append(f3[0])
            f1[3].name = "DQ"
            f3.close()
            
            print f1.info()
            
            f1.writeto(dr + "/subtraction/" + flname.replace("cbcd.fits", "cbcd_merged.fits"), overwrite = True)
            f1.close()

    

    for i in range(len(drs)):
        
        refs = []
        for j in range(len(drs)):
            if j != i:
                if mjds[j] < mjds[i] - 20. or mjds[j] > mjds[i] + 330.:
                    refs.append(j)

        print i, drs[i], refs
    
        wd = drs[i] + "/subtraction/"
        

        f = open("orig_paramfile.txt", 'r')
        lines = f.read()
        f.close()

        fls = []
        first_epoch_count = None
        for ind in [i] + refs:
            glob_results = glob.glob(commands.getoutput("pwd") + "/" + drs[ind] + "/subtraction/*merged.fits")[:max_ims]
            glob_results.sort()
            fls += glob_results
            if first_epoch_count == None:
                first_epoch_count = len(glob_results)



        lines = lines.replace("IIIII", str(fls))
        lines = lines.replace("EEEEE", "[20.]*%i + [1.]*%i" % (first_epoch_count, len(fls) - first_epoch_count))
        lines = lines.replace("PPPPP", "[0]*%i" % len(fls))
        lines = lines.replace("RRRRR", str(   coords[drs[i].split(":")[0]][0]   ))
        lines = lines.replace("DDDDD", str(   coords[drs[i].split(":")[0]][1]   ))
        
        f = open(wd + "/paramfile.txt", 'w')
        f.write(lines)
        f.close()
        
        #assert len(fls) == 90*(len(refs) + 1), str(len(fls)) + " != 90*" + str(len(refs) + 1)
        
        f_all.write("cd " + commands.getoutput("pwd") + "/" + wd + "\n")
        f_all.write("python /home/drubin/new_psf_phot/new_phot.py paramfile.txt > log.txt\n")

f_all.write('echo "Done with forward models" | mailx -s "Done" drubin@stsci.edu\n')
f_all.close()



