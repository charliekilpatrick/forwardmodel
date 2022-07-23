#!/usr/bin/env python
from numpy import *
from astropy.io import fits
import multiprocessing
import sys
from scipy.ndimage.interpolation import map_coordinates
from scipy.interpolate import interp2d, SmoothBivariateSpline
from scipy.stats import scoreatpercentile
from astropy import wcs
from scipy import fftpack as ft
from analysis.DavidsNM import save_img, miniLM_new, miniNM_new
import gzip
import pickle
import time
import json

def parse_line(line):
    parsed = line.split("#")[0]
    parsed = parsed.split(None)
    if len(parsed) > 1:
        parsed = [parsed[0], eval(" ".join(parsed[1:]))]
        return parsed
    else:
        return None

def read_paramfile(parfl):
    """Read in the settings."""

    f = open(parfl)
    lines = f.read().split('\n')
    f.close()

    settings = {}

    for line in lines:
        parsed = parse_line(line)
        if parsed != None:
            settings[parsed[0]] = parsed[1]

    return settings

def finish_settings(settings):
    """Finalize and check settings."""

    assert settings["patch"] % 2 == 1, "patch must be odd!"
    settings["n_img"] = len(settings["images"])

    settings["fitSNoffset"] # Making sure it's there

    settings["epoch_names"] = [0]

    epoch_numbers = []
    for ep in settings["epochs"]:
        if settings["epoch_names"].count(ep) == 0:
            settings["epoch_names"].append(ep)
        epoch_numbers.append(settings["epoch_names"].index(ep))
    settings["epochs"] = array(epoch_numbers)



    settings["n_epoch"] = max(settings["epochs"])
    settings["oversample2"] = int(floor(settings["oversample"]/2.))
    settings["n_coeff"] = reshape_coeffs(coeffs = zeros(100000), radius = settings["splineradius"], return_coeffs_not_size = 0)

    for img in settings["images"]:
        assert img.count("_pam_") == 0
        assert img.count("_pam.") == 0

    assert settings["renormpsf"] == 0

    try:
        settings["iterative_centroid"]
    except:
        settings["iterative_centroid"] = 0

    if settings["iterative_centroid"] and settings["fitSNoffset"]:
        assert 0, "Can't iterate a SN centroid. All images must be fit!"

    for key in ["sciext", "errext", "dqext", "errscale", "pixel_area_map", "bad_pixel_list", "RA0", "Dec0", "psfs"]:
        if type(settings[key]) != list:
            settings[key] = [settings[key]]*settings["n_img"]

    if len(settings["psfs"]) == 1:
        settings["psfs"] = [settings["psfs"][0] for i in range(settings["n_img"])]

    for key in ["RA0", "Dec0"]:
        settings[key] = array(settings[key])

    for key in ["sciext", "errext", "dqext", "errscale", "pixel_area_map", "bad_pixel_list", "images", "epochs", "psfs"]:
        assert len(settings[key]) == settings["n_img"], key + " has wrong length!"

    return settings

def reshape_coeffs(coeffs, radius, return_coeffs_not_size = 1):
    reshaped = zeros([2*radius + 1, 2*radius + 1], dtype=float64)


    ind = 0
    for i in range(2*radius + 1):
        for j in range(2*radius + 1):
           if (i - radius)**2. + (j - radius)**2. < radius**2.:
               reshaped[i,j] = coeffs[ind]
               ind += 1

    if return_coeffs_not_size:
        return reshaped
    else:
        return ind

def unreshape_coeffs(coeffs, radius):
    OneD_coeffs = []

    for i in range(2*radius + 1):
        for j in range(2*radius + 1):
           if (i - radius)**2. + (j - radius)**2. < radius**2.:
               OneD_coeffs.append(coeffs[i,j])
    return array(OneD_coeffs)

def parseCmat(Cmat, settings):

    ind = 0
    ind += settings["n_coeff"]
    ind += settings["n_img"]
    ind += settings["n_img"]
    if settings["fitSNoffset"]:
        ind += 1
        ind += 1

    SNCmat = Cmat[ind:ind+settings["n_epoch"], ind:ind+settings["n_epoch"]]
    ind += settings["n_epoch"]

    return SNCmat

settings = read_paramfile(sys.argv[1])
settings = finish_settings(settings)

print('Getting pickle data')
pkl_data = gzip.open('fit_results.pickle', 'rb')
[all_data, parsed, settings, SNCmat, Cmat] = pickle.load(pkl_data)

print('n_coeff',settings['n_coeff'])
print('n_img',settings['n_img'])
print('fitSNoffset',settings['fitSNoffset'])

SNCmat = parseCmat(Cmat, settings)

print('SNCmat',SNCmat)
print('Cmat',Cmat)

f = open("results.txt", 'w')

for pull_thresh in [5, 10, 20, 50]:
    f.write("Npixels_with_pull_gt_" + str(pull_thresh) + "  " + str('blah') + '\n')


f.write('\n')
for i in range(settings["n_epoch"]):
    f.write("SN_A%s  %f  %f\n" % (settings["epoch_names"][i+1], 7.0, 7.0))

f.write('\n')

print(settings['n_epoch'])

for i in range(settings["n_epoch"]):
    inds = where(settings["epochs"] == i+1)

    f.write("MJD_%s  %f\n" % (settings["epoch_names"][i+1], 7.0))

f.write('\n')
f.write('\nCmat:\n')

for i in range(settings["n_epoch"]):
    for j in range(settings["n_epoch"]):
        f.write(str('blah') + "  ")
    f.write('\n')
