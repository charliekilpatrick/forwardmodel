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
import pickle as pickle
import time
import json

import warnings
warnings.filterwarnings('ignore')

version = 1.33

# version history:
# 1.0 05-01-2018: First release
# 1.01 05-01-2018: Update to handle patches off the edge.
# 1.02 05-02-2018: added SN_centroid_prior_arcsec setting
# 1.03 05-05-2018: added fitSNoffset
# 1.04 05-24-2018: won't crash if an image has no good pixels
# 1.05 05-25-2018: Checks for NaN's in science data (not just error)
# 1.06 05-27-2018: Added apodize setting
# 1.07 06-01-2018: interp2d has a bug. Switched to SmoothBivariateSpline
# 1.08 06-03-2018: Added summary of large pulls, possibly an indicator for bad pixels. Fixed DoF bug.
# 1.09 06-04-2018: RA0, Dec0 can now be image-dependent (for when images are not aligned)
# 1.10 07-12-2018: Better derivative scales on flux
# 1.11 07-13-2018: Fixed map_coords! prefilter = True now
# 1.12 07-13-2018: Fixed map_coords for the galaxy! prefilter = True now
# 1.13 07-16-2018: Identifies SN flux in iteration
# 1.2 07-16-2018: Checks for convergence
# 1.21 07-18-2018: Added option for iterative centroiding, useful for very large numbers of dithers
# 1.3 12-14-2018: Takes one PSF for each image now.
# 1.31 12-23-2018: Doesn't use sparse Jacobian for galaxy-only fit
# 1.32 1-21-2019: Saves model of just point source
# 1.32 01-02-2019: Fixed epochs bug when starting with non-zero epoch
# 1.33 01-16-2021: Dumps json of parsed to result file


print("version: ", version)



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

    print("Final settings:", settings)
    return settings

def robust_index(dat, i1, i2, j1, j2, fill_value = 0):
    sh = dat.shape

    paddat = zeros([sh[0] + (i2 - i1)*2, sh[1] + (j2 - j1)*2], dtype=dat.dtype) + fill_value
    paddat[i2 - i1: i2 - i1 + sh[0], j2 - j1: j2 - j1 + sh[1]] = dat

    return paddat[i1 + (i2 - i1):i2 + (i2 - i1),
                  j1 + (j2 - j1):j2 + (j2 - j1)]


def read_image(im, pam, bad_pix_list, exts, settings, RA0, Dec0):
    """Reads patches from image files."""

    #[badx, bady] = readcol(bad_pix_list, 'ii')
    badx = [] ; bady = []

    patch2 = int(floor(settings["patch"]/2.))

    data = []
    f = fits.open(im)

    try:
        mjd = 0.5*(f[0].header["EXPSTART"] + f[0].header["EXPEND"])
    except:

        try:
            mjd = f[0].header["BMJD_OBS"]
        except:
            print("Couldn't read EXPSTART/EXPEND/BMJD_OBS!")
            mjd = 0.


    w = wcs.WCS(f[exts[0]].header)
    pix_xy = w.all_world2pix([[RA0, Dec0]], 1)[0]

    #print("Found xy", im, pix_xy)
    pix_xy = array(around(pix_xy), dtype=int32)

    subxs = arange(settings["padsize"], dtype=float64)/settings["oversample"]
    subys = arange(settings["padsize"], dtype=float64)/settings["oversample"]

    subxs -= median(arange(settings["patch"]*settings["oversample"], dtype=float64)/settings["oversample"])
    subys -= median(arange(settings["patch"]*settings["oversample"], dtype=float64)/settings["oversample"])

    subxs += pix_xy[0]
    subys += pix_xy[1]

    #print("subxs, subys", subxs, subys)

    subxs, subys = meshgrid(subxs, subys)

    RAs, Decs =  w.all_pix2world(subxs, subys, 1)

    assert Decs.max() - Decs.min() > settings["splinepixelscale"]*(2*settings["splineradius"] + 1)*1.05, "Spline overfills patch!"

    #save_img(RAs, "RAs.fits")
    #save_img(Decs, "Decs.fits")

    subxs = arange(settings["patch"], dtype=float64)
    subys = arange(settings["patch"], dtype=float64)
    subxs -= median(subxs)
    subys -= median(subys)
    subxs += pix_xy[0]
    subys += pix_xy[1]
    subxs, subys = meshgrid(subxs, subys)
    pixel_sampled_RAs, pixel_sampled_Decs =  w.all_pix2world(subxs, subys, 1)

    pixel_sampled_js, pixel_sampled_is = meshgrid(arange(settings["patch"], dtype=float64), arange(settings["patch"], dtype=float64))

    #RADec_to_i = interp2d(x = reshape(pixel_sampled_RAs, settings["patch"]**2),
    #                      y = reshape(pixel_sampled_Decs, settings["patch"]**2),
    #                      z = reshape(pixel_sampled_is, settings["patch"]**2), kind = 'linear')
    #RADec_to_j = interp2d(x = reshape(pixel_sampled_RAs, settings["patch"]**2),
    #                      y = reshape(pixel_sampled_Decs, settings["patch"]**2),
    #                      z = reshape(pixel_sampled_js, settings["patch"]**2), kind = 'linear')


    RADec_to_i = SmoothBivariateSpline(x = reshape(pixel_sampled_RAs, settings["patch"]**2),
                                       y = reshape(pixel_sampled_Decs, settings["patch"]**2),
                                       z = reshape(pixel_sampled_is, settings["patch"]**2), kx = 1, ky = 1)
    RADec_to_j = SmoothBivariateSpline(x = reshape(pixel_sampled_RAs, settings["patch"]**2),
                                       y = reshape(pixel_sampled_Decs, settings["patch"]**2),
                                       z = reshape(pixel_sampled_js, settings["patch"]**2), kx = 1, ky = 1)

    tmp_bad_pix = zeros(f[exts[0]].data.shape, dtype=int32)
    for k in range(len(badx)):
        tmp_bad_pix[bady[k] - 1, badx[k] - 1] = 1

    pixelrange = [pix_xy[1] - patch2 - 1, pix_xy[1] + patch2, pix_xy[0] - patch2 - 1, pix_xy[0] + patch2]

    tmp_bad_pix = robust_index(tmp_bad_pix,
                               pixelrange[0], pixelrange[1], pixelrange[2], pixelrange[3])

    print(im,pixelrange[0], pixelrange[1], pixelrange[2], pixelrange[3])

    for ext in exts:
        data.append(robust_index(array(f[ext].data, dtype=float64),
                                 pixelrange[0], pixelrange[1], pixelrange[2], pixelrange[3]))


        for i in range(settings["patch"]):
            for j in range(settings["patch"]):
                radius = sqrt((i - patch2)**2 + (j - patch2)**2)
                if radius**2. >= (patch2 + 0.5)**2:
                    data[-1][i,j] = 0
                else:
                    if ext == exts[1]: #ERR
                        if settings["apodize"]:
                            data[-1][i,j] /= 1 - (radius/(patch2 + 0.5))**8.


    f.close()


    f = fits.open(pam)
    pixel_area_map = robust_index(array(f[1].data, dtype=float64),
                                  pixelrange[0], pixelrange[1], pixelrange[2], pixelrange[3],
                                  fill_value = 1.)
    f.close()

    return data, RAs, Decs, RADec_to_i, RADec_to_j, pixel_area_map, tmp_bad_pix, mjd, pixel_sampled_RAs, pixel_sampled_Decs, pixelrange


def get_data(settings, all_data):
    """Read in the data."""

    all_data["RADec_to_i"] = []
    all_data["RADec_to_j"] = []
    all_data["pixel_area_map"] = []
    all_data["mjd"] = []
    all_data["pixel_sampled_RAs"] = []
    all_data["pixel_sampled_Decs"] = []
    all_data["pixelranges"] = []

    for i in range(settings["n_img"]):
        print("Image ", i)
        data, RAs, Decs, RADec_to_i, RADec_to_j, pixel_area_map, tmp_bad_pix, mjd, pixel_sampled_RAs, pixel_sampled_Decs, pixelrange = read_image(im = settings["images"][i], pam = settings["pixel_area_map"][i], bad_pix_list = settings["bad_pixel_list"][i],
                                                                                                                                      exts = [settings[key][i] for key in ["sciext", "errext", "dqext"]], settings = settings, RA0 = settings["RA0"][i], Dec0 = settings["Dec0"][i])

        all_data["RADec_to_i"].append(RADec_to_i)
        all_data["RADec_to_j"].append(RADec_to_j)
        all_data["pixel_area_map"].append(pixel_area_map)
        all_data["mjd"].append(mjd)
        all_data["pixel_sampled_RAs"].append(pixel_sampled_RAs)
        all_data["pixel_sampled_Decs"].append(pixel_sampled_Decs)
        all_data["pixelranges"].append(pixelrange)

        DQ = array(data[2], int32)

        for okaydq in settings["okaydqs"]:
            DQ = bitwise_and(DQ, ~ okaydq) # E.g., 111 with okaydq = 2: ~ okaydq = 101, so DQ -> 101. 101 with okaydq = 2: ~ okaydq = 101, so DQ -> 101
        DQ += tmp_bad_pix

        invvars = (data[1]*settings["errscale"][i])**-2. * (DQ == 0)

        bad_inds = where(isinf(invvars) + isnan(invvars) + isnan(data[0]))
        invvars[bad_inds] = 0.
        data[0][bad_inds] = 0.
        assert all(invvars >= 0)

        all_data["invvars"].append(invvars)

        all_data["RAs"].append(RAs)
        all_data["Decs"].append(Decs)
        all_data["scidata"].append(data[0])

    all_data["mjd"] = array(all_data["mjd"])

    return all_data


def get_PSFs(settings):
    """Read in the PSFs. Store FFTs. In the future, don't bother storing FFTs for point sources!!!"""
    all_data = dict(scidata = [], invvars = [], RAs = [], Decs = [])

    all_data["psf_FFTs"] = {}
    all_data["psf_subpixelized"] = {}

    for psf in unique(settings["psfs"]):
        f = fits.open(psf)
        psfdata = array(f[0].data, dtype=float64)
        f.close()

        assert sum(psfdata == psfdata.max()) == 1, "PSF has multiple pixels at the same maximum!"

        print("psf shape ", psfdata.shape)
        model_shape = settings["patch"]*settings["oversample"]
        print("model shape ", model_shape)
        padsize = int(2**ceil(log2(max(max(psfdata.shape), model_shape))))
        print("padsize ", padsize)
        settings["padsize"] = padsize

        psfdata_pad = zeros([padsize]*2, dtype=float64)
        psfdata_pad[:psfdata.shape[0], :psfdata.shape[1]] = psfdata

        psfdata_odd = zeros([max(psfdata.shape) + (max(psfdata.shape) % 2 == 0)]*2, dtype=float64)
        padsize_odd = len(psfdata_odd)
        print("padsize_odd", padsize_odd)
        assert padsize_odd % 2 == 1

        psfdata_odd[:psfdata.shape[0], :psfdata.shape[1]] = psfdata

        if not settings["psf_has_pix"]:
            # Add the pixel convolution
            print("Adding pixel convolution!")

            pixel = zeros([padsize]*2, dtype=float64)
            pixel[:settings["oversample"], :settings["oversample"]] = 1.#/settings["oversample"]**2.
            #assert isclose(sum(pixel), 1.)

            psfdata_pad = array(real(ft.ifft2(ft.fft2(psfdata_pad) * ft.fft2(pixel))), dtype=float64)

            pixel = zeros([padsize_odd]*2, dtype=float64)
            pixel[:settings["oversample"], :settings["oversample"]] = 1.#/settings["oversample"]**2.
            psfdata_odd = array(real(ft.ifft2(ft.fft2(psfdata_odd) * ft.fft2(pixel))), dtype=float64)

        # Now, recenter
        maxinds = where(psfdata_pad == psfdata_pad.max())

        recenter = zeros([padsize]*2, dtype=float64)
        recenter[padsize - maxinds[0][0], padsize - maxinds[1][0]] = 1.

        psfdata_fft = ft.fft2(psfdata_pad) * ft.fft2(recenter)

        psf_test = ft.ifft2(psfdata_fft)
        assert abs(imag(psf_test)).max() < 1e-8

        psf_test = array(real(psf_test), dtype=float64)
        assert psf_test[0,0] == psf_test.max()

        all_data["psf_FFTs"][psf] = psfdata_fft

        maxinds = where(psfdata_odd == psfdata_odd.max())
        print("maxinds ", maxinds)
        recenter = zeros([padsize_odd]*2, dtype=float64)
        recenter[int(floor(psfdata.shape[0]/2.)) - maxinds[0][0], int(floor(psfdata.shape[1]/2.)) - maxinds[1][0]] = 1.

        psfdata = array(real(ft.ifft2(ft.fft2(psfdata_odd) * ft.fft2(recenter))), dtype=float64)

        all_data["psf_subpixelized"][psf] = psfdata

        save_img(all_data["psf_subpixelized"][psf], "psf_subpixelized.fits")

    return all_data

##########################################

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


def parseP(P, settings):
    """P is a vector of parameters."""

    parsed = {}

    ind = 0
    parsed["coeffs"] = reshape_coeffs(P[ind:ind+settings["n_coeff"]], settings["splineradius"])
    ind += settings["n_coeff"]

    parsed["dRA"] = P[ind:ind+settings["n_img"]]
    ind += settings["n_img"]
    parsed["dDec"] = P[ind:ind+settings["n_img"]]
    ind += settings["n_img"]

    parsed["sndRA_offset"] = P[ind]
    ind += 1
    parsed["sndDec_offset"] = P[ind]
    ind += 1

    parsed["SN_ampl"] = P[ind:ind+settings["n_epoch"]]
    ind += settings["n_epoch"]

    parsed["pt_RA"] = settings["RA0"] + parsed["dRA"] + parsed["sndRA_offset"]
    parsed["pt_Dec"] = settings["Dec0"] + parsed["dDec"] + parsed["sndDec_offset"]

    return parsed


def unparseP(parsed, settings):
    """parsed is a dictionary of parameters."""

    P = concatenate((unreshape_coeffs(parsed["coeffs"], settings["splineradius"]),
                     parsed["dRA"], parsed["dDec"], [parsed["sndRA_offset"]], [parsed["sndDec_offset"]], parsed["SN_ampl"]
                 ))
    return P

##########################################


#def indiv_model(args):
#    [RAs, Decs, settings, parsed] = args

def make_pixelized_PSF(parsed, i):

    icoord = all_data["RADec_to_i"][i](parsed["pt_RA"][i], parsed["pt_Dec"][i])[0,0]
    jcoord = all_data["RADec_to_j"][i](parsed["pt_RA"][i], parsed["pt_Dec"][i])[0,0]


    #print "i, icoord, jcoord", i, icoord, jcoord


    j2d, i2d = meshgrid(arange(settings["patch"], dtype=float64)*settings["oversample"],
                        arange(settings["patch"], dtype=float64)*settings["oversample"])

    psfsize = len(all_data["psf_subpixelized"][settings["psfs"][i]])

    i2d -= icoord*settings["oversample"] - floor(psfsize/2.)
    j2d -= jcoord*settings["oversample"] - floor(psfsize/2.)


    i1d = reshape(i2d, settings["patch"]**2)
    j1d = reshape(j2d, settings["patch"]**2)
    coords = array([i1d, j1d])

    pixelized_psf = map_coordinates(all_data["psf_subpixelized"][settings["psfs"][i]], coordinates = coords, order = 2, mode="constant", cval = 0, prefilter = True)
    pixelized_psf = reshape(pixelized_psf, [settings["patch"], settings["patch"]])

    return pixelized_psf


def indiv_model(args):
    [i, parsed, just_pt_flux] = args

    if just_pt_flux:
        pixelized_psf = make_pixelized_PSF(parsed, i)
        if settings["epochs"][i] > 0:
            convolved_model = (pixelized_psf/all_data["pixel_area_map"][i])*(parsed["SN_ampl"][settings["epochs"][i] - 1])
        else:
            convolved_model = pixelized_psf*0.
        return convolved_model



    # map_coordinates numbers starting from 0. E.g., radius = 3 => 0, 1, 2, (3), 4, 5, 6
    xs = (all_data["RAs"][i] - (settings["RA0"][i] + parsed["dRA"][i]))*cos(settings["Dec0"][i]/(180./pi))/settings["splinepixelscale"] + settings["splineradius"]
    ys = (all_data["Decs"][i] - (settings["Dec0"][i] + parsed["dDec"][i]))/settings["splinepixelscale"] + settings["splineradius"]

    #xs = (RAs - settings["RA0"])*cos(settings["Dec0"]/(180./pi))/settings["splinepixelscale"] + settings["splineradius"]
    #ys = (Decs - settings["Dec0"])/settings["splinepixelscale"] + settings["splineradius"]


    #save_img(xs, "RA_frame.fits")
    #save_img(ys, "Dec_frame.fits")

    xs1D = reshape(xs, settings["padsize"]**2)
    ys1D = reshape(ys, settings["padsize"]**2)

    coords = array([xs1D, ys1D])
    #print coords.shape

    #print parsed["coeffs"].shape
    subsampled_model = map_coordinates(parsed["coeffs"], coordinates = coords, order = 2, mode="constant", cval = 0, prefilter = True)
    subsampled_model = reshape(subsampled_model, [settings["padsize"], settings["padsize"]])

    #save_img(subsampled_model, "subsampled_model.fits")
    subsampled_convolved_model = ft.ifft2(ft.fft2(subsampled_model) * all_data["psf_FFTs"][settings["psfs"][i]])
    subsampled_convolved_model = array(real(subsampled_convolved_model), dtype=float64)

    #save_img(subsampled_convolved_model, "subsampled_convolved_model.fits")

    convolved_model = subsampled_convolved_model[settings["oversample2"]::settings["oversample"], settings["oversample2"]::settings["oversample"]]
    convolved_model = convolved_model[:settings["patch"], :settings["patch"]]# + parsed["sky"][i]


    #print "psf ", i, settings["epochs"][i], parsed["SN_ampl"], (parsed["SN_ampl"][settings["epochs"][i] - 1])

    if settings["epochs"][i] > 0:
        pixelized_psf = make_pixelized_PSF(parsed, i)
        convolved_model += (pixelized_psf/all_data["pixel_area_map"][i])*(parsed["SN_ampl"][settings["epochs"][i] - 1])


    if any(all_data["invvars"][i] != 0):
        sky_estimate = sum((all_data["scidata"][i] - convolved_model)*all_data["invvars"][i])/sum(all_data["invvars"][i])
    #print "i, sky_estimate ", i, sky_estimate
    else:
        sky_estimate = 0.

    convolved_model += sky_estimate

    return convolved_model

def modelfn(parsed, im_ind = None, just_pt_flux = 0):#, all_data, settings):
    """Construct the model."""

    #pindiv = partial(indiv_model, all_data = all_data, settings = settings, parsed = parsed)


    #args = [[all_data["RAs"][i], all_data["Decs"][i], settings, parsed] for i in range(settings["n_img"])]

    if im_ind == None:
        im_ind = list(range(settings["n_img"]))

    models = pool.map(indiv_model, [(i, parsed, just_pt_flux) for i in im_ind])

    return array(models)

def default_im_ind(im_ind):
    if im_ind == None:
        im_ind = list(range(settings["n_img"]))
    if type(im_ind) == int:
        im_ind = [im_ind]
    return im_ind

def pull_FN(parsed, im_ind = None):
    im_ind = default_im_ind(im_ind)

    models = modelfn(parsed, im_ind = im_ind)

    pulls = []
    for i, im in enumerate(im_ind):
        pulls.append((all_data["scidata"][im] - models[i])*sqrt(all_data["invvars"][im]))
        #assert 1 - any(isnan(pulls[-1])), str(pulls[-1]) + " " + str(len(pulls)) + " " + str(list(all_data["scidata"][im])) + " " + str(list(all_data["invvars"][im]))

    pulls = array(pulls)
    pulls = reshape(pulls, settings["patch"]**2 * len(im_ind))

    #print parsed["pt_RA"].shape, settings["RA0"]
    #fff
    dRA_arcsec = (parsed["pt_RA"] - settings["RA0"])*cos(settings["Dec0"]/57.2957795)*3600.
    dDec_arcsec = (parsed["pt_Dec"] - settings["Dec0"])*3600.

    pulls = concatenate((pulls, dRA_arcsec/settings["SN_centroid_prior_arcsec"], dDec_arcsec/settings["SN_centroid_prior_arcsec"]))
    return pulls


def pull_FN_wrapper(P, im_ind_wrap, makechi2 = 0):
    """For L-M"""
    im_ind = im_ind_wrap[0]
    assert type(im_ind) == list

    parsed = parseP(P, settings)
    pulls = pull_FN(parsed, im_ind)

    if makechi2:
        return dot(pulls, pulls)
    else:
        return pulls #dot(pulls, pulls)


"""
def lstsq_fit_for_spline(parsed):
    print "Making jacobian for spline..."

    jacob = zeros([settings["patch"]**2 * settings["n_img"], settings["n_coeff"] + settings["n_epoch"]], dtype=float64)

    parsed["coeffs"] *= 0
    parsed["SN_ampl"] *= 0
    models0 = modelfn(parsed)

    reshaped_data = reshape(array(all_data["scidata"]), settings["patch"]**2 * settings["n_img"])
    reshaped_invvars = reshape(array(all_data["invvars"]), settings["patch"]**2 * settings["n_img"])

    for i in range(settings["n_coeff"]):
        if i%10 == 0:
            print i, settings["n_coeff"]
        coeffs = zeros(settings["n_coeff"], dtype=float64)
        coeffs[i] = 1.

        parsed["coeffs"] = reshape_coeffs(coeffs, radius = settings["splineradius"])
        models = modelfn(parsed)

        jacob[:, i] = reshape(models - models0, settings["patch"]**2 * settings["n_img"])*sqrt(reshaped_invvars)

    parsed["coeffs"] *= 0

    for i in range(settings["n_epoch"]):
        parsed["SN_ampl"] = zeros(settings["n_epoch"], dtype=float64)
        parsed["SN_ampl"][i] = 1.

        models = modelfn(parsed)
        jacob[:, i + settings["n_coeff"]] = reshape(models - models0, settings["patch"]**2 * settings["n_img"])*sqrt(reshaped_invvars)

    #save_img(jacob, "jacob.fits")


    coeffs1D = linalg.lstsq(jacob, reshaped_data * sqrt(reshaped_invvars))[0]
    parsed["coeffs"] = reshape_coeffs(coeffs1D[:settings["n_coeff"]], settings["splineradius"])
    parsed["SN_ampl"] = coeffs1D[settings["n_coeff"]: ]

    return parsed
"""

def LM_fit_for_centroids(parsed):
    P = unparseP(parsed, settings)

    pulls = pull_FN(parsed)
    print("chi^2 check before centroid", dot(pulls, pulls))
    assert 1 - isnan(dot(pulls, pulls))


    miniscale_parsed = dict(coeffs = reshape_coeffs(ones(settings["n_coeff"]), radius = settings["splineradius"]),
                            SN_ampl = ones(settings["n_epoch"], dtype=float64)*settings["flux_scale"],
                            sndRA_offset = 0,
                            sndDec_offset = 0,
                            dRA = zeros(settings["n_img"], dtype=float64),
                            dDec = zeros(settings["n_img"], dtype=float64))
    miniscale = unparseP(miniscale_parsed, settings)

    print("Running galaxy+SN-only fit", time.asctime())
    P, F, NA = miniLM_new(ministart = P, miniscale = miniscale, residfn = pull_FN_wrapper, passdata = list(range(settings["n_img"])), verbose = False, maxiter = 1, use_dense_J = True)
    print("Done", time.asctime())

    print("LM chi^2", F)


    #P, F, NA = miniNM_new(ministart = P, miniscale = miniscale, residfn = lambda x, y: pull_FN_wrapper(x, y, makechi2 = 1),
    #                      passdata = range(settings["n_img"]), verbose = True)

    #print "NM chi^2", F

    if settings["iterative_centroid"]:
        assert settings["fitSNoffset"] == 0

        for i in range(settings["n_img"]):
            print("Centroiding ", i, "of", settings["n_img"])
            tmp_dpos = zeros(settings["n_img"], dtype=float64)
            tmp_dpos[i] = 0.1

            miniscale_parsed = dict(coeffs = reshape_coeffs(zeros(settings["n_coeff"]), radius = settings["splineradius"]),
                                    SN_ampl = zeros(settings["n_epoch"], dtype=float64),
                                    sndRA_offset = 0,
                                    sndDec_offset = 0,
                                    dRA = tmp_dpos,
                                    dDec = tmp_dpos)
            miniscale = unparseP(miniscale_parsed, settings)
            P, F, Cmat = miniLM_new(ministart = P, miniscale = miniscale, residfn = pull_FN_wrapper, passdata = [i], verbose = False, maxiter = 3)
    else:

        print("Running centroid-only fit", time.asctime())

        miniscale_parsed = dict(coeffs = reshape_coeffs(zeros(settings["n_coeff"]), radius = settings["splineradius"]),
                                SN_ampl = zeros(settings["n_epoch"], dtype=float64),
                                sndRA_offset = 1.e-1 * settings["fitSNoffset"],
                                sndDec_offset = 1.e-1 * settings["fitSNoffset"],
                                dRA = zeros(settings["n_img"], dtype=float64) + 1e-1,
                                dDec = zeros(settings["n_img"], dtype=float64) + 1e-1)
        miniscale = unparseP(miniscale_parsed, settings)
        P, F, NA = miniLM_new(ministart = P, miniscale = miniscale, residfn = pull_FN_wrapper, passdata = list(range(settings["n_img"])), verbose = False, maxiter = 3, use_dense_J = True)
        print("LM chi^2", F)

        print("Running everything fit", time.asctime())

        miniscale_parsed = dict(coeffs = reshape_coeffs(ones(settings["n_coeff"])*settings["flux_scale"], radius = settings["splineradius"]),
                                SN_ampl = ones(settings["n_epoch"], dtype=float64)*settings["flux_scale"],
                                sndRA_offset = 1.e-1 * settings["fitSNoffset"],
                                sndDec_offset = 1.e-1 * settings["fitSNoffset"],
                                dRA = zeros(settings["n_img"], dtype=float64) + 1.e-1,
                                dDec = zeros(settings["n_img"], dtype=float64) + 1.e-1)
        miniscale = unparseP(miniscale_parsed, settings)
        P, F, Cmat = miniLM_new(ministart = P, miniscale = miniscale, residfn = pull_FN_wrapper, passdata = list(range(settings["n_img"])), verbose = False, maxiter = 3, use_dense_J = True)

        try:
            Cmat[0,0]
        except:
            Cmat = zeros([sum(miniscale != 0)]*2, dtype=float64)

        print("LM chi^2", F)

    parsed = parseP(P, settings)
    print('parsed["sndRA_offset"], parsed["sndDec_offset"]', parsed["sndRA_offset"], parsed["sndDec_offset"])

    pulls = pull_FN(parsed)
    print("iteration ", itr, "dRA, dDec", parsed["dRA"], parsed["dDec"])
    print("chi^2 check after centroid", dot(pulls, pulls))

    return parsed, Cmat



def get_loglike(P, all_data):
    """Compute the log likelihood."""
    return loglike

def show_model(the_model, all_data, settings):
    """Make the plots and images."""
    pass

def initialize_fit(settings, all_data):
    """Come up with a reasonable place to start the fit/MCMC."""
    return P

def do_fit(P, all_data, settings):
    """Run an LBFGS? fit."""
    return P

def do_MCMC(P, all_data, settings):
    """Run HMC? sampling."""
    return samples

def write_results(P, settings):
    """Save results to a results.txt file."""
    pass

def time_to_stop(parsed, last_flux, chi2, last_chi2, settings, itr):
    if itr >= settings["n_iter"]:
        print("Reached maximum iterations!")
        return 1
    if len(parsed["SN_ampl"]) == 0:
        if last_chi2 - chi2 < 0.01:
            print("Chi^2 converged for galaxy-only run!")
            return 1
    else:
        new_SN_ampl_with_floor = abs(parsed["SN_ampl"]) + abs(last_flux.max())
        old_SN_ampl_with_floor = abs(last_flux) + abs(last_flux.max())

        ampl_diff = abs(new_SN_ampl_with_floor/old_SN_ampl_with_floor).max()
        print("ampl_diff ", ampl_diff)

        if abs(ampl_diff - 1.) < 0.0001:
            return 1
    return 0

if __name__ == "__main__":
    pass
settings = read_paramfile(sys.argv[1])
settings = finish_settings(settings)


all_data = get_PSFs(settings)
all_data = get_data(settings, all_data)


settings["flux_scale"] = scoreatpercentile(all_data["scidata"], 99)
print('settings["flux_scale"]', settings["flux_scale"])


pool = multiprocessing.Pool(processes = settings["n_cpu"])

parsed = {}
parsed["coeffs"] = reshape_coeffs(ones(settings["n_coeff"])*settings["flux_scale"], radius = settings["splineradius"])
parsed["dRA"] = zeros(settings["n_img"], dtype=float64)
parsed["dDec"] = zeros(settings["n_img"], dtype=float64)
parsed["sndRA_offset"] = settings["sndRA_offset"]
parsed["sndDec_offset"] = settings["sndDec_offset"]
parsed["SN_ampl"] = ones(settings["n_epoch"], dtype=float64)*settings["flux_scale"]
parsed = parseP(unparseP(parsed, settings), settings)

#models = modelfn(parsed)#parsed, all_data, settings)


save_img(all_data["invvars"], "invvars.fits")
save_img(all_data["scidata"], "scidata.fits")
save_img(all_data["pixel_area_map"], "pixel_area_map_cutout.fits")
#save_img(all_data["RAs"], "RAs.fits")
#save_img(all_data["Decs"], "Decs.fits")
save_img(all_data["pixel_sampled_RAs"], "pixel_sampled_RAs.fits")
save_img(all_data["pixel_sampled_Decs"],"pixel_sampled_Decs.fits")

"""
some_RAs = settings["RA0"] + (random.random(size = 100) - 0.5)/3000.
some_Decs = settings["Dec0"] + (random.random(size = 100) - 0.5)/3000.

from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt


some_is = array([all_data["RADec_to_i"][52](some_RAs[i], some_Decs[i]) for i in range(100)])
print some_is.shape
some_is = some_is[:,0,0]
print some_is.shape



plt.scatter(some_RAs, some_Decs, c = some_is)
plt.xlim(settings["RA0"] - 0.0002, settings["RA0"] + 0.0002)
plt.ylim(settings["Dec0"] - 0.0002, settings["Dec0"] + 0.0002)
plt.colorbar()
plt.savefig("tmp.pdf")
plt.close()
"""

last_flux = zeros(len(parsed["SN_ampl"]), dtype=float64) - 2.
last_chi2 = 1e101
itr = 0
chi2 = 1e100

while time_to_stop(parsed, last_flux, chi2, last_chi2, settings, itr) == 0:
    last_chi2 = chi2
    last_flux = parsed["SN_ampl"]

    #parsed = lstsq_fit_for_spline(parsed)
    parsed, Cmat = LM_fit_for_centroids(parsed)


    pulls = pull_FN(parsed)
    print("iteration ", itr, "dRA, dDec", parsed["dRA"], parsed["dDec"])
    chi2 = dot(pulls, pulls)
    print("chi^2 ", chi2)
    print("SN_ampl", parsed["SN_ampl"])

    itr += 1

assert itr > 0, "No iterations run!"

models = modelfn(parsed)#parsed, all_data, settings)
save_img(models, "models.fits")
save_img([(all_data["scidata"][i] - models[i])*(all_data["invvars"][i] > 0) for i in range(settings["n_img"])], "residuals.fits")
pulls = array([(all_data["scidata"][i] - models[i])*sqrt(all_data["invvars"][i]) for i in range(settings["n_img"])])
save_img(pulls, "pulls.fits")

pt_models = modelfn(parsed, just_pt_flux = 1)#parsed, all_data, settings)
save_img(pt_models, "pt_models.fits")

try:
    SNCmat = parseCmat(Cmat, settings)
    print('Successfully parsed SNCmat')
except:
    SNCmat = zeros([len(parsed["SN_ampl"])]*2)
    print('Populating SNCmat with zeros')

print('Dumping all_data, parsed, settings, SNCmat, Cmat into fit_results.pickle')
pickle.dump([all_data, parsed, settings, SNCmat, Cmat], gzip.open("fit_results.pickle", 'w'))

print("parsed ", parsed)
print("SNCmat ", SNCmat)

f = open("results.txt", 'w')

f.write("version  " + str(version) + '\n')
f.write("chi^2  " + str(chi2) + '\n')
n_pixels = (array(all_data["invvars"]) > 0).sum()
f.write("Npixels  " + str(n_pixels) + '\n')
f.write("DoF  " + str(n_pixels - len(Cmat)) + '\n')
for pull_thresh in [5, 10, 20, 50]:
    f.write("Npixels_with_pull_gt_" + str(pull_thresh) + "  " + str(sum(abs(pulls) > pull_thresh)) + '\n')


f.write('\n')
for i in range(settings["n_epoch"]):
    try:
        f.write("SN_A%s  %f  %f\n" % (settings["epoch_names"][i+1], parsed["SN_ampl"][i], sqrt(SNCmat[i,i])))
    except IndexError:
        print(f'ERROR: SNCmat does not have supernova data')
        continue

f.write('\n')



for i in range(settings["n_epoch"]):
    inds = where(settings["epochs"] == i+1)

    f.write("MJD_%s  %f\n" % (settings["epoch_names"][i+1], mean(array(all_data["mjd"])[inds])))

f.write('\n')
f.write('\nCmat:\n')

for i in range(settings["n_epoch"]):
    for j in range(settings["n_epoch"]):
        f.write(str(SNCmat[i,j]) + "  ")
    f.write('\n')

try:
    SNWmat = linalg.inv(SNCmat)
except:
    SNWmat = SNCmat*0

f.write('\nWmat:\n')
for i in range(settings["n_epoch"]):
    for j in range(settings["n_epoch"]):
        f.write(str(SNWmat[i,j]) + "  ")
    f.write('\n')

f.write("PARSED_JSON_BELOW\n")

for key in parsed:
    try:
        parsed[key] = parsed[key].tolist()
    except:
        pass

f.write(json.dumps(parsed) + '\n')

f.close()

print("Done!")

#P = initialize_fit(settings, all_data)
#P = do_fit(P, all_data, settings)
#write_results(P, settings)
#show_model(modelfn(parseP(P, settings), settings), all_data, settings)

# Explicit jacobian requirements:
# Spitzer: 2000 spline nodes * 1000 pixels * 1000 images * 8 bytes/entry = 16 GB
# WFC3 IR: 2000 spline nodes * 1000 pixels * 100 images * 8 bytes/entry = 1.6 GB

# Data volume for 5x oversampling: 1000 images * 256**2 * 8 bytes/entry * 2 RA/Dec = 1 GB

# Speed test: 9 ms per model per image at 5x. So 18 s per image to build jacobian explictly (5 hours for Spitzer). Clear approach: parallelize over images.
