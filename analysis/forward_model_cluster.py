#!/usr/bin/env python
from numpy import *
from astropy.io import fits
import multiprocessing
import sys
import os
from scipy.ndimage.interpolation import map_coordinates
from scipy.interpolate import interp2d, SmoothBivariateSpline
from scipy.stats import scoreatpercentile
from astropy import wcs
from scipy import fftpack as ft
# For some reason the code sometimes wants full analysis.DavidsNM and other
# times DavidsNM
try:
    from analysis.DavidsNM import save_img, miniLM_new, miniNM_new
    from analysis.forward_model import forward_model
except:
    from DavidsNM import save_img, miniLM_new, miniNM_new
    from forward_model import forward_model
import gzip
import pickle
import time
import json
import ray
ray.init()

import warnings
warnings.filterwarnings('ignore')

def message(msg):
    print('\n\n'+msg+'\n'+'#'*80+'\n'+'#'*80+'\n\n')

# version history:
# 1.00 05-01-2018: First release
# 1.01 05-01-2018: Update to handle patches off the edge.
# 1.02 05-02-2018: added SN_centroid_prior_arcsec setting
# 1.03 05-05-2018: added fitSNoffset
# 1.04 05-24-2018: won't crash if an image has no good pixels
# 1.05 05-25-2018: Checks for NaN's in science data (not just error)
# 1.06 05-27-2018: Added apodize setting
# 1.07 06-01-2018: interp2d has a bug. Switched to SmoothBivariateSpline
# 1.08 06-03-2018: Added summary of large pulls, possibly an indicator for bad
#                  pixels. Fixed DoF bug.
# 1.09 06-04-2018: RA0, Dec0 can now be image-dependent (for when images are not
#                  aligned)
# 1.10 07-12-2018: Better derivative scales on flux
# 1.11 07-13-2018: Fixed map_coords! prefilter = True now
# 1.12 07-13-2018: Fixed map_coords for the galaxy! prefilter = True now
# 1.13 07-16-2018: Identifies SN flux in iteration
# 1.20 07-16-2018: Checks for convergence
# 1.21 07-18-2018: Added option for iterative centroiding, useful for very large
#                  numbers of dithers
# 1.30 12-14-2018: Takes one PSF for each image now.
# 1.31 12-23-2018: Doesn't use sparse Jacobian for galaxy-only fit
# 1.32 01-21-2019: Saves model of just point source
# 1.32 01-02-2019: Fixed epochs bug when starting with non-zero epoch
# 1.33 01-16-2021: Dumps json of parsed to result file
# 1.34 07-24-2022: Refactored to run inside of pipeline and streamline code
# 1.35 08-17-2022: Created forward_model class to run within pipeline
version = 1.35

def parse_line(line):
    parsed = line.split("#")[0]
    parsed = parsed.split(None)
    if len(parsed) > 1:
        parsed = [parsed[0], eval(" ".join(parsed[1:]))]
        return parsed
    else:
        return None

def robust_index(dat, i1, i2, j1, j2, fill_value = 0):
    sh = dat.shape

    paddat = zeros([sh[0] + (i2 - i1)*2, sh[1] + (j2 - j1)*2],
        dtype=dat.dtype) + fill_value
    paddat[i2 - i1: i2 - i1 + sh[0], j2 - j1: j2 - j1 + sh[1]] = dat

    return paddat[i1 + (i2 - i1):i2 + (i2 - i1),
                  j1 + (j2 - j1):j2 + (j2 - j1)]

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
    parsed["coeffs"] = reshape_coeffs(P[ind:ind+settings["n_coeff"]],
        settings["splineradius"])
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

    parsed["pt_RA"] = settings["RA0"]+parsed["dRA"]+parsed["sndRA_offset"]
    parsed["pt_Dec"] = settings["Dec0"]+parsed["dDec"]+parsed["sndDec_offset"]

    return parsed

def unparseP(parsed, settings):
    """parsed is a dictionary of parameters."""

    P = concatenate((unreshape_coeffs(parsed["coeffs"],
        settings["splineradius"]), parsed["dRA"], parsed["dDec"],
        [parsed["sndRA_offset"]], [parsed["sndDec_offset"]], parsed["SN_ampl"]
    ))
    return P


####################################
@ray.remote
def indiv_model(self, args):
        [i, parsed, just_pt_flux] = args

        if just_pt_flux:
            pixelized_psf = self.make_pixelized_PSF(parsed, i)
            if self.settings["epochs"][i] > 0:
                pam=self.all_data["pixel_area_map"][i]
                SN_ampl=parsed["SN_ampl"][self.settings["epochs"][i] - 1]
                convolved_model = (pixelized_psf/pam)*(SN_ampl)
            else:
                convolved_model = pixelized_psf*0.
            return convolved_model

        # map_coordinates numbers starting from 0,
        # e.g., radius=3 => 0,1,2,(3),4,5,6
        dec_scale = cos(self.settings["Dec0"][i]/(180./pi))
        pscale = self.settings["splinepixelscale"]
        rad = self.settings["splineradius"]

        pra = parsed["dRA"][i] ; pdec = parsed["dDec"][i]

        ras = self.all_data["RAs"][i]-self.settings["RA0"][i]-pra
        des = self.all_data["Decs"][i]-self.settings["Dec0"][i]-pdec

        xs = ras * dec_scale / pscale + rad
        ys = des / pscale + rad

        xs1D = reshape(xs, self.settings["padsize"]**2)
        ys1D = reshape(ys, self.settings["padsize"]**2)

        coords = array([xs1D, ys1D])

        subsampled_model = map_coordinates(parsed["coeffs"],
                                           coordinates=coords,
                                           order=2,
                                           mode="constant",
                                           cval=0,
                                           prefilter=True)
        subsampled_model = reshape(subsampled_model,
            [self.settings["padsize"], self.settings["padsize"]])

        psf_fft = self.all_data["psf_FFTs"][self.settings["psfs"][i]]
        scm = ft.ifft2(ft.fft2(subsampled_model) * psf_fft)
        scm = array(real(scm), dtype=float64)

        o1=self.settings["oversample"]
        o2=self.settings["oversample2"]
        convolved_model = scm[o2::o1, o2::o1]
        convolved_model = convolved_model[:self.settings["patch"],
            :self.settings["patch"]]

        if self.settings["epochs"][i] > 0:
            pixelized_psf = self.make_pixelized_PSF(parsed, i)
            pam = self.all_data["pixel_area_map"][i]
            SN_ampl = parsed["SN_ampl"][self.settings["epochs"][i] - 1]
            convolved_model += (pixelized_psf/pam)*(SN_ampl)


        if any(self.all_data["invvars"][i] != 0):
            mod_resid = self.all_data["scidata"][i] - convolved_model
            invvar = self.all_data["invvars"][i]
            sky_estimate = sum(mod_resid*invvar)/sum(invvar)
        else:
            sky_estimate = 0.

        convolved_model += sky_estimate

        return convolved_model

##############################################################################

class forward_cluster_cluster(forward_model):

    def modelfn(self, parsed, im_ind = None, just_pt_flux = 0):
        """Construct the model."""

        if im_ind == None:
            im_ind = list(range(self.settings["n_img"]))

        futures = [indiv_model.remote(self, (i, parsed, just_pt_flux))
            for i in im_ind]
        models = ray.get(futures)

        return array(models)

    def pull_FN(self, parsed, im_ind = None):
        im_ind = self.default_im_ind(im_ind)

        models = self.modelfn(parsed, im_ind = im_ind)

        pulls = []
        for i, im in enumerate(im_ind):
            pulls.append((self.all_data["scidata"][im] - models[i])*\
                sqrt(self.all_data["invvars"][im]))

        pulls = array(pulls)
        pulls = reshape(pulls, self.settings["patch"]**2 * len(im_ind))

        dec_scale = cos(self.settings["Dec0"]/57.2957795)*3600.
        dRA_arcsec = (parsed["pt_RA"] - self.settings["RA0"])*dec_scale
        dDec_arcsec = (parsed["pt_Dec"] - self.settings["Dec0"])*3600.

        pulls = concatenate((pulls,
            dRA_arcsec/self.settings["SN_centroid_prior_arcsec"],
            dDec_arcsec/self.settings["SN_centroid_prior_arcsec"]))

        return pulls

    def do_main_reduction(self, parsed, settings):

        last_flux = zeros(len(parsed["SN_ampl"]), dtype=float64) - 2.
        last_chi2 = 1e101
        itr = 0
        chi2 = 1e100
        Cmat = None
        stop = 0
        max_iter = settings["n_iter"]
        n_img = settings["n_img"]

        while stop==0:

            last_chi2 = chi2
            last_flux = parsed["SN_ampl"]

            message(f"Running LM fit for iteration {itr+1} of {max_iter}")

            parsed, Cmat = self.LM_fit_for_centroids(parsed)

            pulls = self.pull_FN(parsed)
            pRA = parsed["dRA"]
            dDec = parsed["dDec"]
            chi2 = dot(pulls, pulls)
            chi2_fmt = '%11.2f'%chi2
            chi2_fmt = chi2_fmt.strip()

            print(f"chi^2 check for iter {itr+1} after centroid {chi2_fmt}")

            itr += 1

            stop = self.time_to_stop(parsed, last_flux, chi2, last_chi2,
                settings, itr)

        assert itr > 0, "No iterations run!"

        # Parse models, residuals, pulls, and pt_models and save
        models = self.modelfn(parsed)
        save_img(models, os.path.join(self.basedir, "models.fits"))

        residuals = [(self.all_data["scidata"][i] - models[i])*\
            (self.all_data["invvars"][i] > 0) for i in range(n_img)]
        save_img(residuals, os.path.join(self.basedir, "residuals.fits"))

        pulls = array([(self.all_data["scidata"][i] - models[i])*\
            sqrt(self.all_data["invvars"][i]) for i in range(n_img)])
        save_img(pulls, os.path.join(self.basedir, "pulls.fits"))

        pt_models = self.modelfn(parsed, just_pt_flux = 1)
        save_img(pt_models, os.path.join(self.basedir, "pt_models.fits"))

        try:
            SNCmat = parseCmat(Cmat, settings)
            print('Successfully parsed SNCmat')
        except:
            SNCmat = zeros([len(parsed["SN_ampl"])]*2)
            print('Populating SNCmat with zeros')

        self.create_fit_results(self.all_data, parsed,
            settings, SNCmat, Cmat, chi2, pulls)

if __name__ == "__main__":

    fm = forward_model_cluster(sys.argv[1])
    fm.do_main_reduction(fm.parsed, fm.settings)

    print("Done!")
