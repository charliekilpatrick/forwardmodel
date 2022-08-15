import os

analysis_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(analysis_dir, '../data')

data="""oversample      5       # Should match psf
renormpsf       0       # only use if psf is not normalized
psf_has_pix     1       # ePSF containing pixel?

fitSNoffset         0       # 1 => fit offset
patch               QQQQQ   # patch size-- must be odd
fixmodelfromrefs    0       # 1=> fix the model from the references
use_DQ_list         1
n_cpu               32
n_iter              4

sndRA_offset        RRRAOFFSET   # SNRA - GalRA
sndDec_offset       DDECOFFSET  # SNdec - Galdec

RA0                 RRRRR       # Can be single number, or list
Dec0                DDDDD       # Can be single number, or list

images              IIIII

sciext              1
errext              2
dqext               3
okaydqs             [0]
errscale            EEEEE

epochs              PPPPP                           # 0 = reference, 1 = first epoch, ...
psfs                ["{data_dir}/FFFFF_psf_x5_v2.fits"]  # psf for galaxy, SN epoch 1, ...

splineradius        SSSSS
splinepixelscale    0.00015          # default is 2.10e-5

apodize             AAAAA

pixel_area_map      "{data_dir}/blank_pam.fits"
bad_pixel_list      "{data_dir}/blank_bad_pixel_list.txt"

SN_centroid_prior_arcsec    10
iterative_centroid          1
""".format(data_dir=data_dir)
