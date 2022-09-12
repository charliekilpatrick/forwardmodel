import os, glob
from astropy.io import fits
from astropy.table import Table

analysis_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.abspath(os.path.join(analysis_dir, '../data'))

def get_all_prf_data():
    prfs = glob.glob(os.path.join(data_dir, '*prf*.fits'))
    table = Table([['X'*100],[0],[0],['X'*100],['X'*10],['X'*20],['X'*10],[0],
        [0]],
        names=('file','channel','version','origfile','grid','size','type',
            'oversample','iter')).copy()[:0]

    for prf in prfs:

        hdu = fits.open(prf)
        if 'TYPE' not in hdu[0].header.keys():
            continue
        if hdu[0].header['TYPE'] in ['DX','DY']:
            continue
        size = str(hdu[0].header['NAXIS1'])+'x'+str(hdu[0].header['NAXIS2'])
        if 'GRID' in hdu[0].header.keys():
            grid = hdu[0].header['GRID']
        else:
            grid = ''
        table.add_row([prf, int(hdu[0].header['CHANNEL']),
            int(hdu[0].header['VERSION']),
            hdu[0].header['ORIGNAME'], grid, size, hdu[0].header['TYPE'],
            int(hdu[0].header['OVRSAMPL']), int(hdu[0].header['ITER'])])

    return(table)

def get_prf(ch, ver, oversample=5):

    prf_table = get_all_prf_data()
    mask = (prf_table['channel']==int(str(ch).lower().replace('ch',''))) &\
        (prf_table['oversample']==oversample)
    prf_table = prf_table[mask]


    if int(ver) in list(prf_table['version'].data):
        mask = prf_table['version']==int(ver)
        return(prf_table[mask][0]['file'])
    else:
        # Get latest version
        prf_table.sort('version')
        return(prf_table[-1]['file'])

def get_spatially_varying(ch, ver, oversample=5):

    prf = get_prf(ch, ver, oversample=oversample)
    hdu = fits.open(prf, mode='readonly')

    if 'SPATIAL' in hdu[0].header.keys() and hdu[0].header['SPATIAL']==1:
        return(1)
    else:
        return(0)


data="""oversample      OOOOO            # Should match psf
renormpsf       0               # only use if psf is not normalized
psf_has_pix     1               # ePSF containing pixel?

fitSNoffset         0           # 1 => fit offset
patch               QQQQQ       # patch size-- must be odd
fixmodelfromrefs    0           # 1=> fix the model from the references
use_DQ_list         1
n_cpu               NNNNN
n_iter              4

sndRA_offset        RRRAOFFSET  # SNRA - GalRA
sndDec_offset       DDECOFFSET  # SNdec - Galdec

RA0                 RRRRR       # Can be single number, or list
Dec0                DDDDD       # Can be single number, or list

images              IIIII

sciext              1
errext              2
dqext               3
okaydqs             [0]
errscale            EEEEE

epochs              PPPPP       # 0 = reference, 1 = first epoch, ...
psfs                ["VVVVV"]  # psf for galaxy, SN epoch 1, ...

splineradius        SSSSS
splinepixelscale    0.00015     # default is 2.10e-5

apodize             AAAAA

spatial             KKKKK

pixel_area_map      "{data_dir}/blank_pam.fits"
bad_pixel_list      "{data_dir}/blank_bad_pixel_list.txt"

base_dir            BBBBB

SN_centroid_prior_arcsec    10
iterative_centroid          1
""".format(data_dir=data_dir)
