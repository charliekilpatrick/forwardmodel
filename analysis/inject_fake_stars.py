import copy
from astropy.io import fits
from astropy import wcs
from photutils.psf import FittableImageModel
from astropy.convolution import discretize_model

def inject_stars(file, prf, ras, decs, mags, imgext='SCI'):

    ihdu = fits.open(file)
    phdu = fits.open(prf)

    # Create a new data frame from original image data
    data = copy.copy(ihdu[imgext].data)
    hdr = copy.copy(ihdu[imgext].header)
    shape = data.shape

    # Get epsf
    oversample = phdu[0].header['OVRSAMPL']
    epsf = FittableImageModel(phdu[0].data, oversampling=oversample)

    # Get wcs and zero point
    w = wcs.WCS(ihdu[0].header)
    zpt = ihdu[0].header['ZPTMAG']

    # Get x,y coordinates of injected source
    i=0
    for r,d,m in zip(ras, decs, mags):
        i += 1
        pix_xy = w.all_world2pix([[r, d]], 1)[0]

        # Convert magnitude into a total flux in counts
        flux = 10**(-0.4 * (m - zpt))

        setattr(epsf, 'x_mean', pix_xy[1])
        setattr(epsf, 'y_mean', pix_xy[0])
        setattr(epsf, 'x_0', pix_xy[1])
        setattr(epsf, 'y_0', pix_xy[0])
        setattr(epsf, 'flux', flux)

        data += discretize_model(epsf, (0, shape[1]), (0, shape[0]),
            mode='oversample', factor=oversample)

        idx = str(i).zfill(2)
        hdr[f'FS{idx}X']=pix_xy[0]
        hdr[f'FS{idx}Y']=pix_xy[1]
        hdr[f'FS{idx}RA']=r
        hdr[f'FS{idx}DEC']=d
        hdr[f'FS{idx}FLX']=flux
        hdr[f'FS{idx}MAG']=m

    ihdu[imgext].data = data
    ihdu[imgext].header = hdr
    ihdu.writeto(file, overwrite=True, output_verify='silentfix')




