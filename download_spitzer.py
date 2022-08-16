#!/usr/bin/env python
import requests,os,sys,zipfile,shutil
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import utils
from astropy import units as u
import numpy as np

uri = 'http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/DataService'

# Color strings for download messages
green = '\033[1;32;40m'
red = '\033[1;31;40m'
end = '\033[0;0m'

def is_number(num):
    try:
        num = float(num)
    except ValueError:
        return(False)
    return(True)

def parse_coord(ra, dec):
    if (not (is_number(ra) and is_number(dec)) and
        (':' not in ra and ':' not in dec)):
        error = 'ERROR: cannot interpret: {ra} {dec}'
        print(error.format(ra=ra, dec=dec))
        return(None)

    if (':' in ra and ':' in dec):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        error = 'ERROR: Cannot parse coordinates: {ra} {dec}'
        print(error.format(ra=ra,dec=dec))
        return(None)

def post_request(coord, size=0.25):

    ra = coord.ra.degree
    dec = coord.dec.degree

    params = {'RA': ra, 'DEC': dec, 'SIZE': size, 'VERB': 3,
        'DATASET':'ivo%3A%2F%2Firsa.ipac%2Fspitzer.level1'}

    req = requests.get(uri, params=params)
    table = ascii.read(req.text)
    return(table)

# Input a table row from post_request, download the corresponding data product
# and unpack to target directory
def get_file(row, outrootdir = './', download_zip=True, filetype='cbcd'):

    url = row['accessWithAnc1Url']
    if 'NONE' in url or not url:
        return(1)

    unpack_file = os.path.join(outrootdir,
        row['externalname'].replace('bcd.fits', filetype+'.fits'))

    if os.path.exists(unpack_file):
        return(0)

    # This speeds up process if downloading multiple SNe where you might have
    # overlapping files from one object to next.  Also if you have to redo
    # the download process.
    if download_zip and os.path.exists(unpack_file):
        # Move downloaded file to ref_file
        shutil.copyfile(unpack_file, ref_file)
        return(0)

    message = 'Downloading file: {url}'

    if not os.path.exists(outrootdir):
        os.makedirs(outrootdir)

    if download_zip:
        download_file = os.path.join(outrootdir, 'download.zip')
    else:
        download_file = ref_file

    message = f'Downloading file: {unpack_file}'
    sys.stdout.flush()
    dat = utils.data.download_file(url, cache=False,
        show_progress=False, timeout=120)
    shutil.copyfile(dat, download_file)
    #os.chmod(download_file, 0775)
    message = '\r' + message
    message += green+' [SUCCESS]'+end+'\n'
    sys.stdout.write(message.format(url=url))

    # Unpack the file if zip
    if download_zip:
        # Unzip file
        zip_ref = zipfile.ZipFile(download_file, 'r')
        zip_ref.extractall(outrootdir)
        zip_ref.close()

        os.remove(download_file)

    return(0)

def download_from_coord(coord, outdir='./', filetypes=['cbcd','cbunc','bimsk']):

    # Get the table of all Spitzer data
    table = post_request(coord)
    table.sort('scet')

    # Get only IRAC rows of table
    mask = np.array(['IRAC' in row['modedisplayname'] for row in table])

    table = table[mask]
    n = len(table)
    print(f'There are {n} records to download...')

    for row in table:
        for filetype in filetypes:
            get_file(row, outrootdir=outdir, filetype=filetype)

if __name__ == '__main__':

    coord = parse_coord(sys.argv[1], sys.argv[2])

    download_from_coord(coord, outdir=sys.argv[3])
