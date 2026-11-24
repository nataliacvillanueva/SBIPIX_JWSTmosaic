import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
from astropy.table import Table, QTable
import glob
import time
from tqdm import tqdm
import argparse

# set path to data
ssdpath = '/Volumes/dataSSD/GOODSS/'

# grab information from command line arguments
parser = argparse.ArgumentParser(description="A script that accepts options.")
parser.add_argument('--output_fname_prefix', type=str, help='file name prefix (e.g. gsz14_region, include section number if split with _sec#)')
parser.add_argument('--section_number', type=int, help='section number if splitting mosaic (first is 0)', default=None)
parser.add_argument('--tot_sections', type=int, help='total number of sections if splitting mosaic', default=None)
args = parser.parse_args()

# check if want to deal with cutout region of mosaic, and specify the bounds
# sets variable with cutout (or not) segmentation map
if args.tot_sections is not None and args.section_number is not None:
    split_tot_sections = args.tot_sections
    split_number = args.section_number
    seg = np.array_split(fits.open(ssdpath+'segmap_goodss_v0.891.fits')[0].data, int(split_tot_sections), axis=1)[int(split_number)]
else:
    regcheck = input('Cutout region of segmentation map? (cutout or none): ')
    seg = fits.open(ssdpath+'segmap_goodss_v0.891.fits')[0].data
    if regcheck == 'cutout':
        xmin = input('Designate bounding box, xmin: ')
        xmax = input('xmax: ')
        ymin = input('ymin: ')  
        ymax = input('ymax: ')  

        seg = seg[int(ymin):int(ymax), int(xmin):int(xmax)]
        
        # plots segmentation map as a check
        plt.imshow(seg)
        plt.show()
        plt.close('all')
    elif regcheck == 'none':
        pass
    
# specify image and error filepaths (IMPORTANT that these are in the same order as the flux and error column names are specified above)
photfiles = [
'goodss_nrc_f090w_sci_matchf444w.fits',
'goodss_nrc_f115w_sci_matchf444w.fits',
'goodss_nrc_f150w_sci_matchf444w.fits',
'goodss_nrc_f162m_sci_matchf444w.fits',
'goodss_nrc_f182m_sci_matchf444w.fits',
'goodss_nrc_f200w_sci_matchf444w.fits',
'goodss_nrc_f210m_sci_matchf444w.fits',
'goodss_nrc_f250m_sci_matchf444w.fits',
'goodss_nrc_f277w_sci_matchf444w.fits',
'goodss_nrc_f300m_sci_matchf444w.fits',
'goodss_nrc_f335m_sci_matchf444w.fits',
'goodss_nrc_f356w_sci_matchf444w.fits',
'goodss_nrc_f410m_sci_matchf444w.fits',
'goodss_nrc_f430m_sci_matchf444w.fits',
'goodss_nrc_f444w_sci.fits',
'goodss_nrc_f460m_sci.fits',
'goodss_nrc_f480m_sci.fits',
'goodss_acs_f435w_sci_matchf444w.fits',
'goodss_acs_f606w_sci_matchf444w.fits',
'goodss_acs_f775w_sci_matchf444w.fits',
'goodss_acs_f814w_sci_matchf444w.fits',
'goodss_acs_f850l_sci_matchf444w.fits'
]

errfiles = [
'goodss_nrc_f090w_err.fits',
'goodss_nrc_f115w_err.fits',
'goodss_nrc_f150w_err.fits',
'goodss_nrc_f162m_err.fits',
'goodss_nrc_f182m_err.fits',
'goodss_nrc_f200w_err.fits',
'goodss_nrc_f210m_err.fits',
'goodss_nrc_f250m_err.fits',
'goodss_nrc_f277w_err.fits',
'goodss_nrc_f300m_err.fits',
'goodss_nrc_f335m_err.fits',
'goodss_nrc_f356w_err.fits',
'goodss_nrc_f410m_err.fits',
'goodss_nrc_f430m_err.fits',
'goodss_nrc_f444w_err.fits',
'goodss_nrc_f460m_err.fits',
'goodss_nrc_f480m_err.fits',
'goodss_acs_f435w_err.fits',
'goodss_acs_f606w_err.fits',
'goodss_acs_f775w_err.fits',
'goodss_acs_f814w_err.fits',
'goodss_acs_f850l_err.fits',
]

# sets columns for Lazy input catalog, including coordinates to be able to make redshift maps of output fits later
flux_cols = ['f_f090w', 'f_f115w', 'f_f150w', 'f_f162m', 'f_f182m', 'f_f200w', 'f_f210m', 'f_f250m', 'f_f277w', 'f_f300m', 'f_f335m', 'f_f356w', 'f_f410m', 'f_f430m', 'f_f444w', 'f_f460m', 'f_f480m', 'f_f435w', 'f_f606w', 'f_f775w', 'f_f814w', 'f_f850lp']
err_cols = ['e_f090w', 'e_f115w', 'e_f150w', 'e_f162m', 'e_f182m', 'e_f200w', 'e_f210m', 'e_f250m', 'e_f277w', 'e_f300m', 'e_f335m', 'e_f356w', 'e_f410m', 'e_f430m', 'e_f444w', 'e_f460m', 'e_f480m', 'e_f435w', 'e_f606w', 'e_f775w', 'e_f814w', 'e_f850lp']
coords_cols = ['x', 'y']
names = (['ID'] + coords_cols + flux_cols + err_cols)

# import ZP and dust correlation values per pixel
# table is organized in same filter order as loop goes through
zp_dustcorr_tab = Table.read('goodss_zp_mymodifications.txt', format='ascii')
zp = zp_dustcorr_tab['ZP']
dustcorr = zp_dustcorr_tab['dustcorr']

# create empty array to store values, where number of rows corresponds to nonzero values in segmentation map region
emptydata = np.empty((np.count_nonzero(seg), len(names)))
print('Number of pixels to be included in catalog: ', np.count_nonzero(seg))

# create empty Lazy input catalog table, and add ID column for Lazy
cat = Table(emptydata, names=names)
cat['ID'] = np.arange(0, np.count_nonzero(seg))
print('Initialized catalog table')

# run loop to store values in input catalog table
for f in range(len(photfiles)):
    # print diagnostic of which filter is currently running
    print('Now running: ')
    print(photfiles[f])
    print(errfiles[f])
    if args.tot_sections is not None and args.section_number is not None:
        photcut = np.array_split(fits.open(ssdpath+photfiles[f])[0].data, int(split_tot_sections), axis=1)[int(split_number)]
        errcut = np.array_split(fits.open(ssdpath+errfiles[f])[0].data, int(split_tot_sections), axis=1)[int(split_number)]
        print(photcut.shape, errcut.shape)
    else:
        # open image and error data, making cutout if that was specified earlier
        phot = fits.open(ssdpath+photfiles[f])[0].data
        if regcheck == 'cutout':
            photcut = phot[int(ymin):int(ymax), int(xmin):int(xmax)]
        elif regcheck == 'none':
            photcut = phot
        err = fits.open(ssdpath+errfiles[f])[0].data
        if regcheck == 'cutout':
            errcut = err[int(ymin):int(ymax), int(xmin):int(xmax)]
        elif regcheck == 'none':
            errcut = err


    # set counter to access corresponding row of table
    counter = 0
    coord_vals = np.ndenumerate(seg)
    totitems = seg.shape[0] * seg.shape[1]

    for coord, value in tqdm(coord_vals, total=totitems):
        # if value of pixel in segmentation map is 0, then pass the pixel
        if value == 0:
            continue
        # if value is not zero, then save the flux and error values for the pixel
        elif value != 0:
            # save the x and y coordinates, where image has (y,x) coordinates
            cat[counter]['x'] = coord[1]
            cat[counter]['y'] = coord[0]

            # save flux and error values for pixel, where taking f-th element of flux and error column arrays (these have the same filter ordering so f-th corresponds to same filter)
            cat[counter][flux_cols[f]] = photcut[coord[0], coord[1]]*zp[f]*1.0e32*dustcorr[f]
            cat[counter][err_cols[f]] = errcut[coord[0], coord[1]]*zp[f]*1.0e32*dustcorr[f]

            # add to counter to move on to next row
            counter += 1

    # diagnostic check
    print(counter)

print('Cleaning up now...')
# fix bad values in catalog table
# this handles values where data is missing in that filter, or where there was positive flux but nan error (due to convolution to F444W scale) 
# had previously been doing nan values for all, but moved to 1e-4 flux and 1e12 error to be consistent with Steve's catalogs
for f in range(len(flux_cols)):
    # deals with issue from convolution to F444W
    nan_idxs = np.isnan(cat[err_cols[f]])
    cat[flux_cols[f]][nan_idxs] = -0.0001#np.nan
    cat[err_cols[f]][nan_idxs] = 1000000000000.0#np.nan
    # deals with missing data in filter
    zero_flux_idxs = np.where(cat[flux_cols[f]] == 0.0)
    cat[err_cols[f]][zero_flux_idxs] = 1000000000000.0#np.nan
    cat[flux_cols[f]][zero_flux_idxs] = -0.0001#np.nan

print('Saving...')
# set output file name (set this up to work with my automated shell script)
output_fname = './lazymosaic_catalogsections/' + args.output_fname_prefix + '_inlazy.fits'

# write out input catalog
cat.write(output_fname, format='fits', overwrite=True)
print('Saved as: ', output_fname)


