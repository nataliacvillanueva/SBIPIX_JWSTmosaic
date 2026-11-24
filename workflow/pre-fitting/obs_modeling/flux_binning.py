import numpy as np
from astropy.table import Table
import sys
from astropy.io import fits
from tqdm import tqdm
import matplotlib.pyplot as plt


field = input("Enter name of field (e.g. GOODSS): ").strip().upper()

print('Computing 1-sigma mJy depth per pixel for ', field)
def compute_1sig_perpixel(sig5table, field):

    '''
    input: table of 5-sigma depths per filter (ascii format), field name (e.g. GOODSS)
    output: save numpy array of 1-sigma mJy depth per pixel (saved as background_noise_field.npy)
    '''

    # import table with 5-sigma depths per filter
    sig5_tab = Table.read(sig5table, format='ascii')

    # compute mJy fluxes for 5-sigma depths (mJy is units of SBIPIX model training and forward-modeling)
    def mJy_flux(mags):
        'input magnitude 5-sigma depths as entire array, output mJy fluxes as entire array'
        return (3631/1e-6)*(10**(mags/(-2.5)))

    PSF_total_tier_input = input("Enter tier level for PSF_total: ").strip()
    PSF_total_tier = 'PSF_total' + '(' + PSF_total_tier_input + ')'
    sig5_mJy = mJy_flux(sig5_tab[PSF_total_tier])

    # add to running table with 5-sigma depths
    sig5_tab['PSF_total_mJy'] = sig5_mJy

    # compute 1-sigma mJy depth using pixel aperture correction
    def mJy_flux_1sig(mJy_5sig):
        'input mJy 5-sigma depths as entire array, output mJy 1-sigma depths as entire array'
        sig1depth = mJy_5sig/5

        area_circ1 = ((0.1)**2)*np.pi # 1 pixel aperture area in arcsec^2
        area_pix = 0.03*0.03 # pixel area in arcsec^2
        npix_aperture = area_circ1/area_pix # number of pixels in 1 pixel aperture
        sig1depth_perpixel = sig1depth/np.sqrt(npix_aperture) # 3 pixel aperture correction

        return sig1depth_perpixel

    sig1_mJy_perpixel = mJy_flux_1sig(sig5_tab['PSF_total_mJy'])

    # add to running table with 5-sigma depths
    sig5_tab['perpixel_mJy_1sig'] = sig1_mJy_perpixel

    while True:
        print(sig5_tab)
        response = input("\n[y]es continue / [n]o abort / [r]emove filter / [o]rder filters: ").strip().lower()

        if response in ('y', 'yes'):
            break

        elif response in ('n', 'no'):
            print("Aborting.")
            sys.exit(1)

        elif response in ('r', 'remove'):
            filts = input("Enter filter name(s) to remove (comma-separated or one at a time, lowercase format): ").strip().split(',')
            filts = [f.strip() for f in filts]
            mask = np.isin(sig5_tab['Filter'], filts)
            if not mask.any():
                print("None of those filters were found — no change made.")
            else:
                removed = list(sig5_tab['Filter'][mask])
                sig5_tab.remove_rows(np.where(mask)[0])
                print(f"Removed: {removed}")

        elif response in ('o', 'order'):
            new_order = input("Enter filters in desired order (comma-separated, lowercase format): ").strip().split(',')
            new_order = [f.strip() for f in new_order]
            current = list(sig5_tab['Filter'])

            if sorted(new_order) != sorted(current):
                print("Filter list doesn't match current table filters — no change made.")
                print(f"Current filters: {current}")
            else:
                idx = [current.index(f) for f in new_order]
                sig5_tab = sig5_tab[idx]
                print("Reordered.")

        else:
            print("Please enter 'y', 'n', 'r', or 'o'.")

    print("\nProceeding with:")
    print(sig5_tab)

    # save the 1-sigma mJy depth per pixel as a numpy array for use in SBIPIX model training and forward-modeling
    response_fname = input(f"\nSaving 1-sigma mJy depth per pixel as background_noise_{field}_{len(sig5_tab)}filters.npy). Change name? ").strip()
    if response_fname=='y':
        newname = input("Enter new filename (include .npy extension): ").strip()
        np.save(f'{newname}', sig5_tab['perpixel_mJy_1sig'])
        fname = newname
    elif response_fname=='n':
        np.save(f'background_noise_{field}_{len(sig5_tab)}filters.npy', sig5_tab['perpixel_mJy_1sig'])
        fname = f'background_noise_{field}_{len(sig5_tab)}filters.npy'

    return fname, sig5_tab

fname, sig5_tab = compute_1sig_perpixel(f'/Volumes/dataSSD/{field}/Results/{field.lower()}_5sig-depths_v0.95_PSFcorr.txt', field)    
print(f"1-sigma mJy depth per pixel saved as {fname}")

#print('Constructing magnitude bins for ', field)
def compile_magnitudes(field, sig5_tab):

    '''
        input: field name (e.g. GOODSS), 5sigma table that has selected filters
        output: save numpy array of 1-sigma mJy depth per pixel (saved as background_noise_field.npy)
    '''

    # import zero point and dust correction table (provided per filter) to apply later to magnitudes
    with open(f'{field}_v095_zp-dustcorr.txt') as fh:
        lines = fh.readlines()
    colnames = lines[1].lstrip('#').split()  # ['filter', 'zp', 'dustcorr']

    zp_dustcorr_tab = Table.read(
        f'{field}_v095_zp-dustcorr.txt',
        format='ascii',
        data_start=0,
        names=colnames
    )

    zp = dict(zip([str(f).strip() for f in zp_dustcorr_tab['filter']], zp_dustcorr_tab['zp']))
    dustcorr = dict(zip([str(f).strip() for f in zp_dustcorr_tab['filter']], zp_dustcorr_tab['dustcorr']))

    def get_filenames(instrument, filtername, sci_or_err):
        'input instrument (nrc or acs) and filter name (lowercase) and science or error flag (sci or err), output filename'
        if instrument == 'nrc':
            if sci_or_err == 'sci':
                if filtername in ['f444w', 'f460m', 'f480m']:
                    filename = f'{field.lower()}_{instrument}_{filtername}_{sci_or_err}.fits'
                else:
                    filename = f'{field.lower()}_{instrument}_{filtername}_{sci_or_err}_matchf444w.fits'
            elif sci_or_err == 'err':
                filename = f'{field.lower()}_{instrument}_{filtername}_{sci_or_err}.fits'
        elif instrument == 'acs':
            if sci_or_err == 'sci':
                filename = f'{field.lower()}_{instrument}_{filtername}_{sci_or_err}_matchf444w.fits'
            elif sci_or_err == 'err':
                filename = f'{field.lower()}_{instrument}_{filtername}_{sci_or_err}.fits'
        return filename

    # identify all pixels in the field mosaic by using the segmentation map
    print('Opening segmentation map for ', field)
    path = f'/Volumes/dataSSD/{field}/'
    segmap = fits.open(path+f'segmap_{field.lower()}.fits')
    seg = segmap[0].data

    # select pixels that are part of a galaxy and meet the SNR>2.5 cut in F444W
    mask = (seg!=0)
    # get the indices of the pixels that meet the criteria, in flattened array
    mask_idx = np.where(mask.ravel())[0]

    # save length of segmap to have number of pixels handy in array construction
    npix = seg[mask].shape[0]

    # make arrays to store magnitudes and magnitude errors for each filter, for each pixel
    mags = {f: np.full(npix, np.nan) for f in sig5_tab['Filter']}
    mags_err = {f: np.full(npix, np.nan) for f in sig5_tab['Filter']}
    filters = sig5_tab['Filter']

    # function to apply zero point and dust correction to fluxes
    def flux_correction(flux_or_err, zp, dustcorr):
        'input nJy fluxes or errors, zero point, and dust correction, output corrected mJy fluxes or errors'
        corr = flux_or_err*zp*1.0e32*dustcorr
        mJy = corr * 1e-3  # Convert nJy to microJy
        return mJy

    # functions to convert mJy fluxes and errors to magnitudes
    def flux_to_mag(flux_mJy):
        'input mJy fluxes, output magnitudes'
        mag = -2.5 * np.log10(flux_mJy * 1e-6 / 3631)
        return mag

    def err_to_mag(err_mJy, flux_mJy):
        'input mJy flux errors, output magnitude errors'
        mag_err = (2.5 * err_mJy)/(flux_mJy * np.log(10))
        return mag_err

    # function to compile magnitudes and magnitude errors for each filter and pixel into preconstructed arrays
    def compile_magnitudes_filter(filter, instrument, validpix_mask):
        'input filter name and valid pixel mask, output magnitudes and magnitude errors for each pixel in the filter'
        # get the science and error mosaics for the filter
        sci_mosaic = fits.open(path+'data/'+get_filenames(instrument, filter, 'sci'))
        err_mosaic = fits.open(path+'data/'+get_filenames(instrument, filter, 'err'))

        # mask by the valid pixels determined earlier
        sci = sci_mosaic[0].data[validpix_mask]
        err = err_mosaic[0].data[validpix_mask]

        # impose additional restrictions on valid pixels: that flux is positive and that error and flux are finite
        # also impose cut on SNR>1 per pixel (remember the only other cut we imposed was SNR>2.5 in F444W specifically)
        # decided to cut this out ^
        xtra_mask  = (sci > 0) & np.isfinite(err) & np.isfinite(sci)# & (sci/err > 1)

        # apply zero point and dust correction to fluxes and errors
        flux_corr = flux_correction(sci[xtra_mask], zp[filter], dustcorr[filter])
        err_corr = flux_correction(err[xtra_mask], zp[filter], dustcorr[filter])

        # convert corrected fluxes and errors to magnitudes and magnitude errors
        mag = flux_to_mag(flux_corr)
        mag_err = err_to_mag(err_corr, flux_corr)

        mags[filter][xtra_mask] = mag
        mags_err[filter][xtra_mask] = mag_err

    for f in filters:
        if f in ['f435w', 'f606w', 'f775w', 'f814w', 'f850l']:
            instrument = 'acs'
        else:
            instrument = 'nrc'
        print(f"Compiling magnitudes for {f} ({instrument})")
        compile_magnitudes_filter(f, instrument, mask)

    # save the magnitudes and magnitude errors as numpy arrays for use in SBIPIX model training and forward-modeling
    np.save(f'mags_{field}_{len(sig5_tab)}filters.npy', mags)
    np.save(f'mags_err_{field}_{len(sig5_tab)}filters.npy', mags_err)
    print(f"Magnitudes saved as mags_{field}_{len(sig5_tab)}filters.npy")
    print(f"Magnitude errors saved as mags_err_{field}_{len(sig5_tab)}filters.npy")

    return mags, mags_err
#mags, mags_err = compile_magnitudes(field, sig5_tab)

# print('Plotting magnitude distributions for ', field)
def plot_magnitude_distributions(mags, mags_err, sig5_tab):
    '''
        input: magnitudes and magnitude errors for each filter and pixel, 5sigma table that has selected filters
        output: save plots of magnitude distributions for each filter
    '''

    filters = sig5_tab['Filter']

    for f in filters:
        mag = mags[f]
        mag_err = mags_err[f]

        # plot histogram of magnitudes
        plt.figure(figsize=(8,6))
        plt.hist(mag[~np.isnan(mag)], bins=100000)
        plt.xlabel('Magnitude')
        plt.ylabel('Number of Pixels')
        plt.title(f'Magnitude Distribution for {f}')
        plt.savefig(f'magnitude_distribution_{f}.png')
        plt.close()
# plot_magnitude_distributions(mags, mags_err, sig5_tab)