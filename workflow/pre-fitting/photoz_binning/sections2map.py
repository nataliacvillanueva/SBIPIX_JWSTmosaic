from astropy.io import fits
import numpy as np
from tqdm import tqdm
import h5py

'''
This script uses the output HDF5 file from Lazy fits over a given section of the image mosaic.
In the output Lazy photo-z fit, there is a row for every pixel in that section of the mosaic 
that is also in the segmentation map. 
We want to read that output and assign the given photo-z and P(z) to the corresponding
pixel in a redshift map mosaic (precreated, and empty to fill in with this script per section).
'''

print('Starting to read empty precreated mosaics in...')
# let's open the precreated empty photo-z and P(z) mosaics, as well as information on the section splits
f = h5py.File('mosaic.h5', 'a')
sections = f['sections/section_columns'] # information on the sections the mosaic is split into
print('Sections shape:', sections.shape)
print('Sections first entry shape:', sections[0].shape)

segfull = fits.open('/Volumes/dataSSD/GOODSS/segmap_goodss_v0.891.fits')[0].data
print(segfull.shape[0], sections[0].shape[0])

print('Beginning to process sections...')
for i in tqdm(range(len(sections))):
    # print('--------------------------')
    # print('Section ', i, ' out of ', len(sections)-1)
    # print('--------------------------')
    split_number = i

    # let's open the relevant input and output Lazy files for this section
    try:
        lazyin = fits.open('/Volumes/dataSSD/lazymosaic_zcats/lazymosaic_sec' + str(split_number) + '_inlazy.fits')[1].data
        lazyout_Pz = fits.open('/Volumes/dataSSD/lazymosaic_zcats/lazymosaic_sec' + str(split_number) + '_outlazy.fits')[2].data

    except (IndexError, FileNotFoundError):
        print('Output or input Lazy file is empty/nonexistent for section ', split_number, ', skipping...')
        continue

    shapePz = (len(lazyout_Pz[1:]), len(lazyout_Pz['Pz'][0, :1001]))
    shapecoord = (len(lazyout_Pz[1:]), 2)
    with h5py.File('/Volumes/dataSSD/lazymosaic_zcats/lazymosaic_sec' + str(split_number) + '_Pz.hdf5', "w") as f:
        grp_pzmap = f.create_group('Pz')
        grp_coords = f.create_group('coordinates')

        pztab = grp_pzmap.create_dataset(
            "Pz_values",
            shape=shapePz,
            dtype=np.float32,
            chunks=True,
            compression="gzip",
            compression_opts=6,
            shuffle=True,
        )

        coordstab = grp_coords.create_dataset(
            "pixel_coordinates",
            shape=shapecoord,
            dtype=np.int32,
            chunks=True,
            compression="gzip",
            compression_opts=6,
        )

        coordstab[:, 0] = np.array(lazyin['x'], dtype=int)
        coordstab[:, 1] = np.array(lazyin['y'], dtype=int)

        pztab[:, :] = np.array(lazyout_Pz['Pz'][1:, :1001], dtype=np.float32)

        # # Stripe width: pick something safe (e.g., 5 columns)
        # stripe_width = 5

        # for x0 in tqdm(range(0, len(sections[i]), stripe_width)):
        #     #print('Processing stripe starting at x = ', x0)
        #     x1 = min(x0 + stripe_width, len(sections[i]))  # avoid going past the section

        #     # Indices of the pixels in this stripe
        #     mask = (x >= x0) & (x < x1)
        #     x_stripe = x[mask] - x0  # adjust to local stripe coordinates
        #     y_stripe = y[mask]
            
        #     z_stripe = z[mask]
        #     Pz_stripe = Pz[mask, :].T  # keep L × N_stripe

        #     if len(z_stripe) == 0:
        #         print('No pixels in stripe ' + str(x0) + ' to ' + str(x1) + ' skipping...')
        #         continue  # No pixels in this stripe

        #     else:
        #         # Build small blocks
        #         #zmap_block = np.full((seg.shape[0], x1 - x0), np.nan, dtype=np.float64)
        #         Pzmap_block = np.full((Pz.shape[1], segfull.shape[0], x1 - x0), np.nan, dtype=np.float32)
        #         #print(Pzmap_block.shape, Pz_stripe.shape)

        #         # Fill the block
        #         #zmap_block[y_stripe, x_stripe] = z_stripe
        #         Pzmap_block[:, y_stripe, x_stripe] = Pz_stripe

        #         # Write to HDF5
        #         #print('Writing stripe to main mosaic...')
        #         #zmap[:, sections[i][x0:x1]] = zmap_block
        #         pzmap[:, :, x0:x1] = Pzmap_block


print('All sections processed.')
    
