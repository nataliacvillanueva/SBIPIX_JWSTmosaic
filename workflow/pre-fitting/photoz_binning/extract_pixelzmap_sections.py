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

field = input("Enter field name (e.g., GOODSS): ")

print('Starting to read empty precreated mosaics in...')
# let's open the precreated empty per pixel photoz and P(z) mosaics, as well as information on the section splits
f = h5py.File(f'mosaic_{field.lower()}.h5', 'a')
sections = f['sections/section_columns'] # information on the sections the mosaic is split into
pixzmap = f['perpixel_zmap/pixzmap'] # Precreated empty photo-z mosaic
Pzmap = f['Pzmap/Pz_values'] # Precreated empty P(z) mosaic

# create a new group and dataset storing coordinates for reference
# initialize the dataset as empty nan for now
shapecoord = ()
f.create_dataset('coordinates/pixel_coordinates',
    shape=shapecoord,
    dtype=np.int32,
    chunks=True,
    compression="gzip",
    compression_opts=6,
)

seg = fits.open(f'/Volumes/dataSSD/{field}/segmap_{field.lower()}.fits')[0].data
print('Sections shape:', sections.shape)
print('Sections first entry shape:', sections[0].shape)
print('Segmentation map shape:', seg.shape[0])
print('pixzmap shape:', pixzmap.shape)
print('Pzmap shape:', Pzmap.shape)

print('Beginning to process sections...')
for i in range(len(sections)):
    print('--------------------------')
    print('Section ', i, ' out of ', len(sections)-1)
    print('--------------------------')
    split_number = i

    # let's open the relevant input and output Lazy files for this section
    try:
        lazyin = fits.open('/Volumes/dataSSD/lazymosaic_zcats/lazymosaic_sec' + str(split_number) + '_inlazy.fits')[1].data
        with h5py.File('/Volumes/dataSSD/lazymosaic_zcats/lazymosaic_sec' + str(split_number) + '_outlazy.h5', 'r') as f:
            zgrid = f['metadata/zgrid'][:]
            zbest = f['results/zbest'][:]
            pz = f['pz/pz'][:]

            zgrid_shape = zgrid.shape
            zbest_shape = zbest.shape
            pz_shape = pz.shape
        
    except (IndexError, FileNotFoundError):
        print('Output or input Lazy file is empty/nonexistent for section ', split_number, ', skipping...')
        continue

    x = np.array(lazyin['x'], dtype=int)
    y = np.array(lazyin['y'], dtype=int)



    
    z = lazyout_z['z_best']
    #Pz = np.array(lazyout_Pz['Pz'][1:], dtype=np.float32) # null first row

    # Stripe width: pick something safe (e.g., 5 columns)
    stripe_width = 5

    for x0 in tqdm(range(0, len(sections[i]), stripe_width)):
        #print('Processing stripe starting at x = ', x0)
        x1 = min(x0 + stripe_width, len(sections[i]))  # avoid going past the section

        # Indices of the pixels in this stripe
        mask = (x >= x0) & (x < x1)
        x_stripe = x[mask] - x0  # adjust to local stripe coordinates
        y_stripe = y[mask]
        
        z_stripe = z[mask]
        #Pz_stripe = Pz[mask, :].T  # keep L × N_stripe

        if len(z_stripe) == 0:
            print('No pixels in stripe ' + str(x0) + ' to ' + str(x1) + ' skipping...')
            continue  # No pixels in this stripe

        else:
            # Build small blocks
            zmap_block = np.full((seg.shape[0], x1 - x0), np.nan, dtype=np.float64)
            #Pzmap_block = np.full((Pz.shape[1], seg.shape[0], x1 - x0), np.nan, dtype=np.float32)

            # Fill the block
            zmap_block[y_stripe, x_stripe] = z_stripe
            #Pzmap_block[:, y_stripe, x_stripe] = Pz_stripe

            # Write to HDF5
            #print('Writing stripe to main mosaic...')
            zmap[:, sections[i][x0:x1]] = zmap_block
            #Pzmap[:, :, sections[i][x0:x1]] = Pzmap_block
    
#SET SECTIONS INTO MAIN PRECREATED MOSAIC
print('Writing out full updated mosaics...')
f.close()
