from astropy.io import fits
import numpy as np
import h5py
from tqdm import tqdm
from astropy.table import Table
from numba import jit

field = input("Enter field name (e.g., GOODSS): ")

# open segmentation map and photo-z output files for given slice of given field 
segfull = fits.open(f'/Volumes/dataSSD/{field}/segmap_{field.lower()}.fits')[0].data
# am going to hardcode dimensions of Pz array, because it should be the same for all slices (1751) since it's defined by the params file
pzdim = 1751

# set shapes for photo-z and zmap datasets, given length of array in Pz
shapePz = (pzdim, segfull.shape[0], segfull.shape[1])
shapez = (segfull.shape[0], segfull.shape[1])

split_tot_sections = 100 #means we MUST divide GOODSS and GOODSN into 100 sections
# creates splits according to same sections made for Lazy results, so that we can fill these in properly and iteratively
seg = np.array_split(segfull, int(split_tot_sections), axis=1)
col_splits = np.array_split(np.arange(segfull.shape[1]), int(split_tot_sections))

# these sections define the starting column coordinate (in dimension [1] of shape) for the 100 sections
sections = []
for cols in tqdm(col_splits):
    sections.append(cols)

# open zmap (integrated), and F444W data for given field
print('Opening integrated per pixel zmap...')
zmapfull = fits.open(f'/Volumes/dataSSD/{field}/{field.lower()}_intzmap_v0.95.fits')[0].data
print('Opening F444W data and error maps...')
F444full = fits.open(f'/Volumes/dataSSD/{field}/data/{field.lower()}_nrc_f444w_sci.fits')[0].data
F444errfull = fits.open(f'/Volumes/dataSSD/{field}/data/{field.lower()}_nrc_f444w_err.fits')[0].data

# writing mosaic.h5 which we will fill in using other scripts!
print('Starting to write mosaic.h5...')
with h5py.File(f"mosaic_{field.lower()}.h5", "w") as f:
    grp_sec = f.create_group('sections')
    grp_intzmap = f.create_group('integrated_zmap')
    grp_SNRmap = f.create_group('SNRmap')
    grp_pzmap = f.create_group('Pzmap')
    grp_pixzmap = f.create_group('perpixel_zmap')

    # initializing datasets for Pz and per pixel zmap from Lazy results (using failure values -1)
    pzmap = grp_pzmap.create_dataset(
        "Pz_values",
        shape=shapePz,
        dtype=np.float32,
        chunks=(pzdim, 32, 32),
        compression="gzip",
        compression_opts=1,
        fillvalue=np.nan
    )

    pixzmap = grp_pixzmap.create_dataset(
        "pixzmap",
        shape=shapez,
        dtype=np.float32,
        chunks=(64, 64),
        compression="gzip",
        compression_opts=4,
        fillvalue=np.nan
    )
    # don't fill these with any data because we've already initialized empty nan arrays! 

    # filling in integrated zmap
    zmap = grp_intzmap.create_dataset(
        "integrated_z_values",
        shape=shapez,
        dtype=np.float32,
        chunks=(64, 64),
        compression="gzip",
        compression_opts=4,
        fillvalue=np.nan
    )

    zmap[:] = zmapfull

    # filling in SNR F444W map
    SNRmap = grp_SNRmap.create_dataset(
        "SNR_F444W_values",
        shape=shapez,
        dtype=np.float32,
        chunks=(64, 64),
        compression="gzip",
        compression_opts=4,
        fillvalue=np.nan
    )
    
    SNRmap[:] = F444full / F444errfull

    # filling in section columns
    grp_sec.create_dataset("section_columns", data=np.array(sections))

    # rest of sections to do with per pixel z and Pz will be filled in by other scripts, so we don't need to do anything else here

