field = input("Enter name of field (e.g. GOODSS): ").strip().upper()

print('Computing 1-sigma mJy depth per pixel for ', field)
def compute_1sig_perpixel(sig5table, field):

    '''
    input: table of 5-sigma depths per filter (ascii format)
    output: save numpy array of 1-sigma mJy depth per pixel (saved as background_noise_field.npy)
    '''

    # import table with 5-sigma depths per filter
    sig5_tab = Table.read(sig5table, format='ascii')

    # compute mJy fluxes for 5-sigma depths (mJy is units of SBIPIX model training and forward-modeling)
    def mJy_flux(mags):
        'input magnitude 5-sigma depths as entire array, output mJy fluxes as entire array'
        return (3631/1e-6)*(10**(mags/(-2.5)))

    sig5_mJy = mJy_flux(sig5_tab['PSF_total'])

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
            mask = np.isin(sig5_tab['flux'], filts)
            if not mask.any():
                print("None of those filters were found — no change made.")
            else:
                removed = list(sig5_tab['flux'][mask])
                sig5_tab.remove_rows(np.where(mask)[0])
                print(f"Removed: {removed}")

        elif response in ('o', 'order'):
            new_order = input("Enter filters in desired order (comma-separated, lowercase format): ").strip().split(',')
            new_order = [f.strip() for f in new_order]
            current = list(sig5_tab['flux'])

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
    response_fname = input(f"\nSaving 1-sigma mJy depth per pixel as background_noise_{field}.npy). Change name? ").strip()
    if response_fname:
        np.save(f'{response_fname}', sig1_mJy_perpixel)
        fname = response_fname
    else:
        np.save(f'background_noise_{field}.npy', sig1_mJy_perpixel)
        fname = f'background_noise_{field}.npy'

    return fname

fname = compute_1sig_perpixel('table_5sigdepths.txt', field)    
print(f"1-sigma mJy depth per pixel saved as {fname}")

print('Constructing magnitude bins for ', field)
