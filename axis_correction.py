#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 11:35:06 2024

@author: amartinez
"""

# %%
# Modify the size of the axes from the output of SWarp to make them equal
import numpy as np
from astropy.io import fits
import os

pruebas = '/Users/amartinez/Desktop/for_people/for_Herve/pruebas/'

field = 20
for chip in range(2, 3):
    folder = pruebas
    fits_files = [f for f in sorted(os.listdir(folder)) if f.endswith('.fits') and f.startswith('70_pointings_f%s_c%s' % (field, chip))]
    
    for nf, f_file in enumerate(fits_files):
        hdul = fits.open(pruebas + f_file)
        image_data = hdul[0].data
        ejes = [hdul[0].header['NAXIS1'], hdul[0].header['NAXIS2']]
        print(f'Working with {f_file}')
        print(f'Original axes: {ejes}')

        if ejes[0] != ejes[1]:  # Axes are different
            print('Axes are different, padding required...')
            diff = abs(ejes[0] - ejes[1])
            pad_p = diff // 2
            impar = diff % 2

            if ejes[0] < ejes[1]:  # NAXIS1 is smaller, pad horizontally
                new_naxis1 = ejes[1]
                hdul[0].header['NAXIS1'] = new_naxis1
                hdul[0].header['CRPIX1'] += pad_p  # Update reference pixel
                if impar:
                    hdul[0].header['NAXIS1'] += 1
                padded_data = np.pad(image_data, ((0, 0), (pad_p, pad_p + impar)), mode='constant', constant_values=0)
                
            elif ejes[1] < ejes[0]:  # NAXIS2 is smaller, pad vertically
                new_naxis2 = ejes[0]
                hdul[0].header['NAXIS2'] = new_naxis2
                hdul[0].header['CRPIX2'] += pad_p  # Update reference pixel
                if impar:
                    hdul[0].header['NAXIS2'] += 1
                padded_data = np.pad(image_data, ((pad_p, pad_p + impar), (0, 0)), mode='constant', constant_values=0)

            # Update data and save
            hdul[0].data = padded_data
            hdul.writeto(folder + f'PADDED_{f_file}', overwrite=True)
            print(f'New axes: {hdul[0].header["NAXIS1"]}, {hdul[0].header["NAXIS2"]}')
            print(30 * '_')
        
        else:
            print(f'Axes are already equal in {f_file}')
#