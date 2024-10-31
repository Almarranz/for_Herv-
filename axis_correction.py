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
import sys
# pruebas = '/Users/amartinez/Desktop/for_people/for_Herve/pruebas/'
pruebas = '/home/data/alvaro/gns_test/pruebas/'
folder = '/Users/amartinez/Desktop/for_people/for_Herve/gns2/F20/SWarp/outputs/chip1/'

field = 20
# %%
target_size = 2048
for chip in range(1, 2):
    # folder = '/home/data/alvaro/gns_test/F%s/SWarp/chip2%s/'%(field, chip)
    fits_files = [f for f in sorted(os.listdir(folder)) if f.endswith('.fits') and f.startswith('70_pointings_f%s_c%s' % (field, chip))]
    fits_files = [f for f in sorted(os.listdir(folder)) if f.endswith('.fits') and f.startswith('70_pointings_f%s_c%s' % (field, chip)) or f.endswith('weight.fits') and f.startswith('70_pointings_f%s_c%s' % (field, chip))]
    
    for nf, f_file in enumerate(fits_files):
        hdul = fits.open(folder + f_file)
        image_data = hdul[0].data
        ejes = [hdul[0].header['NAXIS1'], hdul[0].header['NAXIS2']]
        print(f'Working with {f_file}')
        print(f'Original axes: {ejes}')
        
        arget_size = 2048
        axis1, axis2 = ejes

        # Calculate crop/padding amounts for each axis
        crop_pad_axis1 = (axis1 - target_size) // 2
        crop_pad_axis2 = (axis2 - target_size) // 2

        # Process for axis 1 if not equal to target size
        if axis1 != target_size:
            if axis1 > target_size:  # Crop
                start_x = crop_pad_axis1
                end_x = start_x + target_size
                image_data = image_data[:, start_x:end_x]
                hdul[0].header['CRPIX1'] -= start_x
            else:  # Pad
                pad_x = abs(crop_pad_axis1)
                image_data = np.pad(image_data, ((0, 0), (pad_x, pad_x)), mode='constant', constant_values=0)
                hdul[0].header['CRPIX1'] += pad_x

        # Process for axis 2 if not equal to target size
        if axis2 != target_size:
            if axis2 > target_size:  # Crop
                start_y = crop_pad_axis2
                end_y = start_y + target_size
                image_data = image_data[start_y:end_y, :]
                hdul[0].header['CRPIX2'] -= start_y
            else:  # Pad
                pad_y = abs(crop_pad_axis2)
                image_data = np.pad(image_data, ((pad_y, pad_y), (0, 0)), mode='constant', constant_values=0)
                hdul[0].header['CRPIX2'] += pad_y

        # Update header axes and save to file
        hdul[0].header['NAXIS1'] = target_size
        hdul[0].header['NAXIS2'] = target_size
        hdul[0].data = image_data
        hdul.writeto(folder + f'PROCESSED_{f_file}', overwrite=True)
        
        print(f'New axes: {hdul[0].header["NAXIS1"]}, {hdul[0].header["NAXIS2"]}')
        print(30 * '_')
    
        
        
# =============================================================================
#         if ejes[0] != ejes[1]:  # Axes are different
#             print('Axes are different, padding required...')
#             diff = abs(ejes[0] - ejes[1])
#             pad_p = diff // 2
#             impar = diff % 2
# 
#             if ejes[0] < ejes[1]:  # NAXIS1 is smaller, pad horizontally
#                 new_naxis1 = ejes[1]
#                 hdul[0].header['NAXIS1'] = new_naxis1
#                 hdul[0].header['CRPIX1'] += pad_p  # Update reference pixel
#                 if impar:
#                     hdul[0].header['NAXIS1'] += 1
#                 padded_data = np.pad(image_data, ((0, 0), (pad_p, pad_p + impar)), mode='constant', constant_values=0)
#                 
#             elif ejes[1] < ejes[0]:  # NAXIS2 is smaller, pad vertically
#                 new_naxis2 = ejes[0]
#                 hdul[0].header['NAXIS2'] = new_naxis2
#                 hdul[0].header['CRPIX2'] += pad_p  # Update reference pixel
#                 if impar:
#                     hdul[0].header['NAXIS2'] += 1
#                 padded_data = np.pad(image_data, ((pad_p, pad_p + impar), (0, 0)), mode='constant', constant_values=0)
# 
#             # Update data and save
#             hdul[0].data = padded_data
#             hdul.writeto(folder + f'PADDED_{f_file}', overwrite=True)
#             print(f'New axes: {hdul[0].header["NAXIS1"]}, {hdul[0].header["NAXIS2"]}')
#             print(30 * '_')
#         
#         else:
#             print(f'Axes are already equal in {f_file}')
# =============================================================================
# %%
# size_ax = np.array([])
# for chip in range(1, 2):
#     # folder = '/home/data/alvaro/gns_test/F%s/SWarp/chip2%s/'%(field, chip)
#     fits_files = [f for f in sorted(os.listdir(folder)) if f.endswith('.fits') and f.startswith('70_pointings_f%s_c%s' % (field, chip))]
#     for nf, f_file in enumerate(fits_files):
#         hdul = fits.open(folder + f_file)
#         image_data = hdul[0].data
#         size_ax = np.append(size_ax,hdul[0].header['NAXIS1'], hdul[0].header['NAXIS2'])
#         size_ax = np.unique(size_ax)
#         ejes = [hdul[0].header['NAXIS1'], hdul[0].header['NAXIS2']]
#         print(f'Working with {f_file}')
#         print(f'Original axes: {ejes}')
# unique_ax = int(np.max(size_ax))
# unique_ax = int(np.max(size_ax))
# sys.exit()






