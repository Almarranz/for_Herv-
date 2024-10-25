#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 11:35:06 2024

@author: amartinez
"""

# Modify the size of the axii from the output of SWarp to make then equal

import numpy as np
from astropy.io import fits
import sys

pruebas = '/Users/amartinez/Desktop/for_people/for_Herve/pruebas/'
slc = 238
hdul = fits.open(pruebas + '70_pointings_f20_c2.%04d.resamp.fits'%(slc))


# Access the image data
image_data = hdul[0].data


ejes = [ hdul[0].header['NAXIS1'],hdul[0].header['NAXIS2'] ]
print(ejes)
if hdul[0].header['NAXIS2']  != hdul[0].header['NAXIS1']:
    print('yommma')
    diff = abs(hdul[0].header['NAXIS1']-hdul[0].header['NAXIS2'])
    if diff%2 !=0:
        diff += 1
    pad_p = int(diff/2)
    # menor = min(hdul[0].header['NAXIS1'],hdul[0].header['NAXIS1'])
    menor = np.where(ejes ==np.min(ejes))[0]
    # ejes[menor] += diff
    if menor == 0:
        hdul[0].header['NAXIS1'] += diff
        hdul[0].header['CRPIX1'] += pad_p
        padded_data = np.pad(image_data, ((0,0), (pad_p, pad_p)), mode='constant', constant_values=0)
    elif menor == 1:    
        hdul[0].header['NAXIS2'] += diff
        hdul[0].header['CRPIX2'] += pad_p
        padded_data = np.pad(image_data, ((pad_p, pad_p), (0, 0)), mode='constant', constant_values=0)
   
    hdul[0].data = padded_data
# sys.exit()
 

# print(hdul[0].header['NAXIS2'] )
# print(hdul[0].header['NAXIS1'] )

# # Pad the image by adding rows of zeros equally at the top and bottom
# padded_data = np.pad(image_data, ((1, 1), (0, 0)), mode='constant', constant_values=0)

# # Update the data in the HDU

# # Modify the header to reflect the padding
# hdul[0].header['NAXIS2'] += 2  # Update the number of pixels along axis 2
# hdul[0].header['CRPIX2'] += 1  # Shift the reference pixel for the added rows
hdul.writeto(pruebas + 'padded_fits_%04d.fits'%(slc), overwrite=True)
