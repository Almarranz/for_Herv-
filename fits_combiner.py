#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 15:40:19 2024

@author: amartinez
"""

import astropy.io.fits as fits
import numpy as np
import os
import sys
# %%

folder = '/Users/amartinez/Desktop/for_people/for_Herve/gns2/F20/'
pruebas = '/Users/amartinez/Desktop/for_people/for_Herve/pruebas/'


field = 20

# %%
from astropy.io import fits

# Define the input files
fits_files = [f for f in os.listdir(folder) if f.endswith('.fits')]
output_file = pruebas + 'combined_cube.fits'

# Create a primary HDU (the first extension should always be the primary HDU)
primary_hdu = fits.PrimaryHDU()
list_of_hdu = [primary_hdu]
dic_header = {} 
dic_header_ = {}
# Loop through the FITS files
for n,fil in enumerate(fits_files):
    # Open the FITS file
    hdulist = fits.open(folder + fil)
    # There are 8 images in each file; we take the first 7
    data_cube = hdulist[0].data
    header = hdulist[0].header  # Get the header from the corresponding extension
    dic_header['h%s'%(n+1)] = header
    for i in range(7):
        # Get the data and header for each image (extension)
        data = data_cube[i]# Use i+1 since the first is the PrimaryHDU
        
        # Create an ImageHDU with the data and the corresponding header
        image_hdu = fits.ImageHDU(data=data, header=header)
        dic_header_['h%s'%(i)] = header
        # Add the new ImageHDU to the HDUList
        list_of_hdu.append(image_hdu)

# Create the HDUList object
combined_hdul = fits.HDUList(list_of_hdu)

# Write to a new file, removing the last slice from each original cube (so 7 from each)
combined_hdul.writeto(pruebas  +'F20_combined.fits', overwrite=True)

print('FITS file with 21 extensions created, each with corresponding headers.')

# %%
# Define the input files
fits_files = [f for f in os.listdir(folder) if f.endswith('.fits')]
output_file = pruebas + 'combined_cube.fits'

for ch in range(1,5):
# Create a primary HDU (the first extension should always be the primary HDU)
    primary_hdu = fits.PrimaryHDU()
    list_of_hdu = [primary_hdu]
    dic_header = {} 
    dic_header_ = {}
# Loop through the FITS files

    for n,fil in enumerate(fits_files):
        # Open the FITS file
        hdulist = fits.open(folder + fil)
        # There are 8 images in each file; we take the first 7
        data_cube = hdulist[0].data
        header = hdulist[0].header  # Get the header from the corresponding extension
        dic_header['h%s'%(n+1)] = header
        for i in range(7):
            # Get the data and header for each image (extension)
            if ch ==1:
                data = data_cube[i][0:2048,0:2048]# Use i+1 since the first is the PrimaryHDU
            if ch == 2:
                data = data_cube[i][0:2048,2048:]
            if ch == 3:
                data = data_cube[i][2048:,2048:]
            if ch == 4:
                data = data_cube[i][2048:,0:2048]
            
            
            # Create an ImageHDU with the data and the corresponding header
            image_hdu = fits.ImageHDU(data=data, header=header)
            dic_header_['h%s'%(i)] = header
            # Add the new ImageHDU to the HDUList
            list_of_hdu.append(image_hdu)
    
    # Create the HDUList object
    combined_hdul = fits.HDUList(list_of_hdu)
    
    # Write to a new file, removing the last slice from each original cube (so 7 from each)
    combined_hdul.writeto(pruebas  +'F20_combined_c%s.fits'%(ch), overwrite=True)

    print('FITS file F%sc%s created.'%(field, ch))

