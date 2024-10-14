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
from astropy.table import Table
from astropy.wcs.utils import fit_wcs_from_points
from astropy.wcs import WCS
from astropy.io import fits
import astroalign as aa
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.wcs.utils import fit_wcs_from_points

# %
# %%
folder = '/Users/amartinez/Desktop/for_people/for_Herve/gns2/F20/'
pruebas = '/Users/amartinez/Desktop/for_people/for_Herve/pruebas/'
sf_folder = '/Users/amartinez/Desktop/for_people/for_Herve/gns2/F20/stars_lists/'
###########IMPORTANT################
# list of raw images is in /home/data/raw/GNS_2/H/Field/20
####################################
field = 20
# chip = 2
# Define the input files
for chip in range(1,2):
    fits_files = [f for f in os.listdir(folder) if f.endswith('.fits') and f.startswith('HAWKI')]
    # fits_files = ['HAWKI.2022-08-02T00:54:39.570.fits','HAWKI.2022-08-02T00:56:08.395.fits']
    primary_hdu = fits.PrimaryHDU()
    list_of_hdu = [primary_hdu]
    for nf,f_file in enumerate(fits_files):
        print(30*'*')
        print(f'Working on file {nf+1}/{len(fits_files)}')
        print(30*'*')
        # Open the file and read all lines
        with open(folder + 'list.txt', 'r') as file:
            lines = file.readlines()
        
        
        # Search for the FITS file in the lines
        for sf_i, line in enumerate(lines, 1):  # enumerate starting from 1 to get line numbers
            if f_file in line:
                gns = Table.read(sf_folder + 'dejitter_stars_%s_%s.txt'%(chip,sf_i), format = 'ascii')
                print(f'Found "{f_file}" on line {sf_i}')
                break
        else:
            print(f'"{f_file}" not found in list.txt')
            sys.exit('YOOOOMAMMAAAAAAA')
       
    
        # %%
        
        
        hdu_list = fits.open(folder + f_file)
        image_data = hdu_list[0].data
        wcs = WCS(hdu_list[0].header, naxis = 2)
        # wcs = WCS(hdu_list[0].header)
        naxis1 = hdu_list[0].header['NAXIS1']
        naxis2 = hdu_list[0].header['NAXIS2']
        x_off = hdu_list[0].header['CRPIX2']
        y_off = hdu_list[0].header['CRPIX1']
        
        # stars_c1 = Table.read(folder + 'stars_calibrated_H_chip1.txt',
        #                  names = ('ra','dec','x','y','f','H','dx','dy','df','dH'), format = 'ascii')
        # stars= Table.read(folder + 'dejitter_common_1_53.txt',
        #                       names = ('x','y'), format = 'ascii')
        vvv = Table.read('/Users/amartinez/Desktop/PhD/Catalogs/VVV/b333/PMS/b333.dat', format = 'ascii')
        # %%
       
        use_idx = vvv['J']<900
        vvv_data = vvv[use_idx]
        
        vvv_ra = vvv_data['ra']
        vvv_dec = vvv_data['dec']
        vvv_gal = SkyCoord(ra = vvv_ra, dec = vvv_dec, unit = 'degree').galactic 
        
        xy = wcs.all_world2pix(np.c_[vvv_data['ra'],vvv_data['dec']],1)
        x = xy[:,0]
        y = xy[:,1]
        
        
        vvv_data.add_columns((xy[:,0],xy[:,1]),names = ('x','y'))
        
        mask = ((x >= -x_off) & (x < naxis1-x_off) & (y >= 0) & (y < naxis2))
        
        # Apply the mask to the vvv.txt data
        vvv_overlap = vvv_data[mask]
        
        
        vvv_overlap.sort('J')
        
        x_vvv = vvv_overlap['x']
        y_vvv = vvv_overlap['y']
        # xh = naxis1/2
        
        
        xh = (max(x_vvv) + min(x_vvv))/2
        yh = (max(y_vvv) + min(y_vvv))/2
        #crop list for each chip
        if (chip == 1):
            idx = np.nonzero((x_vvv < xh) & (y_vvv < yh))
        elif (chip == 2):
            idx = np.nonzero((x_vvv > xh) & (y_vvv < yh))
        elif (chip == 3):
            idx = np.nonzero((x_vvv > xh) & (y_vvv > yh))
        elif (chip == 4):
            idx = np.nonzero((x_vvv < xh) & (y_vvv > yh))
        #
        fig, ax = plt.subplots(1,1)
        ax.scatter(x_vvv,y_vvv)
        ax.scatter(x_vvv[idx],y_vvv[idx])
        ax.axvline(xh, color = 'r')
        # # %%
        
        
        
        vvv_overlap[idx][0:1000].write(pruebas + 'vvv_test_c%s.txt'%(chip), format = 'ascii', overwrite = True)
        # vvv_overlap.write(pruebas + 'vvv_test_all.txt', format = 'ascii', overwrite = True)
        
        #
        
        #This ara the xy coordinates from starfinder on the original cubes
        x_gns = gns['x']
        y_gns = gns['y']
        
        xy_gns = np.c_[x_gns,y_gns]
        
        xy_vvv = np.vstack((vvv_overlap['x'][idx],vvv_overlap['y'][idx])).T
        # Astroalign does not work when repeated value appears in the same array
        xy_test, idx_u, = np.unique(xy_vvv, axis=0, return_index=True)
        xy_vvv_unique = xy_vvv[np.sort(idx_u)]
        # sys.exit(154)
        # fig, ax = plt.subplots(1,1)
        # ax.scatter(xy_vvv[:,0],xy_vvv[:,1], s =20)
        # ax.scatter(xy_vvv_unique[:,0],xy_vvv_unique[:,1], s =1)
        
        # sys.exit(154)
        p, (pos_img, pos_img_t) = aa.find_transform(xy_gns, xy_vvv_unique, max_control_points=200)
        
        # %%
        fig, (ax1,ax2) = plt.subplots(1,2)
        ax1.scatter(pos_img[:,0],pos_img[:,1])
        ax2.scatter(pos_img_t[:,0],pos_img_t[:,1])
        # %%
        # %
        # Find the common VVV stars in the VVV table to retrive the RA and Dec coordinates
        x_vvv = pos_img_t  # 
        
    
        x_overlap = np.array(vvv_overlap['x'][idx])
        y_overlap = np.array(vvv_overlap['y'][idx])
        
    
        indices = []
        
    
        for coord in x_vvv:
            x, y = coord
            # Find the indices where both x and y match in vvv_overlap
            index = np.where((x_overlap == x) & (y_overlap == y))[0]
            
            # If a match is found, append the index
            if index.size > 0:
                indices.append(index[0])  # Assuming one-to-one match
            else:
                indices.append(-1)  # Append -1 for unmatched points
        
        # Convert indices list to numpy array for easier handling
        indices = np.array(indices)
        
        # Print or return the indices
        print(indices)
        vvv_com = vvv_overlap[idx][indices]
        vvv_com.write(pruebas + 'vvv_common_c%s.txt'%(chip), format = 'ascii', overwrite = True)
        
        # %%
        gns_com_xy = Table(pos_img, names = ('x','y'))
        vvv_com_xy = Table(pos_img_t, names = ('x','y'))
        gns_com_xy.write(pruebas + 'gns_c%s_com.txt'%(chip), format = 'ascii',overwrite=True)
        vvv_com_xy.write(pruebas + 'vvv_c%s_com.txt'%(chip), format = 'ascii',overwrite=True)
        # %%
        # vvv_com_ad = np.c_[vvv_com['ra'],vvv_com['dec']]
        vvv_com_ad = SkyCoord(ra = vvv_com['ra'],dec = vvv_com['dec'],unit = 'degree')
        xy_com = np.vstack((pos_img[:,0],pos_img[:,1]))
        wcs_new = fit_wcs_from_points(xy_com, vvv_com_ad, projection="TAN")
        
        # %%
         
        # Create a primary HDU (the first extension should always be the primary HDU)
        
        
        # Loop through the FITS files
        
            # for n,fil in enumerate(fits_files):
    
            # Open the FITS file
        # There are 8 images in each file; we take the first 7
        data_cube = hdu_list[0].data
        header = hdu_list[0].header  # Get the header from the corresponding extension
        for i in range(7):
            # Get the data and header for each image (extension)
            if chip ==1:
                data = data_cube[i][0:2048,0:2048]# Use i+1 since the first is the PrimaryHDU
               
            if chip == 2:
                data = data_cube[i][0:2048,2048:]
                
            if chip == 3:
                data = data_cube[i][2048:,2048:]
               
            if chip == 4:
                data = data_cube[i][2048:,0:2048]
    
            
            
            wcs_header = wcs_new.to_header()
            
            for card in wcs_header.cards:
                # If the key already exists, it will be replaced, otherwise it will be added
                header[card.keyword] = card.value
            # Create an ImageHDU with the data and the corresponding header
            # image_hdu = fits.ImageHDU(data=data, header=wcs_header)
            image_hdu = fits.ImageHDU(data=data, header=header)
            
            # Add the new ImageHDU to the HDUList
            list_of_hdu.append(image_hdu)
    
            # Create the HDUList object
            combined_hdul = fits.HDUList(list_of_hdu)
            
        # Write to a new file, removing the last slice from each original cube (so 7 from each)
    # combined_hdul.writeto(pruebas  +'test_f20_c%s.fits'%(chip), overwrite=True)
    del combined_hdul[35]  
    del combined_hdul[34]  
    del combined_hdul[33]  
    del combined_hdul[14]  # Extension 15 is at index 14
    del combined_hdul[7]   # Extension 8 is at index 7
    del combined_hdul[6]   # Extension 7 is at index 6 (since we just deleted the previous extension 7)

    combined_hdul.writeto(pruebas  +'%s_pointings_f20_c%s.fits'%(len(fits_files),chip), overwrite=True)
    print('FITS file F%sc%s created.'%(field, chip))


# # %%
# # Updated originals headers

# for card in wcs_header.cards:
#     print('NEW')
#     print(card[0],card[1])
    
#     try:
#         print(card[0],header[card[0]])
#         print(10*'_')
#     except:
#         print(f'None {card[0]}')
#         print(10*'_')
# # %%
# for card in wcs_header.cards:
#     # If the key already exists, it will be replaced, otherwise it will be added
#     header[card.keyword] = card.value














