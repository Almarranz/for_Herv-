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

# %

folder = '/Users/amartinez/Desktop/for_people/for_Herve/gns2/F20/'
pruebas = '/Users/amartinez/Desktop/for_people/for_Herve/pruebas/'
chip = 1

field = 20

# Define the input files
fits_files = [f for f in os.listdir(folder) if f.endswith('.fits')]
output_file = pruebas + 'combined_cube.fits'


hdu_list = fits.open(folder + fits_files[0])
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

# # corners = np.array([[0, 0], [naxis1, 0], [0, naxis2], [naxis1, naxis2]])
# corners = np.array([[-x_off, 0], [-x_off,naxis1], [naxis2-x_off,0], [naxis2-x_off, naxis1]])
# # corners = np.array([[min(x), min(y)], [max(x), min(y)], [min(x), max(y)], [max(x), max(y)]])
# world_corners = wcs.all_pix2world(corners, 1)

# # %
# # Extract RA and Dec from world_corners
# ra_corners = world_corners[:, 0]
# dec_corners = world_corners[:, 1]
# # %
# # color = 'red'
# # with open(pruebas+ 'corners.reg', 'w') as f:
# #     f.write('# Region file format: DS9 version 4.1'+"\n"+'global color=%s dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'%(color)+"\n"+'fk5'+'\n')
# #     f.close


# # for i in range(len(ra_corners)):
# #     if i%1 == 0:
# #         with open(pruebas+ 'corners.reg',  'a') as f:
# #             f.write('\n'.join(('point(%s,%s) # point=circle'%(ra_corners[i],dec_corners[i]),'\n')))   
# #             # print('ssss')    
# #         f.close
# # %
# # Convert RA/Dec to Galactic Coordinates using SkyCoord
# skycoord_corners = SkyCoord(ra=ra_corners*u.degree, dec=dec_corners*u.degree, frame='icrs')
# galactic_corners = skycoord_corners.galactic


# # Create a mask to filter vvv.txt data
# mask = (vvv_gal.l >= min(galactic_corners.l)) & (vvv_gal.l <= max(galactic_corners.l)) & (vvv_gal.b >= min(galactic_corners.b)) & (vvv_gal.b <= max(galactic_corners.b))

# Apply the mask to the vvv.txt data
# vvv_overlap = vvv_data[mask]

vvv_overlap.sort('J')

x_vvv = vvv_overlap['x']
y_vvv = vvv_overlap['y']
# xh = naxis1/2
yh = naxis1/2

xh = (max(x_vvv) + min(x_vvv))/2
yh = (max(y_vvv) + min(y_vvv))/2
chip= 4
#crop list for each chip
if (chip == 1):
    idx = np.nonzero((x_vvv < xh) & (y_vvv < yh))
elif (chip == 2):
    idx = np.nonzero((x_vvv > xh) & (y_vvv < yh))
elif (chip == 3):
    idx = np.nonzero((x_vvv > xh) & (y_vvv > yh))
elif (chip == 4):
    idx = np.nonzero((x_vvv < xh) & (y_vvv > yh))
# if (chip == 1):
#     idx = np.nonzero((x_vvv < xh) & (y_vvv < yh))
# elif (chip == 2):
#     idx = np.nonzero((x_vvv > xh) & (y_vvv < yh))
# elif (chip == 3):
#     idx = np.nonzero((x_vvv > xh) & (y_vvv > yh))
# elif (chip == 4):
#     idx = np.nonzero((x_vvv < xh) & (y_vvv > yh))
# # %%
fig, ax = plt.subplots(1,1)
ax.scatter(x_vvv,y_vvv)
ax.scatter(x_vvv[idx],y_vvv[idx])
ax.axvline(xh, color = 'r')
# # %%



vvv_overlap[idx][0:1000].write(pruebas + 'vvv_test_c%s.txt'%(chip), format = 'ascii', overwrite = True)
# vvv_overlap.write(pruebas + 'vvv_test_all.txt', format = 'ascii', overwrite = True)

# %%

# gns = Table.read(folder + 'dejitter_stars_1_53.txt', format = 'ascii')

# x_gns = gns['x']
# y_gns = gns['y']

# xy_gns = np.c_[x_gns,y_gns]

# xy_vvv = np.vstack((vvv_overlap['x'],vvv_overlap['y'])).T

# p, (pos_img, pos_img_t) = aa.find_transform(xy_gns, xy_vvv, max_control_points=200)






