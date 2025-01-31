#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 12:06:39 2022

@author: amartinez
"""

# =============================================================================
# Creates a .reg file in order to check cluster members in DS9
# =============================================================================




import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import QTable
from matplotlib import rcParams
import os
import glob
import sys
import math
from astropy.table import Table

# pruebas = '/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_off/pruebas/'
# pruebas = '/Users/amartinez/Desktop/PhD/HAWK/GNS_pm_scripts/clusters_regions/'

folder = '/Users/amartinez/Desktop/for_people/for_Herve/gns2/F20/'
pruebas = '/Users/amartinez/Desktop/for_people/for_Herve/pruebas/'

# name = 'vvv_all_xy.reg'

chip = 3
# name = 'vvv_cut_c%s_ad.reg'%(chip)
# name = 'vvv_test_c%s_ad.reg'%(chip)
# name = ' vvv_c%s_com.reg'%(chip)
# name = ' vvv_c%s_com_table.reg'%(chip)
name = 'stars_c%s_cut.reg'%(chip)
color = 'blue'

# stars = Table.read(pruebas + 'vvv_test_c%s.txt'%(chip), format = 'ascii')
# stars = Table.read(pruebas +  'vvv_c%s_com.txt'%(chip), format = 'ascii')
# stars = Table.read(pruebas + 'vvv_test_all.txt', format = 'ascii')
stars= Table.read(folder + 'dejitter_stars_%s_53.txt'%(chip), format = 'ascii')
# stars = Table.read(folder + 'stars_1.txt',format = 'ascii')
x,y = stars['x'][0:500], stars['y'][0:500]
# x,y = stars['x'], stars['y']
# ra,dec = stars['ra'], stars['dec']
# ra,dec = stars['ra'][0:500], stars['dec'][0:500]
# %%
# with open(pruebas+ name, 'w') as f:
#     f.write('# Region file format: DS9 version 4.1'+"\n"+'global color=%s dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'%(color)+"\n"+'fk5'+'\n')
#     f.close


# for i in range(len(ra)):
#         with open(pruebas+ name, 'a') as f:
#             f.write('\n'.join(('point(%s,%s) # point=x'%(ra[i],dec[i]),'\n')))   
#             # print('ssss')    
#         f.close
    

with open(pruebas+ name, 'w') as f:
    f.write('# Region file format: DS9 version 4.1'+"\n"+'global color=%s dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'%(color)+"\n"+'physical'+'\n')
    f.close

for i in range(len(x)):
# for i in range(50):
    if i%1 == 0:
        with open(pruebas+ name, 'a') as f:
            f.write('\n'.join(('point(%s,%s) # point=x'%(x[i],y[i]),'\n')))   
            # print('ssss')    
        f.close
    
    
    
    