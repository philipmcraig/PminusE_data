# -*- coding: utf-8 -*-
"""
Created on Mon May  4 15:17:29 2015

@author: np838619

Code for reading in monthly_mean E-P data from monthly_means.nc and plotting
each month
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset, MFDataset
import os
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.colors import Normalize
from matplotlib import ticker
import scipy.ndimage as ndimage
from read_ncfiles import MidpointNormalize

ncfile = Dataset('monthly_means.nc','r')
lat = ncfile.variables['lat']
lon = ncfile.variables['lon']
EmP = ncfile.variables['EmP']
#ncfile.close()

lon_0 = lon.mean()
lat_0 = lat.mean()

pl.subplot(621) #6 rows, 2 columns, 1st plot
m = Basemap(projection='merc',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20)
