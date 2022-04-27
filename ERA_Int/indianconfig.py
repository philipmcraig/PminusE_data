# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:24:28 2015

@author: np838619

Indian Ocean configuration, indianconfig.py. To be read by click2.py.

bounds array contains the northern and southern boundaries of the required ocean
mask. Contains one latitude and two longitudes in each sub-array.

The basemap is also defined, using the Miller projection.
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid

# 1st sub-array: northlat, northlon1, northlon2
# 2nd sub-array: southlat, southlon1, southlon2
bounds = pl.array([[30.5,32.0,50.0],[-35.0,20.0,117.0]])

m = Basemap(llcrnrlon=15.,llcrnrlat=-40.,urcrnrlon=145.,urcrnrlat=35.,
            lon_0 = lons.mean(),lat_0=lats.mean(),lat_ts=20.,
            resolution='l',projection='mill',suppress_ticks=True)
m.drawcoastlines()