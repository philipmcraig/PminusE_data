# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:24:28 2015

@author: np838619

Baltic Sea configuration, balticconfig.py. To be read by click2.py.

bounds array contains the northern and southern boundaries of the required ocean
mask. Contains one latitude and two longitudes in each sub-array.

The basemap is also defined, using the Miller projection
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid

bounds = pl.array([[66.0,21.5,26.0],[53.5,9.0,20.0]])

m = Basemap(llcrnrlon=0.,llcrnrlat=50.,urcrnrlon=37.,urcrnrlat=67.,
            lon_0 = lons.mean(),lat_0=lats.mean(),lat_ts=20.,
resolution='l',projection='mill',suppress_ticks=True)
m.drawcoastlines()