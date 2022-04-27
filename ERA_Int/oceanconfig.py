# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:24:28 2015

@author: np838619

Ocean configurations
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid


m = Basemap(llcrnrlon=15.,llcrnrlat=-40.,urcrnrlon=145.,urcrnrlat=35.,lon_0 = lons.mean(),lat_0=lats.mean(),lat_ts=20.,resolution='l',projection='mill',suppress_ticks=True)
m.drawcoastlines()
#newmask, lons = shiftgrid(180.0, newmask, lons, start=False)