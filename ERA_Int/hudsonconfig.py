# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:24:28 2015

@author: np838619

Hudson Bay configuration, hudsoncconfig.py. To be read by click2.py.

bounds array contains the northern and southern boundaries of the required ocean
mask. Contains one latitude and two longitudes in each sub-array.

The basemap is also defined, using the Miller projection
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid

bounds = pl.array([[65.5,-96.0,-64.0],[51.0,-96.0,-64.0]])

m = Basemap(lon_0 = lons.mean(),lat_0=lats.mean(),lat_ts=20.,
            llcrnrlon=-100.,llcrnrlat=45.,urcrnrlon=-55.,urcrnrlat=70.,
            resolution='l',projection='mill',suppress_ticks=True)
m.drawcoastlines()