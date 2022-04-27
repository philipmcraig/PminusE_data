# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:24:28 2015

@author: np838619

Arctic Ocean configuration, arcticcconfig.py. To be read by click2.py.

bounds array contains the northern and southern boundaries of the required ocean
mask. Contains one latitude and two longitudes in each sub-array.

The basemap is also defined, using the Miller projection
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid

bounds = pl.array([[47.5,28.0,40.5],[30.0,13.0,22.0]])

m = Basemap(lon_0 = 270.,lat_0=lats.mean(),lat_ts=20.,
            llcrnrlon=-180.,llcrnrlat=50.,urcrnrlon=180.,urcrnrlat=90.,
            resolution='l',projection='mill',suppress_ticks=True)
m.drawcoastlines()