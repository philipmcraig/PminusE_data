# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:24:28 2015

@author: np838619


White Sea configuration, whitecconfig.py. To be read by arctic_click.py.
The White Sea is part of the Arctic Ocean by Archangel. This needs to be defined
as it is south of 65.5N.

bounds array contains the northern and southern boundaries of the required ocean
mask. Contains one latitude and two longitudes in each sub-array.

The basemap is also defined, using the Miller projection
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid

bounds = pl.array([[65.5,32.0,41.0],[63.7,32.0,41.0]])

m = Basemap(lon_0 = lons.mean(),lat_0=lats.mean(),lat_ts=20.,
            llcrnrlon=27.,llcrnrlat=62.,urcrnrlon=46.,urcrnrlat=70.,
            resolution='l',projection='mill',suppress_ticks=True)
m.drawcoastlines()