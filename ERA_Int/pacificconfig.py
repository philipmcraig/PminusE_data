# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:24:28 2015

@author: np838619

Pacific Ocean configuration, pacificconfig.py. To be read by click2.py.

bounds array contains the northern and southern boundaries of the required ocean
mask. Contains one latitude and two longitudes in each sub-array.

The basemap is also defined, using the rotated pole projection.
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid

# 1st sub-array: northlat,northlon1, northlon2
# 2nd sub-array: southlat, southlon1, southlon2
bounds = pl.array([[65.5,187.0,-166.0],[-30.0,150.0,-70.0]])

m = Basemap(llcrnrlon=90.,llcrnrlat=-35.,urcrnrlon=-65.,urcrnrlat=70.,
            lon_0 = -145.,o_lat_p=90.,o_lon_p=180.,lat_ts=20., # o_*_p only used in rotpole
            resolution='l',projection='rotpole',suppress_ticks=True)
m.drawcoastlines()