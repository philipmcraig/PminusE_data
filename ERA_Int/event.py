# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 09:26:16 2015

@author: np838619
"""

import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.colors import Normalize
from netCDF4 import Dataset
from read_ncfiles import MidpointNormalize
from functions import *


def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    #print 'x = %d, y = %d'%(
    #    ix, iy)

    global coords
    coords.append((ix, iy))
    

    if len(coords) == 6:
        fig.canvas.mpl_disconnect(cid)
        for i in range(len(coords)):
            print coords[i]

    return #coords


#x = pl.arange(-10,10)
#y = x**2

#take longitude and latitude from an ERA-Interim netCDF file
nc1 = Dataset('tcdq_ann_mean_ERAI.nc','r')
lons = nc1.variables['lon'][:]
lats = nc1.variables['lat'][:]
mask = nc1.variables['LSM'][:]
#tcdq = ncfile.variables['tcdq'][:]
nc1.close()

fig = pl.figure()
#ax = fig.add_subplot(111)
#ax.plot(x,y)

#coords = []


#cid = fig.canvas.mpl_connect('button_press_event', onclick)



lon_0 = lons.mean()
lat_0 = lats.mean()

#pl.figure()
# change lat & lon extents manually for each map
m = Basemap(llcrnrlon=90.,llcrnrlat=-35.,urcrnrlon=-65.,urcrnrlat=50.,lon_0=-145.,
                o_lat_p=90.,o_lon_p=180., # remove these for normal projections
                lat_ts=20.,resolution='l',projection='rotpole',suppress_ticks=True)
mask, lons = shiftgrid(180.0, mask, lons, start=False)
lon, lat = pl.meshgrid(lons,lats)
X, Y = m(lon,lat)
cmap = pl.get_cmap('Purples_r')
#cs=m.pcolormesh(X,Y,pl.squeeze(newmask),cmap=cmap)
m.drawcoastlines(color='k')
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='white')
m.drawparallels(lats,labels=[0,0,0,0],linewidth=0.5,ax=None)
m.drawmeridians([-90,0,90],labels=[1,0,0,1],linewidth=0,ax=None)

# specify the northern & southern boundaries required
# change them manually, not sure if I can be bothered making this interactive
nrth_bndry_lat = [lats[NearestIndex(lats,47.00)],lats[NearestIndex(lats,47.00)]]
nrth_bndry_lon = [lons[NearestIndex(lons,136.00)],lons[NearestIndex(lons,-122.00)]]
sth_bndry_lat =[lats[NearestIndex(lats,-30.00)],lats[NearestIndex(lats,-30.00)]]
sth_bndry_lon =[lons[NearestIndex(lons,150.00)],lons[NearestIndex(lons,-71.00)]]

# plot the northern & southern boundaries
x1,y1 = m(nrth_bndry_lon,nrth_bndry_lat)
m.plot(x1,y1,'D-',color='r',linewidth=2)
x2,y2 = m(sth_bndry_lon,sth_bndry_lat)
m.plot(x2,y2,'D-',color='r',linewidth=2)


coords = []
cid = fig.canvas.mpl_connect('button_press_event', onclick)
boundary = pl.zeros([len(coords),2])
boundary[:] = coords[:]
ilon, ilat = m(boundary[:,0],boundary[:,1],inverse=True)
#print 'This is finished'
    

