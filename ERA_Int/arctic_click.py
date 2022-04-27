# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 16:08:31 2015

@author: np838619

Code to extract moisture flux divergence from the Arctic Ocean. Everything north
of 65.5N, plus everything north of 60N in the Atlantic sector, Hudson Bay and a
little bit of the Arctic Ocean near Archangel.

Last updated: 23/09/2015 4:15PM 3rd commit
"""

import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.colors import Normalize
from netCDF4 import Dataset
from functions import *


#take longitude and latitude from an ERA-Interim netCDF file
nc1 = Dataset('tcdq_ann_mean_ERAI.nc','r')
lons = nc1.variables['lon'][:]
lats = nc1.variables['lat'][:]
mask = nc1.variables['LSM'][:]
#tcdq = ncfile.variables['tcdq'][:]
nc1.close()

# annual mean moisture flux divergence:
nc2 = Dataset('tcdq_mean.nc','r')
tcdq = nc2.variables['tcdq'][:]
nc2.close()

newmask = 1-mask

# Before selecting additional regions with clicking first extract everything
# north of 65.5N:
ind65 = NearestIndex(lats,65.5) # Nearest Index in the ERA-Interim latitude array to 65.5N
arc65_tcdq = tcdq[0,:ind65,:]

# Define longitude boundaries & extract North Atlantic data between Canada & Norway
# Remember that the 0 meridian is in the way so the data must be extracted in 2 parts
# Longitude of Canadian coastline
canlon = 360 - 65.7 # ERA-Interim longitudes range from -180 to 180, Basemap goes from 0-360
canind = NearestIndex(lons,canlon) # Nearest index in the ERA-Interim longitude array
# LOngitude of Norwegian coastline:
norlon = 12.5
norind = NearestIndex(lons,norlon)
ind60 = NearestIndex(lats,60) # Nearest Index in the ERA-Interim latitude array to 60N

# Extract data for the N Atlantic 60-65.5N from Canada to the zero meridian:
naw_tcdq = tcdq[0,ind65:ind60,canind:]
# Extract data for the N Atlantic 60-65.5N from the zero meridian to Norway:
nae_tcdq = tcdq[0,ind65:ind60,:norind+1]

# set up empty array for arctic data:
arc_tcdq = pl.zeros_like(tcdq[0])

arc_tcdq[:ind65,:] = arc65_tcdq[:,:]
arc_tcdq[ind65:ind60,canind:] = naw_tcdq[:,:]
arc_tcdq[ind65:ind60,:norind+1] = nae_tcdq[:,:]

pl.figure(1)
exec(open('arcticconfig.py').read())
arctic, lons = shiftgrid(180.0, arc_tcdq*(1-mask[0]), lons, start=False)
lon, lat = pl.meshgrid(lons,lats)
X, Y = m(lon,lat)
cmap = pl.get_cmap('seismic')
norm = MidpointNormalize(midpoint=0)
cs = m.pcolormesh(X,Y,arctic,cmap=cmap,norm=norm)


#elif select == 6:
#    ocean = 'arctic'
#    flag = False

ocean = ''
flag = True
while flag == True:
    select = input("Which map is required?\n\nEnter number to select basin:\n\n " + 
                    "1 - Hudson Bay\n 2 - White Sea\n\n Enter anything else to exit\n\n")
    nmbr = isinstance(select,int) # is the input an integer?
    if nmbr == True:
        if select == 1:
            ocean = 'hudson'
            flag = False
        elif select == 2:
            ocean = 'white'
            flag = False
        else:
            ocean = 'none'
            flag = False
            # entering anything other than 1 or 2 ends the program 
    
ncfile = Dataset('tcdq_ann_mean_ERAI.nc','r') # really irritating that I need to do this
lons = ncfile.variables['lon'][:]
lats = ncfile.variables['lat'][:]
ncfile.close()
pl.figure(2)
exec(open(ocean + 'config.py').read())
nm, lons = shiftgrid(180.0, newmask, lons, start=False) # shifting the grid works for some reason
lon, lat = pl.meshgrid(lons,lats)
X, Y = m(lon,lat)
m.drawparallels(lats,labels=[0,0,0,0],linewidth=0.5,ax=None)
m.drawmeridians([-90,0,90],labels=[1,0,0,1],linewidth=0,ax=None)
pl.subplots_adjust(right=0.99,left=0.01,top=1.0,bottom=0.0)# use as much of the figure window as possible
    
# specify the northern & southern boundaries required from bounds array in config file:
nrth_bndry_lat = [lats[NearestIndex(lats,bounds[0,0])],lats[NearestIndex(lats,bounds[0,0])]]
nrth_bndry_lon = [lons[NearestIndex(lons,bounds[0,1])],lons[NearestIndex(lons,bounds[0,2])]]
sth_bndry_lat =[lats[NearestIndex(lats,bounds[1,0])],lats[NearestIndex(lats,bounds[1,0])]]
sth_bndry_lon =[lons[NearestIndex(lons,bounds[1,1])],lons[NearestIndex(lons,bounds[1,2])]]


# plot the northern & southern boundaries
x1,y1 = m(nrth_bndry_lon,nrth_bndry_lat)
m.plot(x1,y1,'D-',color='r',linewidth=2)
x2,y2 = m(sth_bndry_lon,sth_bndry_lat)
m.plot(x2,y2,'D-',color='r',linewidth=2)

#print 'Enter how many points are required'
#no_of_points = raw_input()
no_of_points = len(lats[NearestIndex(lats,bounds[0,0]):NearestIndex(lats,bounds[1,0])+1])
#no_of_points = NearestIndex(lats,-35.00) - NearestIndex(lats,45.00)
print 'Click ' + str(2*no_of_points+2) + ' points \nRight click to cancel last input, click the middle button to stop.'
coords = pl.ginput(n=2*no_of_points,timeout=0) # right click cancels last input


lonlat = MapCoords(coords,m)

a,b = m(lonlat[:,0],lonlat[:,1])
m.plot(a,b,linewidth=0,color='r',marker='.') # plot each click as a red dot


# find the nearest grid point indices on the ERA-Interim grid to the clicked points:
grid = GridPoints(lonlat,lons,lats)


# split the grid array into the western & eastern boundaries of the mask:
grid_west, grid_east = SplitGrid(grid,lonlat)

#write to file
#f = open('Atlantic_boundary.txt','w')
#f.write('Longitde, Latitude co-ordinates defining a boundary around the Atlantic Ocean basin. \n Used to extract E-P data. \n \n')
#pl.savetxt(f,lonlat,fmt='%9.5f')
#f.close()

# make array with points inside mask=1 & outside mask = 0:
eragrid = MaskOutsideGrid(grid_west,grid_east,lons,lats)


if int(select) == 1:
    # Need to multiply eragrid by the tcdq array
    # all required points as 1 X moisture flux divergence X weighting of each point
    hud_tcdq = eragrid.transpose()*tcdq[0]*newmask[0] 
    # make an array where everything required=1 & the rest is nan, USEFUL FOR WHEN SAVING IN NETCDF!!!!!!!
    # this is done before plotting to avoid shiftgrid issues & no need to store a further array at start of code
    newarray = BasinOnes(grid_west,grid_east,eragrid)
    # Make the final tcdq array
    # multiplying byarray of 1 & nan avoids problem of zeros which might interfere with calculations
    hudnan = hud_tcdq*newarray
    # Make the final mask for the ocean basin
    hud_fnlmask = BasinMask(hudnan)
elif int(select) == 2:
    whi_tcdq = eragrid.transpose()*tcdq[0]*newmask[0]
    newarray = BasinOnes(grid_west,grid_east,eragrid)
    whinan = whi_tcdq*newarray
    whi_fnlmask = BasinMask(whinan)





pl.figure(3)
ax2 = pl.subplot(111)
ncfile = Dataset('tcdq_ann_mean_ERAI.nc','r') # really irritating that I need to do this
lons = ncfile.variables['lon'][:]
lats = ncfile.variables['lat'][:]
ncfile.close()
# latitude & longitude extents should match those from figure 1
exec(open(ocean + 'config.py').read())
if int(select) == 1:
    hud2_tcdq, lons = shiftgrid(180., hud_tcdq, lons, start=False)
elif int(select) == 2:
    whi2_tcdq, lons = shiftgrid(180., whi_tcdq, lons, start=False)
lon, lat = pl.meshgrid(lons,lats)
X, Y = m(lon,lat)
m.fillcontinents(color='white',lake_color='white')
cmap = pl.get_cmap('seismic')
norm = MidpointNormalize(midpoint=0)
if int(select) == 1:
    cs = m.pcolormesh(X,Y,hud2_tcdq,cmap=cmap,norm=norm)
elif int(select) == 2:
    cs = m.pcolormesh(X,Y,whi2_tcdq,cmap=cmap,norm=norm)
m.drawcoastlines(color='k')
pl.subplots_adjust(right=0.99,left=0.01,hspace=-1.35)


# Uncomment and run this to show everything:
#arc2_tcdq = arc_tcdq + hud_tcdq + whi_tcdq
#pl.figure(4)
#ncfile = Dataset('tcdq_ann_mean_ERAI.nc','r')
#lons = ncfile.variables['lon'][:]
#lats = ncfile.variables['lat'][:]
#ncfile.close()
#exec(open('arcticconfig.py').read())
#A2, lons = shiftgrid(180.0, arc2_tcdq*(1-mask[0]), lons, start=False) # shifting the grid works for some reason
#lon, lat = pl.meshgrid(lons,lats)
#X, Y = m(lon,lat)
#cs = m.pcolormesh(X,Y,A2,cmap=cmap,norm=norm)

eragrid = pl.zeros_like(tcdq[0])
q = hud_fnlmask + whi_fnlmask
p = pl.where(hud_fnlmask==1.); l = pl.where(whi_fnlmask==1.)
eragrid[p[0][-1]+1:,:] = pl.float64('nan')
for i in range(ind65,p[0][-1]):
    for j in range(q.shape[1]):
        if q[i,j] == 0:
            eragrid[i,j] = pl.float64('nan')
for i in range(ind65,ind60):
    eragrid[i,p[1][-1]:] = pl.float64('nan')
pl.imshow(eragrid)