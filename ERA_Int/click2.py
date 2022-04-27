# -*- coding: utf-8 -*-
"""
Created on Tue May 05 16:23:07 2015

This code can be used to define boundaries for ocean basins NOT watersheds.
Click on map and the boundary array returns the lon-lat co-ordinates of each
point. Nearest indices on ERA-Interim grid of each co-ordinate determined & 
data within each region extracted, plotted & can be stored in a netCDF file.

@author: Philip Craig

Last updated: 24/01/16 5:56PM 18th commit
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

nc2 = Dataset('tcdq_mean.nc','r')
tcdq = nc2.variables['tcdq'][:]
nc2.close()

#newmask = pl.zeros_like(mask)
#for i in range(mask[0].shape[0]):
#    for j in range(mask[0].shape[1]):
#        if mask[0,i,j] != 0:
#            newmask[0,i,j] = 1
#nm2 = newmask # copy newmask array so it can be used again later as shiftgrid will mess it up
m2 = mask
# since 1 represents a cell 100% covered by land & 0 represents a cell 100% covered by water,
# 1-mask gives the correct weighting with which to multiply the moisture flux divergence
# might be ok to get rid of evrything .0.8 at some point
newmask = 1-mask
nm = newmask # duplicate newmask because shiftgrid will move move it 180degrees out of alignment with ERA-Interim longitudes


#----------------------------CALL UP MAP FOR CLICKING--------------------------

#lon_0 = lons.mean()
#lat_0 = lats.mean()

# Which map is required?
#print 'Which map is required?\n\nEnter number to select ocean:\n\n '
#print ' 0 - end\n 1 - Atlantic\n 2 - Pacific\n 3 - Indian'
#select = raw_input("Which map is required?\n\nEnter number to select ocean:\n\n " + 
#"0 - end\n 1 - Atlantic\n 2 - Pacific\n 3 - Indian\n\n")

ocean = ''
flag = True
while flag == True:
    select = input("Which map is required?\n\nEnter number to select ocean:\n\n " + 
                    "1 - Atlantic\n 2 - Pacific\n 3 - Indian\n 4 - Baltic\n" + 
                    " 5 - Mediterranean & Black\n\n Enter anything else to exit\n\n")
    nmbr = isinstance(select,int) # is the input an integer?
    #if nmbr == False:
    #    print 'Incorrect input! '
    if nmbr == True:
        if select == 1:
            ocean = 'atlantic'
            flag = False
        elif select == 2:
            ocean = 'pacific'
            flag = False
        elif select == 3:
            ocean = 'indian'
            flag = False
        elif select == 4:
            ocean = 'baltic'
            flag = False
        elif select == 5:
            ocean = 'med'
            flag = False
        else:
            ocean = 'none'
            #break # selecting 0 breaks out of the while loop and ends the program as a result
            flag = False

pl.figure(1)
ax1 = pl.subplot(111)
exec(open(ocean+'config.py').read())
newmask, lons = shiftgrid(180.0, newmask, lons, start=False) # shifting the grid works for some reason
lon, lat = pl.meshgrid(lons,lats)
X, Y = m(lon,lat)
#cmap = pl.get_cmap('Purples_r')
#cs=m.pcolormesh(X,Y,pl.squeeze(newmask),cmap=cmap)
#m.drawcoastlines(color='k')
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='white',zorder=0)
m.drawparallels(lats,labels=[0,0,0,0],linewidth=0.5,ax=None)
m.drawmeridians([-90,0,90],labels=[1,0,0,1],linewidth=0,ax=None)
pl.subplots_adjust(right=0.99,left=0.01,top=1.0,bottom=0.0)# use as much of the figure window as possible


# specify the northern & southern boundaries required
# change them manually, not sure if I can be bothered making this interactive
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
#NEED SOME WAY OF KEEPING COUNT OF CLICKS!!!!!
#print 'you clicked:',x
#pl.show()

#coords = []
#odds = [x for x in range(2*no_of_points) if x % 2 != 0]
#latbound = lats[NearestIndex(lats,47.00):NearestIndex(lats,-30.00)+1]
#for i in range(2*no_of_points):
#    
#    if i = no_of_points:
#        m.drawparallels([latbound[i-odds[i-no_of_points]-1]],linewidth=2,color='w',ax=None,zorder=15)
#        m.drawparallels([latbound[i-odds[i-no_of_points]]],linewidth=0.5,color='g',ax=None,zorder=20)
#    elif i > no_of_points:
#        m.drawparallels([latbound[i-odds[i-no_of_points]+1]],linewidth=2,color='w',ax=None,zorder=15)
#        m.drawparallels([latbound[i-odds[i-no_of_points]]],linewidth=0.5,color='g',ax=None,zorder=20)
#    else:
#        m.drawparallels([latbound[i-1]],linewidth=2,color='w',ax=None,zorder=10)
#        m.drawparallels([latbound[i]],linewidth=0.5,ax=None,zorder=0)
#    coords.append(pl.ginput(1,timeout=0,show_clicks=True))
#    time.sleep(0.25)


# convert the list of x,y figure co-ordinates to map long-lat co-ordinates:
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
    atl_tcdq = eragrid.transpose()*tcdq[0]*nm[0] 
    # make an array where everything required=1 & the rest is nan, USEFUL FOR WHEN SAVING IN NETCDF!!!!!!!
    # this is done before plotting to avoid shiftgrid issues & no need to store a further array at start of code
    newarray = Zero2Nan(eragrid.T)
    # Make the final tcdq array
    # multiplying byarray of 1 & nan avoids problem of zeros which might interfere with calculations
    atlnan = atl_tcdq*newarray
    # Make the final mask for the ocean basin
    atl_fnlmask = BasinMask(atlnan)
elif int(select) == 2:
    pac_tcdq = eragrid.transpose()*tcdq[0]*nm[0]
    newarray = Zero2Nan(eragrid.T)
    pacnan = pac_tcdq*newarray
    pac_fnlmask = BasinMask(pacnan)
elif int(select) == 3:
    ind_tcdq = eragrid.transpose()*tcdq[0]*nm[0]
    newarray = Zero2Nan(eragrid.T)
    indnan = ind_tcdq*newarray
    ind_fnlmask = BasinMask(indnan)
elif int(select) == 4:
    bal_tcdq = eragrid.transpose()*tcdq[0]*nm[0]
    newarray = Zero2Nan(eragrid.T)
    balnan = bal_tcdq*newarray
    bal_fnlmask = BasinMask(balnan)
elif int(select) == 5:
    med_tcdq = eragrid.transpose()*tcdq[0]*nm[0]
    newarray = Zero2Nan(eragrid.T)
    mednan = med_tcdq*newarray
    med_fnlmask = BasinMask(mednan)



#------------------------PLOT TCDQ FROM AREA OF INTEREST-----------------------

pl.figure(2)
ax2 = pl.subplot(111)
#pl.imshow(H)
ncfile = Dataset('tcdq_ann_mean_ERAI.nc','r') # really irritating that I need to do this
lons = ncfile.variables['lon'][:]
lats = ncfile.variables['lat'][:]
ncfile.close()
# latitude & longitude extents should match those from figure 1
exec(open(ocean+'config.py').read())
if int(select) == 1:
    atl_shft, newlons = shiftgrid(180., atl_tcdq, lons, start=False)
elif int(select) == 5:
    med_shft, newlons = shiftgrid(180., med_tcdq, lons, start=False)
else:
    newlons = lons
lon, lat = pl.meshgrid(newlons,lats)
X, Y = m(lon,lat)
m.fillcontinents(color='white',lake_color='white')
cmap = pl.get_cmap('seismic')
norm = MidpointNormalize(midpoint=0)
if int(select) == 1:
    cs = m.pcolormesh(X,Y,atl_shft,cmap=cmap,norm=norm)
elif int(select) == 2:
    cs = m.pcolormesh(X,Y,pac_tcdq,cmap=cmap,norm=norm)
elif int(select) == 3:
    cs = m.pcolormesh(X,Y,ind_tcdq,cmap=cmap,norm=norm)
elif int(select) == 4:
    cs = m.pcolormesh(X,Y,bal_tcdq,cmap=cmap,norm=norm)
elif int(select) == 5:
    cs = m.pcolormesh(X,Y,med_shft,cmap=cmap,norm=norm)
m.drawcoastlines(color='k')
pl.subplots_adjust(right=0.99,left=0.01,hspace=-1.35)




#--------------------!!!SAVE MASKED DATA IN NETCDF FILE!!!---------------------

newnc = Dataset('tcdq_basins_60.nc',mode='a',format='NETCDF4')

#lat_dim = newnc.createDimension('lat', 256)
#lon_dim = newnc.createDimension('lon', 512)
#lat_in = newnc.createVariable('lat', pl.float32, ('lat',))
#lat_in.units = 'degrees_north'
#lat_in.long_name = 'latitude'
#lon_in = newnc.createVariable('lon', pl.float32, ('lon',))
#lon_in.units = 'degrees_east'
#lon_in.long_name = 'longitude'
#lat_in[:] = lats # straight from ERA-Interim nc file
#lon_in[:] = lons # straight from ERA-Interim nc file
#
##time_dim = newnc.createDimension('time', None)
##time = newnc.createVariable('time', pl.float64, ('time',))
##time.long_name = 'time'
#
#Atlantic_TCDQ = newnc.createVariable('tcdq_Atl',pl.float64,('lat','lon'))
#Atlantic_TCDQ.units = 'kg m**-2 s**-1'
#Atlantic_TCDQ.standard_name = 'Atlantic Ocean (35S-60N) moisture flux divergence'
#Atlantic_TCDQ[:,:] = atlnan
#
#Atlantic_mask = newnc.createVariable('maskA',pl.float64,('lat','lon'))
#Atlantic_mask.units = '0 or 1'
#Atlantic_mask.standard_name = 'Atlantic Ocean (35S-60N) mask'
#Atlantic_mask[:,:] =  atl_fnlmask
#
#Indian_TCDQ = newnc.createVariable('tcdq_Ind',pl.float64,('lat','lon'))
#Indian_TCDQ.units = 'kg m**-2 s**-1'
#Indian_TCDQ.standard_name = 'Indian Ocean (>35S) moisture flux divergence'
#Indian_TCDQ[:,:] = indnan
#
#Indian_mask = newnc.createVariable('maskI',pl.float64,('lat','lon'))
#Indian_mask.units = '0 or 1'
#Indian_mask.standard_name = 'Indian Ocean (>35S) mask'
#Indian_mask[:,:] = ind_fnlmask 
#
#Pacific_TCDQ = newnc.createVariable('tcdq_Pac',pl.float64,('lat','lon'))
#Pacific_TCDQ.units = 'kg m**-2 s**-1'
#Pacific_TCDQ.standard_name = 'Pacific Ocean (30S-47N) moisture flux divergence'
#Pacific_TCDQ[:,:] = pacnan
#
#Pac_mask = newnc.createVariable('maskP',pl.float64,('lat','lon'))
#Pac_mask.units = '0 or 1'
#Pac_mask.standard_name = 'Pacific Ocean (30S-Bering Strait) mask'
#Pac_mask[:,:] = pac_fnlmask
#
#Baltic_TCDQ = newnc.createVariable('tcdq_Bal',pl.float64,('lat','lon'))
#Baltic_TCDQ.units = 'kg m**-2 s**-1'
#Baltic_TCDQ.standard_name = 'Baltic Sea moisture flux divergence'
#Baltic_TCDQ[:,:] = balnan
#
#Bal_mask = newnc.createVariable('maskB',pl.float64,('lat','lon'))
#Bal_mask.units = '0 or 1'
#Bal_mask.standard_name = 'Baltic Sea mask'
#Bal_mask[:,:] = bal_fnlmask

Med_TCDQ = newnc.createVariable('tcdq_Med',pl.float64,('lat','lon'))
Med_TCDQ.units = 'kg m**-2 s**-1'
Med_TCDQ.standard_name = 'Mediterannean and Black Seas moisture flux divergence'
Med_TCDQ[:,:] = mednan

Med_mask = newnc.createVariable('maskM',pl.float64,('lat','lon'))
Med_mask.units = '0 or 1'
Med_mask.standard_name = 'Mediterannean and Black Seas mask'
Med_mask[:,:] = med_fnlmask

newnc.close()