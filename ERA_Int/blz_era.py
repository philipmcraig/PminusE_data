# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 17:12:47 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs
from functions import *

def SeasonalCyc(variable_array):
    """Function to calculate the climatological seasonal means (DJF,MAM,JJA,SON)
    of a variable.
    
    Args:
        variable_array (array): variable of which the climatological seasonal 
                                means are required
    
    Returns:
        variable_seasons (array): climatological seasonal means of the input 
                                  variable
    """
    # Calculate the mean of each trio of months for each year, missing out the 
    # first year as there is no data for December the year before the data starts
    MAM = variable_array[:,2:5] # March, April, May
    JJA = variable_array[:,5:8] # June, July, August
    SON = variable_array[:,8:11] # September, October, November
    DJF = pl.zeros_like(MAM) # December, January, February
    DJF[:,0] = variable_array[:,-1]; DJF[:,1:] = variable_array[:,:2]
    
    # Calculate the climatological mean of each season:  
    MAM_mn = pl.mean(MAM,axis=(0,1)) # only need axis=0 for profiles!!!!
    JJA_mn = pl.mean(JJA,axis=(0,1))
    SON_mn = pl.mean(SON,axis=(0,1))
    DJF_mn = pl.mean(DJF,axis=(0,1))
    
    # Stick all seasons in one array before outputting:
    variable_seasons = pl.array([DJF_mn,MAM_mn,JJA_mn,SON_mn])
    
    return variable_seasons

pl.close('all')
eradir = '/panfs/jasmin/era/era-in/netc/monthly_means/'

nc2 = Dataset('/home/np838619/Watershed/ggis198901010000.nc','r')
eralon = nc2.variables['longitude'][:]
eralat = nc2.variables['latitude'][:]
geop = nc2.variables['Z'][:]
#lsm = nc2.variables['LSM'][:]
nc2.close()
geop = pl.squeeze(geop)#; lsm = pl.squeeze(lsm)

# calculate orography field by dividing geopotential by gravity:
orog = geop/9.8066

years = pl.linspace(2010,2014,5)
year_input = [str(int(i)) for i in years]

# empty array for filenames
filenames = pl.zeros([len(years),96],dtype='S17') # Only need every 4th file!

#loop over years:
for year in range(len(years)):
    #path = path to mm folder + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + str(year_input[year])
    #filenames[year] = PrintFiles(path,type)
    filenames[year] = PrintFiles(path,'hgfs')

lessfiles = pl.zeros([filenames.shape[0],filenames.shape[1]/4],dtype='S17')
for year in range(filenames.shape[0]):
    filelist = []
    for name in range(filenames.shape[1]):
        if '12.nc' in filenames[year,name]:
            filelist.append(filenames[year,name])
    for i in range(len(filelist)):
        lessfiles[year,i] = filelist[i]

#for yr in range(years.size):
#    filenames[yr] = PrintFiles(eradir+str(int(years[yr]))+'/','ggfs')
#    filenames[yr] = pl.sort(filenames[yr])

blz = pl.zeros([filenames.shape[0],filenames.shape[1],1,1,256,512])
#
for yr in range(years.size):
    for mn in range(filenames.shape[1]):
        ncfile = Dataset(eradir+str(int(years[yr]))+'/'+filenames[yr,mn],'r')
        blz[yr,mn] = ncfile.variables['BLH'][:]
        ncfile.close()

blz = pl.squeeze(blz)

x=range(filenames.shape[1]); x = x[::8]
blz_dm = pl.zeros([blz.shape[0],blz.shape[1]/8,eralat.size,eralon.size])
#
for yr in range(years.size):
    for mn in range(len(x)):
        for i in range(eralat.size):
            for j in range(eralon.size):
                blz_dm[yr,mn,i,j] = max(blz[yr,x[mn]:x[mn]+8,i,j])
#
blz_mn = pl.mean(blz_dm,axis=(0,1))
#blz_sns = SeasonalCyc(blz_dm)
##
eralon[-1] = 360
norm = pl.Normalize(0,2500)
levels = pl.linspace(0,2500,11)#[0,200,400,600,800,1000]
ax = pl.subplot(projection=ccrs.PlateCarree(central_longitude=0))
lons,lats = pl.meshgrid(eralon,eralat)
cs = ax.contourf(lons,lats,blz_mn,cmap='inferno_r',transform=ccrs.PlateCarree(),
                 norm=norm,levels=levels,extend='max')
ax.gridlines(); ax.coastlines()
cb = pl.colorbar(cs,orientation='horizontal',pad=0.03)
cb.set_label('metres',fontsize=18,labelpad=-1)
cb.ax.tick_params(labelsize=16)
#
#pl.tight_layout()


#fig, ax = pl.subplots(2,2,figsize=(10,10))
#for i in range(4):
#    axx = pl.subplot(2,2,i+1,projection=ccrs.Robinson(central_longitude=0))
#    cs = axx.contourf(lons,lats,blz_sns[i],cmap='plasma',transform=ccrs.PlateCarree(),
#                      norm=norm,levels=levels,extend='max')
#    axx.gridlines(); axx.coastlines()
#
#pl.tight_layout()