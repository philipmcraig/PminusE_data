# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 14:19:46 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from functions import *
import cartopy
import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib as mp

pl.close('all')

ncfile = Dataset('/home/np838619/PminusE_data/ERA_Int/evapprec_36yrs.nc','r')
lon = ncfile.variables['lon'][:]
lat = ncfile.variables['lat'][:]
ncfile.close()

#---------READ ALL THE ERA-INTERIM NETCDF FILES HERE USING PRINTFILES----------
# need array for years, empty array for file names, tcdq & tcwv
years = pl.linspace(1979,2014,36).astype(int)
year_input = [str(i) for i in years] # convert to string for including in path
filenames = pl.zeros([len(years),12],dtype='S17')
for year in range(len(year_input)):
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + year_input[year]
    filenames[year] = PrintFiles(path,'ggaw')
# need empty arrays for tcdq & tcwv
tcdq = pl.zeros([len(year_input),len(filenames[0]),1,1,256,512])
tcwv = pl.zeros_like(tcdq); tcuq = pl.zeros_like(tcdq); tcvq = pl.zeros_like(tcdq)
# LOOP OVER YEARS AND EXTRACT TCDQ & TCWV FROM NC FILES USING FILENAMES ARRAY
for year in range(len(year_input)):
    for name in range(len(filenames[year])):
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                    str(year_input[year]) + '/' + str(filenames[year,name]),'r')
        tcdq[year,name] = ncfile.variables['TCDQ'][:]
        tcwv[year,name] = ncfile.variables['TCWV'][:]
        tcuq[year,name] = ncfile.variables['TCUQ'][:]
        tcvq[year,name] = ncfile.variables['TCVQ'][:]
        ncfile.close()


tcdq = pl.squeeze(tcdq); tcwv = pl.squeeze(tcwv)
tcuq = pl.squeeze(tcuq); tcvq = pl.squeeze(tcvq)

# Use moisture budget to get monthly E-P arrays
tcwv_mnths_dt = pl.zeros_like(tcwv)
for i in range(tcwv.shape[0]):
    tcwv_mnths_dt[i] = dtcwv_dt(tcwv[i])
EmP = tcwv_mnths_dt + tcdq # E-P = d(TCWV)/dt + TCDQ

# Calculate the seasonal mean E-P:
MAM = pl.mean(EmP[:,2:5,:,:],axis=1)
JJA = pl.mean(EmP[:,5:8,:,:],axis=1)
SON = pl.mean(EmP[:,8:11,:,:],axis=1)
DJF = pl.zeros_like(MAM)
DJF[0,:,:] = pl.float64('nan')
for i in range(1,DJF.shape[0]):
    DJF[i,:,:] = (EmP[i-1,-1,:,:]+EmP[i,0,:,:]+EmP[i,1,:,:])/3

# Combine them all into 1 array to make things easier later:
EmP_sns = pl.array([DJF,MAM,JJA,SON])

JJA_mn = pl.mean(JJA,axis=0)
JJA_an = JJA - JJA_mn

# seasonal means of tcuq and tcvq:
MAM_u = pl.mean(tcuq[:,2:5,:,:],axis=1); MAM_v = pl.mean(tcvq[:,2:5,:,:],axis=1)
JJA_u = pl.mean(tcuq[:,5:8,:,:],axis=1); JJA_v = pl.mean(tcvq[:,5:8,:,:],axis=1)
SON_u = pl.mean(tcuq[:,8:11,:,:],axis=1); SON_v = pl.mean(tcvq[:,8:11,:,:],axis=1)
DJF_u = pl.zeros_like(MAM_u); DJF_v = pl.zeros_like(MAM_v)
DJF_u[0,:,:] = pl.float64('nan'); DJF_v[0,:,:] = pl.float64('nan')
for i in range(1,DJF_u.shape[0]):
    DJF_u[i,:,:] = (tcuq[i-1,-1,:,:]+tcuq[i,0,:,:]+tcuq[i,1,:,:])/3
    DJF_v[i,:,:] = (tcvq[i-1,-1,:,:]+tcvq[i,0,:,:]+tcvq[i,1,:,:])/3

JJAu_mn = pl.mean(JJA_u,axis=0)
JJAu_an = JJA_u - JJAu_mn
JJAv_mn = pl.mean(JJA_v,axis=0)
JJAv_an = JJA_v - JJAv_mn


fig,ax = pl.subplots(2,3)
proj = ccrs.LambertCylindrical()
norm=pl.Normalize(-9,9); cmap = 'seismic'
levels=[-11,-9,-7,-5,-3,-1,1,3,5,7,9]
inds = [16,26,31]

for i in range(3):
    axx = pl.subplot(2,3,i+1,projection=proj)
    axx.set_extent([-120,0,-30,30],crs=ccrs.PlateCarree()); axx.coastlines()
    lons,lats = pl.meshgrid(lon,lat)
    axx.contourf(lons,lats,JJA_an[inds[i]]*86400,norm=norm,cmap=cmap,
                levels=levels,extend='both',transform=proj)
#pl.tight_layout()

#fig,ax = pl.subplots(1,3)
for i in range(3,6):
    axx = pl.subplot(2,3,i+1)
    m = Basemap(projection='cyl',llcrnrlon=-120,llcrnrlat=-30,urcrnrlon=0,
                urcrnrlat=30)
    m.drawcoastlines()
    uproj,vproj,xx,yy = m.transform_vector(pl.flipud(JJAu_an[inds[i-3]]),
                                       pl.flipud(JJAv_an[inds[i-3]]),lon,
                                        lat[::-1],nx=25,ny=25,returnxy=True)
    qp = mp.pyplot.quiver(xx,yy,uproj,vproj,scale=750,width=0.001,headwidth=8,
                          headlength=8)
pl.tight_layout()  