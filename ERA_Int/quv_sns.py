# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 15:22:14 2015

@author: np838619

Code for plotting DJF & JJA mean 1000 hPa q and winds from ERA-Interim

Last updated: 08/09/2015 10:48AM 2nd commit
"""

from __future__ import division
import pylab as pl
import matplotlib as mp
from netCDF4 import Dataset
import os
from mpl_toolkits.basemap import Basemap, shiftgrid
import scipy.ndimage as scifilt
from functions import *

pl.close('all')
# set up years array here
years = pl.linspace(1979,2014,36)
year_input = [str(int(i)) for i in years]

# set up empty filenames array here
filenames = pl.zeros([len(year_input),12],dtype='S13')
# loop over years to stick file names into array using PrintFiles
for year in range(len(year_input)):
        path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + year_input[year]
        filenames[year] = PrintFiles(path,'ggaw') #ggap for u &v on levels, ggas for u10 & v10

# set up empty q, u & v arrays here
tcwv = pl.zeros([len(year_input),len(filenames[0]),1,1,256,512])
tcuq = pl.zeros_like(tcwv); tcvq = pl.zeros_like(tcwv)
tcdq = pl.zeros_like(tcwv)
# loop over years and filenames and extract data, stick into arrays
for year in range(len(year_input)):
    for name in range(len(filenames[year])):
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + 
                str(year_input[year]) + '/' + str(filenames[year,name]),'r')
        # change the second index in the second square brackets to change pressure level
        tcwv[year,name] = ncfile.variables['TCWV'][:]
        tcuq[year,name] = ncfile.variables['TCUQ'][:]
        tcvq[year,name] = ncfile.variables['TCVQ'][:]
        tcdq[year,name] = ncfile.variables['TCDQ'][:]
        ncfile.close()

tcwv = pl.squeeze(tcwv); tcuq = pl.squeeze(tcuq); tcvq = pl.squeeze(tcvq)
tcdq = pl.squeeze(tcdq)

tcwv_sns = SeasonalMeans(tcwv); tcuq_sns = SeasonalMeans(tcuq)
tcvq_sns = SeasonalMeans(tcvq); tcdq_sns = SeasonalMeans(tcdq)

tcdq_sns = scifilt.gaussian_filter(tcdq_sns,sigma=(0,2.5,1),order=0)
#tcdq_sns = tcdq_sns*(86400*365)/1000

#ws = pl.sqrt(u_DJF_mn**2+v_DJF_mn**2)

nc1 = Dataset('/home/np838619/Watershed/ggis198901010000.nc','r')
lon = nc1.variables['longitude'][:]; lat = nc1.variables['latitude'][:]
lsm = nc1.variables['LSM'][:]
nc1.close()

seasons = ['(a) DJF','(b) MAM','(c) JJA','(d) SON']
locs = (5,-80)

# plot the means
fig, ax = pl.subplots(2,2,figsize=(12,10))
lons, lats = pl.meshgrid(lon,lat)
norm = pl.Normalize(0,50); Q = []

for i in range(tcwv_sns.shape[0]):
    axx = pl.subplot(2,2,i+1)

    m = Basemap(projection='cyl',lon_0=180.)#,llcrnrlat=-30,llcrnrlon=200,
#                                            urcrnrlat=60,urcrnrlon=340)
    m.drawcoastlines()
    if i in (0,2):
        m.drawparallels([-60,-30,0,30,60],labels=[1,0,0,1],ax=None,fontsize=13,linewidth=0)
        #axx.set_yticklabels([-45,-30,-15,0,15,30,45])
    else:
        m.drawparallels([-60,-30,0,30,60],labels=[0,1,1,0],ax=None,fontsize=13,linewidth=0)
    if i in (2,3):
        m.drawmeridians([0,90,180,270],labels=[0,0,0,1],ax=None,
                        fontsize=13,yoffset=2.5,linewidth=0)
    else:
        m.drawmeridians([0,90,180,270],labels=[0,0,1,0],ax=None,fontsize=13,
                        linewidth=0,yoffset=0.5)
    X, Y = m(lons,lats)
    
    cq = m.pcolormesh(X,Y,tcwv_sns[i],cmap='GnBu',norm=norm)
    #pl.colorbar(cq,orientation='horizontal',pad=0.01)
    uproj,vproj,xx,yy = m.transform_vector(pl.flipud(tcuq_sns[i]),
                                       pl.flipud(tcvq_sns[i]),lon,
                                        lat[::-1],nx=40,ny=40,returnxy=True)
    qp = mp.pyplot.quiver(xx,yy,uproj,vproj,scale=12000,width=0.001,headwidth=8,
                          headlength=9)
    Q.append(qp)
    #axx.annotate(seasons[i],locs[i],xycoords='figure fraction',fontsize=20)
    axx.text(locs[0],locs[1],seasons[i],bbox={'facecolor':'white'},size=20)

pl.subplots_adjust(hspace=-0.38,wspace=0.12,left=0.05,right=0.97,top=1,bottom=0.04)

mp.pyplot.quiverkey(Q[2],1.8,-0.12,300,'300 kg/m/s',labelpos='E',linewidth=1.5,
                     fontproperties={'size':16})

f = pl.gcf()
colax = f.add_axes([0.11,0.11,0.76,0.05])
clb = pl.colorbar(cq, cax=colax,orientation='horizontal',extend='max')
clb.set_label('kg m$^{-2}$',fontsize=18)
clb.ax.tick_params(labelsize=14)
pl.tight_layout()
pl.subplots_adjust(left=0.05,right=0.95,wspace=0.04)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/Amr_tcwv_sns1014.png')
###############################################################################

fig, ax = pl.subplots(2,2,figsize=(12,10))
norm = pl.Normalize(-0.08,0.08); Q = []
levels = [-0.08,-0.04,-0.02,-0.01,-0.005,0.005,0.01,0.02,0.04,0.08]
locs = (222,-22)#[(119,-25),(0.85,0.60),(0.38,0.23),(0.85,0.23)]

for i in range(tcdq_sns.shape[0]):
    axx = pl.subplot(2,2,i+1)

    m = Basemap(projection='cyl',lon_0=290.,llcrnrlat=-30,llcrnrlon=200,
                                            urcrnrlat=60,urcrnrlon=340)
    m.drawcoastlines()
    if i in (0,2):
        m.drawparallels([-30,-15,0,15,30,45,60],labels=[1,0,0,1],ax=None,fontsize=13)
        #axx.set_yticklabels([-45,-30,-15,0,15,30,45])
    else:
        m.drawparallels([-30,-15,0,15,30,45,60],labels=[0,0,0,0],ax=None)
    if i in (2,3):
        m.drawmeridians(pl.arange(200,350,20),labels=[0,0,0,1],ax=None,
                        fontsize=13,yoffset=2.5)
        #m.drawmeridians([0],labels=[0,0,0,1],ax=None,fotnsize=13,yoffset=2.5)
    else:
        m.drawmeridians(pl.arange(200,350,20),labels=[0,0,0,0,0],ax=None)
    X, Y = m(lons,lats)
    
    cd = m.contourf(X,Y,tcdq_sns[i]*(10**3),cmap='seismic',norm=norm,levels=levels,
                    extend='both')
    uproj,vproj,xx,yy = m.transform_vector(pl.flipud(tcuq_sns[i]),
                                       pl.flipud(tcvq_sns[i]),lon,
                                        lat[::-1],nx=40,ny=40,returnxy=True)
    qp = mp.pyplot.quiver(xx,yy,uproj,vproj,scale=15000,width=0.001,headwidth=8,
                          headlength=8)
    Q.append(qp)
    #axx.annotate(seasons[i],locs[i],xycoords='figure fraction',fontsize=20)
    axx.text(locs[0],locs[1],seasons[i],bbox={'facecolor':'white'},size=20)

pl.subplots_adjust(hspace=-0.38,wspace=0.12,left=0.05,right=0.97,top=1,bottom=0.04)

mp.pyplot.quiverkey(Q[2],0.1,2.15,300,'300 kg/m/s',labelpos='E',linewidth=1.5,
                    fontproperties={'size':16})
f = pl.gcf()
colax = f.add_axes([0.11,0.08,0.76,0.05])
clb = pl.colorbar(cd, cax=colax,orientation='horizontal',extend='max')
clb.set_label('10$^3 $kg m$^{-2}$ s$^{-1}$',fontsize=18)
clb.ax.tick_params(labelsize=14)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/Amr_tcdq_sns1014.png')

f=(tcwv[-1,-1]-tcwv[0,0])/(5*365*86400)
pl.figure(3)
m = Basemap(projection='cyl',lon_0=180.)
m.drawcoastlines()
levels = pl.linspace(-1,1,6)
cm=m.contourf(X,Y,f*(10**7),norm=pl.Normalize(-1,1),cmap='PiYG',extend='both',
              levels=levels)
pl.colorbar(cm,orientation='horizontal',pad=0.05)

masks = Dataset('/home/np838619/PminusE_data/ERA_Int/tcdq_basins_60.nc')
pacmask = masks.variables['maskP'][:]
indmask = masks.variables['maskI'][:]
atlmask = masks.variables['maskA'][:]
masks.close()

lsm = pl.squeeze(lsm)

stor_pac = f*pacmask
stor_ind = f*indmask
stor_atl = f*atlmask

# Convert lat & lon arrays to radians
lat_rad = pl.radians(lat[:])
lon_rad = pl.radians(lon[:])

#--------------------CREATE LATITUDE HALF-GRID ARRAY HERE----------------------
#----------------------see HalfGrid function for details-----------------------
lat_half = HalfGrid(lat_rad)

#--------------------------CALCULATE DELTA LAMBDA HERE-------------------------
# delta_lambda = 2*pi/nlon
nlon = lon_rad.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon

#--------------calculate cell areas here, use function from above--------------
# set up empty array for area, size lat_half X lon_half
areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
radius = 6.37*(10**6)
# loop over latitude and longitude
for i in range(lat_half.shape[0]-1): # loops over 256
    latpair = (lat_half[i],lat_half[i+1])
    for j in range(lon.shape[0]): # loops over 512
        areas[i,j] = AreaFullGaussianGrid(radius,delta_lambda,latpair)

stor_pac = pl.sum(stor_pac*areas); print stor_pac/(10**9)
stor_ind = pl.sum(stor_ind*areas); print stor_ind/(10**9)
stor_atl = pl.sum(stor_atl*areas); print stor_atl/(10**9)