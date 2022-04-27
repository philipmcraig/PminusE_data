# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 15:22:14 2015

@author: np838619

Code for plotting DJF & JJA mean 1000 hPa q and winds from ERA-Interim

Last updated: 08/09/2015 10:48AM 2nd commit
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import os
from mpl_toolkits.basemap import Basemap, shiftgrid
import glob
import scipy.ndimage as scifilt
from matplotlib import gridspec
from functions import *

pl.close()
homedir = '/home/users/qx911590/np838619/'
# tcdq_36years_means.nc
# /home/users/np838619/Watershed/wvfluxes_7914.nc

# set up years array here
#years = pl.linspace(2010,2014,5)
#year_input = [str(int(i)) for i in years]

# set up empty filenames array here
#filenames = pl.zeros([len(year_input),12],dtype='S13')
# loop over years to stick file names into array using PrintFiles
#for year in range(len(year_input)):
#        path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + year_input[year]
#        filenames[year] = PrintFiles(path,'ggaw') #ggap for u &v on levels, ggas for u10 & v10

# set up empty q, u & v arrays here
#tcwv = pl.zeros([len(year_input),len(filenames[0]),1,1,256,512])
#tcuq = pl.zeros_like(tcwv)
#tcvq = pl.zeros_like(tcwv)
## loop over years and filenames and extract data, stick into arrays
#for year in range(len(year_input)):
#    for name in range(len(filenames[year])):
#        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + 
#                str(year_input[year]) + '/' + str(filenames[year,name]),'r')
#        # change the second index in the second square brackets to change pressure level
#        tcwv[year,name] = ncfile.variables['TCWV'][:]
#        tcuq[year,name] = ncfile.variables['TCUQ'][:]
#        tcvq[year,name] = ncfile.variables['TCVQ'][:]
#        ncfile.close()

ncfile = Dataset(homedir+'PminusE_data/ERA_Int/tcdq_36years_means.nc','r')
tcdq = ncfile.variables['tcdq'][:]
lat = ncfile.variables['lat'][:]
lon = ncfile.variables['lon'][:]
ncfile.close()

ncfile = Dataset(homedir+'Watershed/wvfluxes_7914.nc','r')
tcuq = ncfile.variables['tcuq'][:]
tcvq = ncfile.variables['tcvq'][:]
ncfile.close()

tcdq = pl.squeeze(tcdq); tcuq = pl.squeeze(tcuq); tcvq = pl.squeeze(tcvq)
# calculate means from 1980-2014

tcdq_mn = pl.mean(tcdq,axis=0); tcdq_mn = -1*tcdq_mn
tcuq_mn = pl.mean(tcuq,axis=0)
tcvq_mn = pl.mean(tcvq,axis=0)
#

tcdq_mn = scifilt.gaussian_filter(tcdq_mn, sigma=(2.5,1), order=0)
tcdq_mn = tcdq_mn*((86400*365)/1000)
#
#J = pl.mean(tcuq_mn,axis=1)
#K=tcuq_mn - J[:,None]

# plot the means

#latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
#                                    year_input[0] + '/' + filenames[0,0],'r')
#lon = latlon.variables['longitude'][:]
#lat = latlon.variables['latitude'][:]
#latlon.close()

zm = pl.mean(tcdq_mn,axis=1)

#
fig = pl.figure(figsize=(12.5,7.05))
#fig, (ax1, ax2) = pl.subplots(1,2,figsize=(12,7))
gs = gridspec.GridSpec(5,4)
ig = [gs[1:4,0],gs[:,1:]]

ax1 = pl.subplot(ig[0])#pl.subplot(121)
ax1.plot(zm,lat,lw=2,zorder=2,color='k')
ax1.set_yticks([-60,-30,0,30,60])
ax1.grid(axis='y',ls='--')
pl.axvline(x=0,ls='--',color='grey',zorder=0)
pl.ylim(-80,80); pl.xlim(-1.0,1.1)
ax1.tick_params(direction='in'); pl.xticks(fontsize=12); pl.yticks(fontsize=12)
pl.xlabel('$[\overline{p-e}\,]$ (m yr$^{-1}$)',fontsize=12)
pl.ylabel('latitude',fontsize=12,labelpad=-10)
ax1.annotate('(a)',(-0.95,70),fontsize=15)

ax2 = pl.subplot(ig[1])#pl.subplot(122)
m = Basemap(projection='cyl',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
        llcrnrlon=-180.,urcrnrlon=179.5,lat_ts=10)
tcuq_mn, newlons = shiftgrid(180.0, tcuq_mn, lon, start=False)
tcvq_mn, newlons = shiftgrid(180.0, tcvq_mn, lon, start=False)
tcdq_mn, newlons = shiftgrid(180.0, tcdq_mn, lon, start=False)
m.drawcoastlines(linewidth=0.5)
lons, lats = pl.meshgrid(newlons,lat)
X, Y = m(lons,lats)
#
##norm = MidpointNormalize(midpoint=0)
#
cmap = pl.get_cmap('seismic_r')
levels=[-4,-2,-1,-0.5,-0.25,-0.125,0.125,0.25,0.5,1,2,4]
cs = m.contourf(X,Y,tcdq_mn,cmap=cmap,norm=pl.Normalize(-4,4,clip=False),
                levels=levels,extend='max',alpha=0.5)
m.fillcontinents(color='whitesmoke')
#cbar = m.colorbar(cs,location='bottom', pad = "2%",extend="min",
#                  ticks=[-4,-1,-0.25,0,0.25,1,4])
#cbar.set_label('m yr$^{-1}$',fontsize=16)
#cbar.ax.tick_params(labelsize=14); cbar.update_ticks()
#ugrid,newlons2 = shiftgrid(180.,pl.flipud(u_JJA_mn[0,0]),lon,start=False)
#vgrid,newlons2 = shiftgrid(180.,pl.flipud(v_JJA_mn[0,0]),lon,start=False)
#m.streamplot(X,Y,pl.squeeze(ugrid),pl.squeeze(vgrid))

uproj,vproj,xx,yy = m.transform_vector(pl.flipud(tcuq_mn[:]),pl.flipud(tcvq_mn[:]),newlons,lat[::-1],
                                                   nx=50,ny=50,returnxy=True)
Q = m.quiver(xx,yy,uproj,vproj,width=0.001,scale_units='dots',scale=7.5,zorder=5)
pl.matplotlib.pyplot.quiverkey(Q,0.73,0.025,200,'200 kg m$^{-1}$ s$^{-1}$',
                               labelpos='E',fontproperties={'size':16})

m.drawmeridians([-90,0,90],dashes=[1,2],linewidth=0.5,labels=[0,0,1,0],
                yoffset=0.001,fontsize=12)
m.drawparallels([-60,-30,0,30,60],dashes=[1,2],linewidth=0.5,labels=[0,1,0,0],
                xoffset=1.5,fontsize=12)

#m.drawparallels([-35,60],linewidth=0.7,labels=[0,0,0,0],color='grey',dashes=[1,1])
pl.axhline(y=60,color='k',linewidth=1.7)
pl.axhline(y=-35,color='k',linewidth=1.7)

ax2.annotate('(b)',(-169,69),fontsize=15,backgroundcolor='w',
             bbox=dict(facecolor='w',alpha=1,BoxStyle=round))

f = pl.gcf()
colax = f.add_axes([0.32,0.185,0.6,0.03])                   
cb = pl.colorbar(cs,orientation='horizontal',cax=colax,
                 ticks=[-4,-1,-0.25,0,0.25,1,4])
cb.ax.tick_params(labelsize=14); cb.update_ticks()
cb.set_label('m yr$^{-1}$',fontsize=16,labelpad=-2)
#fig.set_size_inches(12,7)
pl.tight_layout()
pl.subplots_adjust(right=0.96,left=0.05)
#m.contour(X,Y,pl.squeeze(u_JJA_mn),colors='b')
#m.contour(X,Y,pl.squeeze(v_JJA_mn),colors='r')