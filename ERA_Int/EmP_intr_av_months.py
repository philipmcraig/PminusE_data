# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 09:48:30 2015

@author: np838619

Code to calculate the inter-annual E-P variability of each month from ERA-Interim.
E-P is calculated from monthly mean total column water vapour (TCWV) & vertically
integrated moisture flux divergence (TCDQ). The moisture budget is used to
calculate E-P: d(TCWV)dt + TCDQ = E-P, where the time derivative term is a 
centred difference, calculated in the function dtcwv_dt in functions.py.

Time series for each basin and latitude bands therein are produced and standard
standard deviations are calculated for each month for basins & bands.
"""

# import modules up here: pylab, netCDF4, functions
from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from functions import *


# read netcdf files here: tcdq_basins2.nc, tcdq_basins3.nc & file with ERA_I Mask
basinfile = Dataset('tcdq_basins_60.nc','r')
atlmask = basinfile.variables['maskA'][:] # mask for Atlantic Ocean
balmask = basinfile.variables['maskB'][:] # mask for Baltic Sea
pacmask = basinfile.variables['maskP'][:] # mask for Pacific Ocean
indmask = basinfile.variables['maskI'][:] # mask for Indian Ocean
lon = basinfile.variables['lon'][:] # ERA-Interim longitude
lat = basinfile.variables['lat'][:] # ERA-Interim latitude
basinfile.close()

#IPbasinfile = Dataset('tcdq_basins3.nc','r')
#pacmask = IPbasinfile.variables['maskP'][:] # mask for Pacific Ocean
#indmask = IPbasinfile.variables['maskI'][:] # mask for Indian Ocean
#IPbasinfile.close()

maskfile = Dataset('/panfs/jasmin/era/era-in/netc/ggis/1989/jan1989/ggis198901010000.nc','r')
eramask = maskfile.variables['LSM'][:] # ERA-Interim land-sea mask
maskfile.close()


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
tcwv = pl.zeros_like(tcdq)
# LOOP OVER YEARS AND EXTRACT TCDQ & TCWV FROM NC FILES USING FILENAMES ARRAY
for year in range(len(year_input)):
    for name in range(len(filenames[year])):
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                    str(year_input[year]) + '/' + str(filenames[year,name]),'r')
        tcdq[year,name] = ncfile.variables['TCDQ'][:]
        tcwv[year,name] = ncfile.variables['TCWV'][:]
        ncfile.close()


# get rid of the unneccessary axes in tcdq, tcwv & eramask arrays
eramask = pl.squeeze(eramask,axis=(0,1))
tcdq = pl.squeeze(tcdq,axis=(2,3))
tcwv = pl.squeeze(tcwv,axis=(2,3))


# Use moisture budget to get monthly E-P arrays
tcwv_mnths_dt = pl.zeros_like(tcwv)
for i in range(tcwv.shape[0]):
    tcwv_mnths_dt[i] = dtcwv_dt(tcwv[i])
EmP = tcwv_mnths_dt + tcdq # E-P = d(TCWV)/dt + TCDQ


# specify bounding latitudes & convert to nearest indices
# which latitudes will be required for splitting the basins into subregions?
atl_lats = pl.array([60,45,24,-16,-35])
pac_lats = pl.array([65.5,47,24,-30])
ind_lats = pl.array([32,-8,-20,-35])

# find the indices of the closest latitudes on the ERA-Interim grid:
atl_ind = pl.zeros_like(atl_lats)
for i in range(atl_lats.shape[0]):
    atl_ind[i] = NearestIndex(lat,atl_lats[i])

pac_ind = pl.zeros_like(pac_lats)
for i in range(pac_lats.shape[0]):
    pac_ind[i] = NearestIndex(lat,pac_lats[i])
 
ind_ind = pl.zeros_like(ind_lats)
for i in range(ind_lats.shape[0]):
    ind_ind[i] = NearestIndex(lat,ind_lats[i])


# get the half-grid & calculate grid cell areas
# convert lat & lon array to radians:
lat_rad = pl.radians(lat[:])
lon_rad = pl.radians(lon[:])

# calculate half-grid using HalfGrid function
lat_half = HalfGrid(lat_rad)

nlon = lon_rad.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon

# calculate grid cell areas for full Gaussian grid
# set up empty array for area, size lat_half X lon_half
areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
radius = 6.37*(10**6) # radius of Earth, in metres
for i in range(lat_half.shape[0]-1):
    latpair = (lat_half[i],lat_half[i+1])
    for j in range(lon.shape[0]):
        areas[i,j] = AreaFullGaussianGrid(radius,delta_lambda,latpair)


# add the Atlantic and Baltic masks together:
Amask = atlmask + balmask



#-------------CALCULATE NET E-P FOR EACH MONTH OF EACH YEAR--------------------
rho = 10**3 # density of water, kg/m^3
# set up empty arrays for net E-P, size 36X12:
atlEmP_mnths = pl.zeros([EmP.shape[0],EmP.shape[1]])
pacEmP_mnths = pl.zeros_like(atlEmP_mnths)
indEmP_mnths = pl.zeros_like(atlEmP_mnths)
# loop over each year then each month:
for i in range(EmP.shape[0]):
    #for j in range(EmP.shape[1]):
    atlEmP_mnths[i] = EmPyears(EmP[i],Amask,eramask,areas,rho)
    pacEmP_mnths[i] = EmPyears(EmP[i],pacmask,eramask,areas,rho)
    indEmP_mnths[i] = EmPyears(EmP[i],indmask,eramask,areas,rho)

# Calculate standard deviations for each month for each basin
atl_mnths_sd = pl.std(atlEmP_mnths,axis=0)
pac_mnths_sd = pl.std(pacEmP_mnths,axis=0)
ind_mnths_sd = pl.std(indEmP_mnths,axis=0)


#-------------------------TIME SERIES FOR EACH MONTH---------------------------
months = pl.array(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
clrs = ['black','grey','red','sienna','blue','skyblue','green','lime','yellow',
        'goldenrod','purple','orchid']
# ATLANTIC OCEAN 35S-45N MONTHLY INTER-ANNUAL E-P VARIABILITY
fig, ax1= pl.subplots(1,1) # specifying 1,1 makes every plot show for some reason
for i in range(EmP.shape[1]):
    pl.plot(pl.linspace(0,35,36),atlEmP_mnths[:,i],label=months[i],color=clrs[i])
ax1.set_ylabel('$E-P$  (Sv)',fontsize=16)
ax1.set_xticklabels(pl.linspace(1979,2014,8).astype(int)) # display xtick labels as year integers
ax1.set_xlabel('Years',fontsize=16)
pl.ylim(0.6,1.9) # extend the y limit to fit the legend
pl.title('Atlantic inter-annual variability of each month')
pl.legend(loc=2,ncol=6,columnspacing=0.4)
pl.show()

# PACIFIC OCEAN 30S-47N MONTHLY INTER-ANNUAL E-P VARIABILITY
fig, ax2 = pl.subplots(1,1)
for i in range(EmP.shape[1]):
    pl.plot(pl.linspace(0,35,36),pacEmP_mnths[:,i],label=months[i],color=clrs[i])
pl.axhline(linewidth=1,color='k',ls='--') # add a black dashed line for E-P=0
ax2.set_ylabel('$E-P$  (Sv)',fontsize=16)
ax2.set_xticklabels(pl.linspace(1979,2014,8).astype(int))
ax2.set_xlabel('Years',fontsize=16)
pl.ylim(-0.6,1.2)
pl.title('Pacific inter-annual variability of each month')
ax2.legend(loc=4,ncol=6,columnspacing=0.4)
pl.show()

# INDIAN OCEAN >35S MONTHLY INTER-ANNUAL E-P VARIABILITY
fig, ax3 = pl.subplots(1,1)
for i in range(EmP.shape[1]):
    pl.plot(pl.linspace(0,35,36),indEmP_mnths[:,i],label=months[i],color=clrs[i])
pl.axhline(linewidth=1,color='k',ls='--')
ax3.set_ylabel('$E-P$  (Sv)',fontsize=16)
ax3.set_xticklabels(pl.linspace(1979,2014,8).astype(int))
ax3.set_xlabel('Years',fontsize=16)
pl.ylim(-0.4,1.2)
pl.title('Indian inter-annual $E-P$  variability of each month')
ax3.legend(loc=3,ncol=6,columnspacing=0.4)
pl.show()


#----------------------------LATITUDE BANDS------------------------------------
# need empty arrays 3X36X12 (2 for Pacific)
atl_mnths_bnds = pl.zeros([EmP.shape[0],EmP.shape[1],atl_ind.shape[0]-1])
pac_mnths_bnds = pl.zeros([EmP.shape[0],EmP.shape[1],pac_ind.shape[0]-1])
ind_mnths_bnds = pl.zeros([EmP.shape[0],EmP.shape[1],ind_ind.shape[0]-1])

# multiply EmP array by basins masks & 1-ERA-Interim mask:
atl_EmP = EmP*(1-eramask)*atlmask
pac_EmP = EmP*(1-eramask)*pacmask
ind_EmP = EmP*(1-eramask)*indmask

# loop over years and months:
for i in range(EmP.shape[0]):
    for j in range(EmP.shape[1]):
        atl_mnths_bnds[i,j] = NetEmPBands(atl_ind,areas,atl_EmP[i,j],rho)
        pac_mnths_bnds[i,j] = NetEmPBands(pac_ind,areas,pac_EmP[i,j],rho)
        ind_mnths_bnds[i,j] = NetEmPBands(ind_ind,areas,ind_EmP[i,j],rho)

# Calculate standard deviations for each band for each month (12X3):
atl_bnds_sd = pl.std(atl_mnths_bnds,axis=0)
pac_bnds_sd = pl.std(pac_mnths_bnds,axis=0)
ind_bnds_sd = pl.std(ind_mnths_bnds,axis=0)


#-------------------TIME SERIES FOR EACH BAND FOR EACH MONTH-------------------
# list of latitude band boundaries for subplot titles:
abands_titles = ['45$^\circ$n-60$^\circ$N','24$^\circ$N-45$^\circ$N',
                 '16$^\circ$S-24$^\circ$N','35$^\circ$S-16$^\circ$S']
# ATLANTIC OCEAN (24N-45N,16S-24N,35S-16S) MONTHLY INTER-ANNUAL E-P VARIABILITY
pl.figure(4,figsize=(10,10))
for p in range(atl_ind.shape[0]-1): # loop over the number of bands
    ax = pl.subplot(atl_ind.shape[0]-1,1,p+1) # nrows=no. bands, ncols=1, plot_no.=atl_ind index +1
    # Python uses 0 indexing but subplot starts at 1!
    for i in range(EmP.shape[1]): #  loop over months
        pl.plot(pl.linspace(0,35,36),atl_mnths_bnds[:,i,p],label=months[i],color=clrs[i])
    pl.axhline(linewidth=1,color='k',ls='--')
    # use 3 ticks and set the upper & lower ticks to the max/min values for each band:
    ax.set_yticks(pl.linspace(atl_mnths_bnds[:,:,p].min(),atl_mnths_bnds[:,:,p].max(),3))
    ax.set_ylabel('$E-P$  (Sv)',fontsize=16,labelpad=1.75)
    ax.xaxis.set_major_formatter(pl.NullFormatter()) # remove default numbering from x-axes
    pl.title('Atlantic (' + abands_titles[p] + ') inter-annual variability of each month')
pl.subplots_adjust(hspace=0.3,left=0.15) # more space between subplots, resize a bit to fit y-label
ax.set_xticklabels(pl.linspace(1979,2014,8).astype(int)) # label the x-axis of the bottom subplot
ax.set_xlabel('Years',fontsize=16)
ax.legend(loc=(0,-0.85),ncol=6,columnspacing=2.75,fontsize=14,labelspacing=1,
          handletextpad=0.3,borderaxespad=0.1,mode="expand") # mess around with legend to stick it at side
pl.subplots_adjust(top=0.93,bottom=0.16)

pbands_titles = ['47$^\circ$N- Bering Strait','24$^\circ$N-47$^\circ$N','30$^\circ$S-24$^\circ$N']
pl.figure(5,figsize=(10,10))
for p in range(pac_ind.shape[0]-1):
    ax = pl.subplot(pac_ind.shape[0]-1,1,p+1)
    for i in range(EmP.shape[1]):
        pl.plot(pl.linspace(0,35,36),pac_mnths_bnds[:,i,p],label=months[i],color=clrs[i])
    pl.axhline(linewidth=1,color='k',ls='--')
    ax.set_yticks(pl.linspace(pac_mnths_bnds[:,:,p].min(),pac_mnths_bnds[:,:,p].max(),3))
    #pl.matplotlib.ticker.FormatStrFormatter('%.1g')
    ax.set_ylabel('$E-P$  (Sv)',fontsize=16,labelpad=1.75)
    ax.xaxis.set_major_formatter(pl.NullFormatter())
    pl.title('Pacifc (' + pbands_titles[p] + ') inter-annual variability of each month')
pl.subplots_adjust(hspace=0.3,left=0.15)
ax.set_xticklabels(pl.linspace(1979,2014,8).astype(int))
ax.set_xlabel('Years',fontsize=16)
ax.legend(loc=(0,-0.6),ncol=6,columnspacing=2.75,fontsize=14,labelspacing=1,
          handletextpad=0.3,borderaxespad=0.1,mode="expand")
pl.subplots_adjust(top=0.93,bottom=0.16)

ibands_titles = ['>8$^\circ$S','20$^\circ$S-8$^\circ$S','35$^\circ$S-20$^\circ$S']
pl.figure(6)
for p in range(ind_ind.shape[0]-1):
    ax = pl.subplot(ind_ind.shape[0]-1,1,p+1)
    for i in range(EmP.shape[1]):
        pl.plot(pl.linspace(0,35,36),ind_mnths_bnds[:,i,p],label=months[i],color=clrs[i])
    pl.axhline(linewidth=1,color='k',ls='--')
    ax.set_yticks(pl.linspace(ind_mnths_bnds[:,:,p].min(),ind_mnths_bnds[:,:,p].max(),3))
    #pl.matplotlib.ticker.FormatStrFormatter('%.1g')
    ax.set_ylabel('$E-P$  (Sv)',fontsize=16)
    ax.xaxis.set_major_formatter(pl.NullFormatter())
    pl.title('Indian (' + ibands_titles[p] + ') inter-annual variability of each month')
pl.subplots_adjust(hspace=0.3,left=0.15)
ax.set_xticklabels(pl.linspace(1979,2014,8).astype(int))
ax.set_xlabel('Years',fontsize=16)
ax.legend(loc=(1.01,0.05),ncol=1,columnspacing=0.5,fontsize=10,labelspacing=2,
          handletextpad=0.3,borderaxespad=0.1)
pl.show()


#----------------------WRITE STANDARD DEVIATIONS TO A FILE---------------------

"""f = open('inter-ann_mm_stdev60.txt','w')
f.write('Standard deviations of monthly means (inter-annual variability) of E-P' +
' from ERA-Interim for the Atlantic (35S-60N), Pacific (30S-BS) and Indian (>35S)' +
' Oceans. E-P is from the moisture budget. Also presented are standard deviations' +
' for various latitude bands within each basin.\n\nAtlantic E-P (35S-60N) monthly' +
' standard deviations:\n')
months.tofile(f,sep="   ")
f.write('\n')
atl_mnths_sd.tofile(f,sep="   ",format="%.3f")
f.write('\n\nPacific E-P (30S-BS) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
pac_mnths_sd.tofile(f,sep="   ",format="%.3f")
f.write('\n\nIndian E-P (>35S) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
ind_mnths_sd.tofile(f,sep="   ",format="%.3f")
f.write('\n\nAtlantic E-P (24N-60N) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
AM_bnds_sd[:,0].tofile(f,sep="   ",format="%.3f")
f.write('\n\nAtlantic E-P (16S-24N) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
AM_bnds_sd[:,1].tofile(f,sep="   ",format="%.3f")
f.write('\n\nAtlantic E-P (35S-16S) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
AM_bnds_sd[:,2].tofile(f,sep="   ",format="%.3f")
f.write('\n\nPacific E-P (24N-47N) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
pac_bnds_sd[:,0].tofile(f,sep="   ",format="%.3f")
f.write('\n\nPacific E-P (30S-24N) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
pac_bnds_sd[:,1].tofile(f,sep="   ",format="%.3f")
f.write('\n\nIndian E-P (>8S) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
ind_bnds_sd[:,0].tofile(f,sep="   ",format="%.3f")
f.write('\n\nIndian E-P (20S-8S) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
ind_bnds_sd[:,1].tofile(f,sep="   ",format="%.3f")
f.write('\n\nIndian E-P (35S-20S) monthly standard deviations: \n')
months.tofile(f,sep="   ")
f.write('\n')
ind_bnds_sd[:,2].tofile(f,sep="   ",format="%.3f")
f.close()"""