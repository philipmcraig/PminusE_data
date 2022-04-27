# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 13:47:47 2015

@author: np838619

Code to calculate climatological monthly mean evaporation & precipitation
from ERA-Interim for each ocean basin and show time series of the intra-annual
variability.

Last updated: 3/03/2016 5:05PM 7th commit
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from functions import *

pl.close('all')

def BandPlusSeaScaled(areas,net_flux,sea_flux,basin_mask1,basin_mask2,lsm,indices):
    """
    """
    band_areas = areas*basin_mask1*(1-lsm)
    sea_areas = areas*basin_mask2*(1-lsm)
    band_sum = pl.sum(band_areas[indices[0]:indices[1]+1])
    sea_sum = pl.sum(sea_areas)
    scaled_flux = (net_flux+sea_flux)/(band_sum+sea_sum)
    
    return scaled_flux

# read in nc files with basin masks
# extract lat & lon
# also need ERA-Interim mask
basinfile = Dataset('tcdq_basins_60.nc','r')
atlmask = basinfile.variables['maskA'][:] # mask for Atlantic Ocean
balmask = basinfile.variables['maskB'][:] # mask for Baltic Sea
medmask = basinfile.variables['maskM'][:] # mask for Med sea
pacmask = basinfile.variables['maskP'][:] # mask for Pacific Ocean
indmask = basinfile.variables['maskI'][:] # mask for Indian Ocean
lon = basinfile.variables['lon'][:] # ERA-Interim longitude
lat = basinfile.variables['lat'][:] # ERA-Interim latitude
basinfile.close()

maskfile = Dataset('tcdq_ann_mean_ERAI.nc','r')
lsm = maskfile.variables['LSM'][:]
maskfile.close()

# add the Atlantic & Baltic masks together:
ABMmask = atlmask + balmask + medmask

# list of years
years = pl.linspace(1979,2014,36)
# need next part to get rid of decimal point
year_input = [str(int(i)) for i in years]

# empty array for filenames
filenames = pl.zeros([len(years),96],dtype='S17') # Only need every 4th file!

#loop over years:
for year in range(len(years)):
    #path = path to mm folder + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + str(year_input[year])
    #filenames[year] = PrintFiles(path,type)
    filenames[year] = PrintFiles(path,'hafs')

lessfiles = pl.zeros([filenames.shape[0],filenames.shape[1]/4],dtype='S17')
for year in range(filenames.shape[0]):
    filelist = []
    for name in range(filenames.shape[1]):
        if '12.nc' in filenames[year,name]:
            filelist.append(filenames[year,name])
    for i in range(len(filelist)):
        lessfiles[year,i] = filelist[i]

#empty array for evaporation
evap = pl.zeros([len(years),lessfiles.shape[1],1,1,256,512])
#empty array for precipitation
prec = pl.zeros_like(evap)

lessfiles = pl.sort(lessfiles)

#loop over years:
for year in range(len(years)):
    #loop over filenames:
    for name in range(lessfiles.shape[1]):
        #load ncfile
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                str(year_input[year]) + '/' + str(lessfiles[year,name]),'r')
        #extract E & TP data
        evap[year,name] = ncfile.variables['E'][:]
        prec[year,name] = ncfile.variables['TP'][:]
        ncfile.close()
        
# why are there -ve values of TP? Nothing to worry about.

evap_tot = pl.zeros([evap.shape[0],evap.shape[1]/2,evap.shape[2],evap.shape[3],
                     evap.shape[4],evap.shape[5]])
prec_tot = pl.zeros_like(evap_tot)
# loop over number of years:
for year in range(evap.shape[0]):
    # loop over number of months:
    for month in range(int(evap.shape[1]/2)):
        evap_tot[year,month] = -1*(evap[year][2*month] + evap[year][2*month+1])
        prec_tot[year,month] = (prec[year][2*month] + prec[year][2*month+1])


# which latitudes will be required for splitting the basins into subregions?
atl_lats = pl.array([60,35,15,-15,-35])
pac_lats = pl.array([65.5,35,15,-15,-30])
ind_lats = pl.array([32,8,-15,-35])

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
    
    

# Calculate the climatological monthly means for E & P:
evap_mnths = pl.mean(evap_tot,axis=0)
prec_mnths = pl.mean(prec_tot,axis=0)

# get rid of the one-dimensional axes:
evap_mnths = pl.squeeze(evap_mnths)
prec_mnths = pl.squeeze(prec_mnths)

# extract the relevant ocean basins:

ABM_evap = evap_mnths*ABMmask*(1-lsm)
ABM_prec = prec_mnths*ABMmask*(1-lsm)

atl_evap = evap_mnths*atlmask*(1-lsm)
atl_prec = prec_mnths*atlmask*(1-lsm)

bal_evap = evap_mnths*balmask*(1-lsm)
bal_prec = prec_mnths*balmask*(1-lsm)

med_evap = evap_mnths*medmask*(1-lsm)
med_prec = prec_mnths*medmask*(1-lsm)

pac_evap = evap_mnths*pacmask*(1-lsm)
pac_prec = prec_mnths*pacmask*(1-lsm)

ind_evap = evap_mnths*indmask*(1-lsm)
ind_prec = prec_mnths*indmask*(1-lsm)


#---------------------CALCULATE AREAS OF GRID CELLS----------------------------

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


# Use NetEmPCalc but need to include a rescaling factor:
rho = 1000 # density in kg/m^3
factor = rho/86400 # rescaling factor

# set up empty arrays:
atl_evapint = pl.zeros_like(atl_evap[:,0,0])
atl_precint = pl.zeros_like(atl_evapint)
pac_evapint = pl.zeros_like(atl_evapint)
pac_precint = pl.zeros_like(atl_evapint)
ind_evapint = pl.zeros_like(atl_evapint)
ind_precint = pl.zeros_like(atl_evapint)

for m in range(atl_evapint.shape[0]):
    atl_evapint[m] = NetEmPCalc(areas*(1-lsm[0])*ABMmask,ABM_evap[m]*factor,rho)
    atl_precint[m] = NetEmPCalc(areas*(1-lsm[0])*ABMmask,ABM_prec[m]*factor,rho)
    pac_evapint[m] = NetEmPCalc(areas*(1-lsm[0])*pacmask,pac_evap[m]*factor,rho)
    pac_precint[m] = NetEmPCalc(areas*(1-lsm[0])*pacmask,pac_prec[m]*factor,rho)
    ind_evapint[m] = NetEmPCalc(areas*(1-lsm[0])*indmask,ind_evap[m]*factor,rho)
    ind_precint[m] = NetEmPCalc(areas*(1-lsm[0])*indmask,ind_prec[m]*factor,rho)


mnth_lbls = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
fig, ax  = pl.subplots(3,1,figsize=(10,10))

ax1 = pl.subplot(311)
pl.plot(pl.linspace(0,11,12),atl_evapint,color='r')
pl.plot(pl.linspace(0,11,12),atl_precint,color='b')
pl.xticks(pl.linspace(0,11,12))
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.yticks(pl.linspace(2.0,3.6,5))
ax1.set_ylabel('Sv',fontsize=18)
pl.title('Atlantic',fontsize=18)

ax2 = pl.subplot(312)
pl.plot(pl.linspace(0,11,12),pac_evapint,color='r')
pl.plot(pl.linspace(0,11,12),pac_precint,color='b')
pl.xticks(pl.linspace(0,11,12))
ax2.xaxis.set_major_formatter(pl.NullFormatter())
pl.yticks(pl.linspace(5.8,7.0,5))
ax2.set_ylabel('Sv',fontsize=18)
pl.title('Pacific',fontsize=18)

ax3 = pl.subplot(313)
pl.plot(pl.linspace(0,11,12),ind_evapint,color='r',label='$E$')
pl.plot(pl.linspace(0,11,12),ind_precint,color='b',label='$P$')
pl.xticks(pl.linspace(0,11,12))
ax3.set_xticklabels(mnth_lbls)
ax3.set_xlabel('months',fontsize=16)
pl.yticks(pl.linspace(1.3,2.8,4))
ax3.set_ylabel('Sv',fontsize=18)
pl.title('Indian',fontsize=18)
ax3.legend(loc=0)

fig.suptitle('Intra-annual variability of ERA-Interim evaporation\n and precipitation 1979-2014',
             fontsize=19)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/E_TP_intra_basins.png')
             

#--------------------------------SCALE BY AREA---------------------------------

scl_fctr = (10**6)*1000*86400 # factor to rescale from Sv/m^2 to mm/day

# setup empty arrays for scaled fluxes:
atl_evapscl = pl.zeros_like(atl_evapint)
atl_precscl = pl.zeros_like(atl_evapscl)
pac_evapscl = pl.zeros_like(atl_evapscl)
pac_precscl = pl.zeros_like(atl_evapscl)
ind_evapscl = pl.zeros_like(atl_evapscl)
ind_precscl = pl.zeros_like(atl_evapscl)
for m in range(atl_evapscl.shape[0]):
    atl_evapscl[m] = FluxScaled(areas*(1-lsm[0]),ABMmask,atl_evapint[m])*scl_fctr
    atl_precscl[m] = FluxScaled(areas*(1-lsm[0]),ABMmask,atl_precint[m])*scl_fctr
    pac_evapscl[m] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_evapint[m])*scl_fctr
    pac_precscl[m] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_precint[m])*scl_fctr
    ind_evapscl[m] = FluxScaled(areas*(1-lsm[0]),indmask,ind_evapint[m])*scl_fctr
    ind_precscl[m] = FluxScaled(areas*(1-lsm[0]),indmask,ind_precint[m])*scl_fctr

pl.rcParams['xtick.major.pad']='8'    
fig, ax  = pl.subplots(3,1,figsize=(10,10))

ax1 = pl.subplot(311)
pl.plot(pl.linspace(0,11,12),atl_evapscl,color='r',linewidth=3,marker='o',mew=0,ms=8)
pl.plot(pl.linspace(0,11,12),atl_precscl,color='b',linewidth=3,marker='s',mew=0,ms=8)
pl.xticks(pl.linspace(0,11,12)); pl.xlim(0,11)
ax1.xaxis.set_major_formatter(pl.NullFormatter())
#ax1.annotate('(a)',(0,1.02),xycoords='axes fraction',size=25)
pl.yticks(pl.linspace(2,6,5)); pl.grid(axis='y')
pl.setp(ax1.get_yticklabels(),fontsize=18)
ax1.set_ylabel('mm/day',fontsize=20)
pl.title('(a) Atlantic',fontsize=21)

ax2 = pl.subplot(312)
pl.plot(pl.linspace(0,11,12),pac_evapscl,color='r',linewidth=3,label='$\\bar{e}$',
        marker='o',mew=0,ms=8)
pl.plot(pl.linspace(0,11,12),pac_precscl,color='b',linewidth=3,label='$\\bar{p}$',
        marker='s',mew=0,ms=8)
pl.xticks(pl.linspace(0,11,12)); pl.xlim(0,11)
ax2.xaxis.set_major_formatter(pl.NullFormatter())
#ax2.annotate('(b)',(0,1.02),xycoords='axes fraction',size=25)
pl.yticks(pl.linspace(2,6,5)); pl.grid(axis='y')
pl.setp(ax2.get_yticklabels(),fontsize=18)
ax2.set_ylabel('mm/day',fontsize=20)
pl.title('(b) Pacific',fontsize=21)
ax2.legend(loc=4,fontsize=25,ncol=2)

ax3 = pl.subplot(313)
pl.plot(pl.linspace(0,11,12),ind_evapscl,color='r',linewidth=3,marker='o',mew=0,ms=8)
pl.plot(pl.linspace(0,11,12),ind_precscl,color='b',linewidth=3,marker='s',mew=0,ms=8)
pl.xticks(pl.linspace(0,11,12)); pl.xlim(0,11)
ax3.set_xticklabels(mnth_lbls)
pl.setp(ax3.get_xticklabels(),fontsize=18)
ax3.set_xlabel('months',fontsize=18)
#ax3.annotate('(c)',(0,1.02),xycoords='axes fraction',size=25)
pl.yticks(pl.linspace(2,6,5)); pl.grid(axis='y')
pl.setp(ax3.get_yticklabels(),fontsize=18)
ax3.set_ylabel('mm/day',fontsize=20)
pl.title('(c) Indian',fontsize=21)

pl.tight_layout()
#fig.suptitle('Intra-annual variability of ERA-Interim evaporation\n and precipitation 1979-2014 scaled by area',
#             fontsize=19)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/E_TP_intra_basins_scl.png')
             

#----------------CALCULATED FLUXES SCALED BY AREA FOR SUBREGIONS---------------

# Net fluxes in bands:
atl_evp_net = pl.zeros([len(mnth_lbls),atl_ind.shape[0]-1])
atl_prc_net = pl.zeros_like(atl_evp_net)

for m in range(atl_evp_net.shape[0]):
    atl_evp_net[m] = NetEmPBands(atl_ind,areas*(1-lsm[0])*atlmask,atl_evap[m]*factor,rho)
    atl_prc_net[m] = NetEmPBands(atl_ind,areas*(1-lsm[0])*atlmask,atl_prec[m]*factor,rho)

atl_evp_bnds = pl.zeros([len(mnth_lbls),atl_ind.shape[0]-1])
atl_prc_bnds = pl.zeros_like(atl_evp_bnds)

bal_evp_mnths = pl.zeros([12]); bal_prc_mnths = pl.zeros_like(bal_evp_mnths);
med_evp_mnths = pl.zeros([12]); med_prc_mnths = pl.zeros_like(bal_evp_mnths)

for m in range(bal_evp_mnths.shape[0]):
    bal_evp_mnths[m] = NetEmPCalc(areas*balmask*(1-lsm[0]),bal_evap[m]*factor,rho)
    bal_prc_mnths[m] = NetEmPCalc(areas*balmask*(1-lsm[0]),bal_prec[m]*factor,rho)
    med_evp_mnths[m] = NetEmPCalc(areas*medmask*(1-lsm[0]),med_evap[m]*factor,rho)
    med_evp_mnths[m] = NetEmPCalc(areas*medmask*(1-lsm[0]),med_evap[m]*factor,rho)

bal_evp_scl = pl.zeros_like(bal_evp_mnths); bal_prc_scl = pl.zeros_like(bal_evp_mnths);
med_evp_scl = pl.zeros_like(bal_evp_mnths); med_prc_scl = pl.zeros_like(bal_evp_mnths)

for m in range(bal_evp_mnths.shape[0]):
    bal_evp_scl[m] = FluxScaled(areas*(1-lsm[0]),balmask,bal_evp_mnths[m])*scl_fctr
    bal_prc_scl[m] = FluxScaled(areas*(1-lsm[0]),balmask,bal_prc_mnths[m])*scl_fctr
    med_evp_scl[m] = FluxScaled(areas*(1-lsm[0]),balmask,med_evp_mnths[m])*scl_fctr
    med_prc_scl[m] = FluxScaled(areas*(1-lsm[0]),balmask,med_prc_mnths[m])*scl_fctr

for m in range(atl_evp_bnds.shape[0]):
    atl_evp_bnds[m] = FluxBandScaled(areas*(1-lsm[0]),atlmask,atl_ind,
                                                        atl_evp_net[m])*scl_fctr
    atl_evp_bnds[m,0] = BandPlusSeaScaled(areas,atl_evp_net[m,0],bal_evp_mnths[m],
                        atlmask,balmask,lsm[0],(atl_ind[0],atl_ind[1]))*scl_fctr
    atl_evp_bnds[m,2] = BandPlusSeaScaled(areas,atl_evp_net[m,0],med_evp_mnths[m],
                        atlmask,medmask,lsm[0],(atl_ind[1],atl_ind[2]))*scl_fctr
    atl_prc_bnds[m] = FluxBandScaled(areas*(1-lsm[0]),atlmask,atl_ind,
                                                        atl_prc_net[m])*scl_fctr
    atl_prc_bnds[m,0] = BandPlusSeaScaled(areas,atl_prc_net[m,0],bal_prc_mnths[m],
                        atlmask,balmask,lsm[0],(atl_ind[0],atl_ind[1]))*scl_fctr
    atl_prc_bnds[m,2] = BandPlusSeaScaled(areas,atl_prc_net[m,2],med_prc_mnths[m],
                        atlmask,medmask,lsm[0],(atl_ind[1],atl_ind[2]))*scl_fctr

# Recalculate the two northern bands to include the Baltic Sea & Med/Black Seas:
# Baltic Sea:


# Med/Black Seas



clrs = ['k','b','r','g']
lbls = ['35N-60N','15N-35N','15S-15N','35S-15S']
fig, ax = pl.subplots(2,1,figsize=(8,8))    
ax1 = pl.subplot(211)             
for p in range(atl_evp_bnds.shape[1]):
    pl.plot(pl.linspace(0,11,12),atl_evp_bnds[:,p],color=clrs[p],label=lbls[p],
            linewidth=2,marker='o',mew=0,ms=8)
pl.xticks(pl.linspace(0,11,12)); pl.xlim(0,11)
ax1.xaxis.set_major_formatter(pl.NullFormatter())
ax1.annotate('(a) Evaporation',(0.05,0.90),xycoords='axes fraction',size=18)
pl.ylim(0,6); pl.grid(axis='y'); pl.yticks(fontsize=15)
#pl.yticks(pl.arange(40,200,40))
ax1.set_ylabel('mm/day',fontsize=20)
pl.title('Atlantic Ocean annual cycle',fontsize=18)


ax2 = pl.subplot(212)                  
for p in range(atl_prc_bnds.shape[1]):
    pl.plot(pl.linspace(0,11,12),atl_prc_bnds[:,p],color=clrs[p],label=lbls[p],
            linewidth=2,marker='s',mew=0,ms=8)
pl.xticks(pl.linspace(0,11,12)); pl.xlim(0,11)
ax2.set_xticklabels(mnth_lbls)
pl.rcParams['xtick.labelsize'] = 16 
pl.yticks(fontsize=16)#pl.rcParams['ytick.labelsize'] = 18 
ax2.set_xlabel('months',fontsize=18)
pl.ylim(0,6); pl.grid(axis='y')
#pl.yticks(pl.arange(40,200,40))
ax2.set_ylabel('mm/day',fontsize=20)
#pl.title('(b) Precipitation')
ax2.legend(loc=1,ncol=2,columnspacing=0.5,fontsize=15)
ax2.annotate('(b) Precipitation',(0.05,0.90),xycoords='axes fraction',size=18)
#fig.suptitle('Atlantic Ocean annual cycle',fontsize=18,y=0.98)
pl.tight_layout()
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/E_TP_atl_anncyc.png')

###############################################################################
pac_evp_bnds = pl.zeros([len(mnth_lbls),pac_ind.shape[0]-1])
pac_prc_bnds = pl.zeros_like(pac_evp_bnds)

for m in range(pac_evp_bnds.shape[0]):
    pac_evp_bnds[m] = FluxBandScaled(areas*(1-lsm[0]),pacmask,pac_ind,
                    NetEmPBands(pac_ind,areas*(1-lsm[0]),pac_evap[m]*factor,rho))*scl_fctr
    pac_prc_bnds[m] = FluxBandScaled(areas*(1-lsm[0]),pacmask,pac_ind,
                    NetEmPBands(pac_ind,areas*(1-lsm[0]),pac_prec[m]*factor,rho))*scl_fctr


lbls = ['35N-BS','15N-35N','15S-15N','30S-15S']
fig, ax = pl.subplots(2,1,figsize=(8,8))    
ax1 = pl.subplot(211)             
for p in range(pac_evp_bnds.shape[1]):
    pl.plot(pl.linspace(0,11,12),pac_evp_bnds[:,p],color=clrs[p],label=lbls[p],
            linewidth=2,marker='o',mew=0,ms=8)
pl.xticks(pl.linspace(0,11,12)); pl.xlim(0,11)
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(0,6); pl.grid(axis='y')
#pl.yticks(pl.arange(40,200,40))
ax1.set_ylabel('mm/day',fontsize=20)
ax1.annotate('(a) Evaporation',(0.01,0.05),xycoords='axes fraction',size=18)
ax1.legend(loc=4,ncol=1,columnspacing=0.1,fontsize=14)
pl.title('Pacific Ocean annual cycle',fontsize=18)


ax2 = pl.subplot(212)                  
for p in range(pac_prc_bnds.shape[1]):
    pl.plot(pl.linspace(0,11,12),pac_prc_bnds[:,p],color=clrs[p],label=lbls[p],
            linewidth=2,marker='s',mew=0,ms=8)
pl.xticks(pl.linspace(0,11,12)); pl.xlim(0,11)
ax2.set_xticklabels(mnth_lbls)
ax2.set_xlabel('months',fontsize=18)
ax2.annotate('(b) Precipitation',(0.01,0.05),xycoords='axes fraction',size=18)
pl.ylim(0,6); pl.grid(axis='y')
#pl.yticks(pl.arange(40,200,40))
ax2.set_ylabel('mm/day',fontsize=20)
#pl.title('Precipitation')

pl.tight_layout()
#fig.suptitle('Pacific Ocean annual cycle',fontsize=18)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/E_TP_pac_anncyc.png')

###############################################################################
ind_evp_bnds = pl.zeros([len(mnth_lbls),ind_ind.shape[0]-1])
ind_prc_bnds = pl.zeros_like(ind_evp_bnds)

for m in range(ind_evp_bnds.shape[0]):
    ind_evp_bnds[m] = FluxBandScaled(areas*(1-lsm[0]),indmask,ind_ind,
                    NetEmPBands(ind_ind,areas*(1-lsm[0]),ind_evap[m]*factor,rho))*scl_fctr
    ind_prc_bnds[m] = FluxBandScaled(areas*(1-lsm[0]),pacmask,ind_ind,
                    NetEmPBands(ind_ind,areas*(1-lsm[0]),ind_prec[m]*factor,rho))*scl_fctr


lbls = ['>8N','15S-8N','35S-15S']
fig, ax = pl.subplots(2,1,figsize=(8,8))    
ax1 = pl.subplot(211)             
for p in range(ind_evp_bnds.shape[1]):
    pl.plot(pl.linspace(0,11,12),ind_evp_bnds[:,p],color=clrs[p],label=lbls[p])
pl.xticks(pl.linspace(0,11,12))
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(0,6)
#pl.yticks(pl.linspace(90,210,5))
ax1.set_ylabel('mm/day',fontsize=16)
pl.title('Evaporation')
ax1.legend(loc=4,ncol=1,columnspacing=0.5,fontsize=13)

ax2 = pl.subplot(212)                  
for p in range(ind_prc_bnds.shape[1]):
    pl.plot(pl.linspace(0,11,12),ind_prc_bnds[:,p],color=clrs[p],label=lbls[p])
pl.xticks(pl.linspace(0,11,12))
ax2.set_xticklabels(mnth_lbls)
ax2.set_xlabel('months',fontsize=16)
pl.ylim(0,6)
#pl.yticks(pl.arange(0,240,40))
ax2.set_ylabel('mm/day',fontsize=16)
pl.title('Precipitation')

fig.suptitle('Indian Ocean annual cycle',fontsize=18)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/E_TP_ind_anncyc.png')
