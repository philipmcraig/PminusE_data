# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 15:19:39 2015

Code to calculate net E-P(-R) for 1- degree latitude bands for each basin and
plot as bar charts.

Last updated: 2/03/2016 5:57PM 5th commit

@author: np838619
"""

# import stuff up here
from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import pandas
from functions import *

def BandPlusSeaScaled(areas,net_flux,sea_flux,basin_mask1,basin_mask2,lsm,indices):
    """
    """
    band_areas = areas*basin_mask1*(1-lsm)
    sea_areas = areas*basin_mask2*(1-lsm)
    band_sum = pl.sum(band_areas[indices[0]:indices[1]+1])
    sea_sum = pl.sum(sea_areas)
    scaled_flux = (net_flux+sea_flux)/(band_sum+sea_sum)
    
    return scaled_flux

pl.close('all')

# read in data here
# basin masks are in netcdf files
maskfile = Dataset('tcdq_basins_60.nc','r')
Amask = maskfile.variables['maskA'][:] # atlantic basin mask
Bmask = maskfile.variables['maskB'][:] # baltic basin mask
Mmask = maskfile.variables['maskM'][:] # med sea mask
Pmask = maskfile.variables['maskP'][:] # pacific basin mask
Imask = maskfile.variables['maskI'][:] # indian basin mask
maskfile.close()

#atlmedfile = Dataset('tcdq_basin2.nc','r')
#Amask = atlmedfile.variables['maskA'][:]
#Mmask = atlmedfile.variables['maskMB'][:]
#atlmedfile.close()

# E, P & E-P should be in a netcdf file, I've not done that yet though
fluxfile = Dataset('fluxes_ann_mean.nc','r')
EmP = fluxfile.variables['E-P'][:] # E-P in m/day
evap = fluxfile.variables['E'][:] # E in m/day
prec = fluxfile.variables['P'][:] # P in m/day
lat = fluxfile.variables['lat'][:] # ERA-Interim latitudes
lon = fluxfile.variables['lon'][:] # ERA-Interim longitudes
fluxfile.close()

maskfile = Dataset('tcdq_ann_mean_ERAI.nc','r')
lsm = maskfile.variables['LSM'][:] # ERA-Interim land-sea mask
maskfile.close()

# runoff data:
runoff = pandas.read_csv('DT03_runoff.csv',sep=',',header=1,skiprows=0,
    names=['row','lat','Atl','Pac','Ind','Glo'])
# to get the values of a column type runoff.name e.g. runoff.Pac
run_atl = runoff.Atl
run_pac = runoff.Pac
run_ind = runoff.Ind

# Make a new mask for Atlantic, Mediterranean & Baltic:
#ABmask = Amask + Bmask + Mmask

# Multiply each flux array by each mask & lsm to extract correct flux:
atl_EmP = EmP*Amask*(1-lsm); atl_evap = evap*Amask*(1-lsm); atl_prec = prec*Amask*(1-lsm)
bal_EmP = EmP*Bmask*(1-lsm); bal_evap = evap*Bmask*(1-lsm); bal_prec = prec*Bmask*(1-lsm)
med_EmP = EmP*Mmask*(1-lsm); med_evap = evap*Mmask*(1-lsm); med_prec = prec*Mmask*(1-lsm)
pac_EmP = EmP*Pmask*(1-lsm); pac_evap = evap*Pmask*(1-lsm); pac_prec = prec*Pmask*(1-lsm)
ind_EmP = EmP*Imask*(1-lsm); ind_evap = evap*Imask*(1-lsm); ind_prec = prec*Imask*(1-lsm)

# define arrays for 5 degree latitude bands for each basin here
bands = pl.linspace(60,-30,10)
# need to convert these to nearest index on ERA-Interim latitude grid:
bnd_ind = pl.zeros_like(bands)
for i in range(bands.shape[0]):
    bnd_ind[i] = NearestIndex(lat,bands[i])

#-----------------------CALCULATE AREAS OF GRID CELLS--------------------------

# Convert lat & lon arrays to radians
lat_rad = pl.radians(lat[:])
lon_rad = pl.radians(lon[:])

#--------------------CREATE LATITUDE HALF-GRID ARRAY HERE----------------------
#----------------------see HalfGrid function for details-----------------------

lat_half = HalfGrid(lat_rad)


#------------------------CALCULATE DELTA PHI ARRAY HERE------------------------
# set up empty array for delta phi, same size as lat
# loop over delta phi:
    # delta_phi = lat_half[i] - lat_half[i+1]


#--------------------------CALCULATE DELTA LAMBDA HERE-------------------------
# delta_lambda = 2*pi/nlon
nlon = lon_rad.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon


#---------calculate cell areas here, use AreaFullGaussianGrid function---------
# set up empty array for area, size lat_half X lon_half
areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
radius = 6.37*(10**6)
# loop over latitude and longitude
for i in range(lat_half.shape[0]-1): # loops over 256
    latpair = (lat_half[i],lat_half[i+1])
    for j in range(lon.shape[0]): # loops over 512
        #lonpair = (lon_half[i],lon_half[i+1])
        areas[i,j] = AreaFullGaussianGrid(radius,delta_lambda,latpair)


# Net E-P functions expect a flux in units kg/m^2s but these units are m/day,
# so they need to be rescaled in order to produce helpful results
rho = 1000 # density in kg/m^3
factor = rho/86400 # density/no. of seconds in 1 day


# First calculate the NET fluxes for the Mediterranean & Blatic Seas:
# Mediterranean & Black Seas:
med_EmPnet = NetEmPCalc(areas,med_EmP*factor,rho)
med_evapnet = NetEmPCalc(areas,med_evap*factor,rho)
med_precnet = NetEmPCalc(areas,med_prec*factor,rho)
# Baltic Sea:
bal_EmPnet = NetEmPCalc(areas,bal_EmP*factor,rho)
bal_evapnet = NetEmPCalc(areas,bal_evap*factor,rho)
bal_precnet = NetEmPCalc(areas,bal_prec*factor,rho)

# Calculate runoff in Med & Baltic Seas:
    # Difference between global runoff & sum of Atl, Pac & Ind runoff
API = run_atl + run_pac + run_ind
rundiff = runoff.Glo - API
    # Find indices of max/min latitudes of Med & Baltic Seas:
med_nth = NearestIndex(runoff.lat,47.3); med_sth = NearestIndex(runoff.lat,30.5)
bal_nth = NearestIndex(runoff.lat,66.0); bal_sth = NearestIndex(runoff.lat,53.8)
    # calculate runoff
med_run = pl.sum(rundiff[med_nth:med_sth+1])
bal_run = pl.sum(rundiff[bal_nth:bal_sth+1])

# Calculate E-P-R for Med and Baltic Seas:
med_EPRnet = med_EmPnet - med_run
bal_EPRnet = bal_EmPnet - bal_run

# Use NetEmPBands functions and perhaps EPRBands & TotalRunoff to calculate net
# fluxes in each band:

# Atlantic Ocean:
atl_EmPbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(atl_EmP)*factor,rho)
atl_evapbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(atl_evap)*factor,rho)
atl_precbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(atl_prec)*factor,rho)

atl_run10deg = pl.zeros_like(atl_EmPbnds)
for i in range(len(bnd_ind)-1):
    latpair = (bnd_ind[i],bnd_ind[i+1])
    atl_run10deg[i] = TotalRunoff(runoff.lat,lat,latpair,run_atl)


atl_EPRbnds = EPRbands(atl_EmPbnds,bnd_ind,runoff.lat,lat,run_atl)

# Baltic Sea enters the North Sea in 50N-60N band, add bal_EPRnet to this band:
"""atl_EPRbnds[0] = atl_EPRbnds[0] + bal_EPRnet
atl_evapbnds[0] = atl_evapbnds[0] + bal_evapnet
atl_precbnds[0] = atl_precbnds[0] + bal_precnet

# Mediterranean enters Atlantic in 30N-40N band, add med_EPRnet to this band:
atl_EPRbnds[2] = atl_EPRbnds[2] + med_EPRnet
atl_evapbnds[2] = atl_evapbnds[2] + med_evapnet
atl_precbnds[2] = atl_precbnds[2] + med_precnet"""

# Pacific Ocean:
pac_EmPbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(pac_EmP)*factor,rho)
pac_evapbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(pac_evap)*factor,rho)
pac_precbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(pac_prec)*factor,rho)

pac_run10deg = pl.zeros_like(pac_EmPbnds)
for i in range(len(bnd_ind)-1):
    latpair = (bnd_ind[i],bnd_ind[i+1])
    pac_run10deg[i] = TotalRunoff(runoff.lat,lat,latpair,run_pac)

pac_EPRbnds = EPRbands(pac_EmPbnds,bnd_ind,runoff.lat,lat,run_pac)



# Indian Ocean:
ind_EmPbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(ind_EmP)*factor,rho)
ind_evapbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(ind_evap)*factor,rho)
ind_precbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(ind_prec)*factor,rho)

ind_run10deg = pl.zeros_like(ind_EmPbnds)
for i in range(len(bnd_ind)-1):
    latpair = (bnd_ind[i],bnd_ind[i+1])
    ind_run10deg[i] = TotalRunoff(runoff.lat,lat,latpair,run_ind)

ind_EPRbnds = EPRbands(ind_EmPbnds,bnd_ind,runoff.lat,lat,run_ind)




#---------------------------SCALE EVERYTHING BY AREA---------------------------

factor = (10**6)*100*86400*365 # factor to rescale from Sv/m^2 to cm/yr

# Atlantic Ocean

atl_EPRscl = FluxBandScaled(areas,Amask*(1-lsm[0]),bnd_ind,atl_EPRbnds)*factor
# Recalculate 50N-60N band to include Baltic Sea:
#atl_EPRscl[0] = BandPlusSeaScaled(areas,atl_EPRbnds[0],bal_EPRnet,Amask,Bmask,
#                                        lsm[0],(bnd_ind[0],bnd_ind[1]))*factor
# Recalculate 30N-40N band to include Med Sea:
#atl_EPRscl[2] = BandPlusSeaScaled(areas,atl_EPRbnds[2],med_EPRnet,Amask,Mmask,
#                                        lsm[0],(bnd_ind[2],bnd_ind[3]))*factor

atl_evapscl = FluxBandScaled(areas,Amask*(1-lsm[0]),bnd_ind,atl_evapbnds)*factor
# Recalculate 50N-60N band to include Baltic Sea:
#atl_evapscl[0] = BandPlusSeaScaled(areas,atl_evapbnds[0],bal_evapnet,Amask,Bmask,
#                                        lsm[0],(bnd_ind[0],bnd_ind[1]))*factor
# Recalculate 30N-40N band to include Med Sea:
#atl_evapscl[2] = BandPlusSeaScaled(areas,atl_evapbnds[2],med_evapnet,Amask,Mmask,
#                                       lsm[0],(bnd_ind[2],bnd_ind[3]))*factor

atl_precscl = FluxBandScaled(areas,Amask*(1-lsm[0]),bnd_ind,atl_precbnds)*factor
# Recalculate 50N-60N band to include Baltic Sea:
#atl_precscl[0] = BandPlusSeaScaled(areas,atl_precbnds[0],bal_precnet,Amask,Bmask,
#                                        lsm[0],(bnd_ind[0],bnd_ind[1]))*factor
# Recalculate 30N-40N band to include Med Sea:
#atl_precscl[2] = BandPlusSeaScaled(areas,atl_precbnds[2],med_precnet,Amask,Mmask,
#                                        lsm[0],(bnd_ind[2],bnd_ind[3]))*factor

atl_runscl = FluxBandScaled(areas,Amask*(1-lsm[0]),bnd_ind,atl_run10deg)*factor
# Recalculate 50N-60N band to include Baltic Sea:
#atl_runscl[0] = BandPlusSeaScaled(areas,atl_run10deg[0],bal_run,Amask,Bmask,
#                                        lsm[0],(bnd_ind[0],bnd_ind[1]))*factor
# Recalculate 30N-40N band to include Med Sea:
#atl_runscl[2] = BandPlusSeaScaled(areas,atl_run10deg[2],med_run,Amask,Mmask,
#                                        lsm[0],(bnd_ind[2],bnd_ind[3]))*factor

# Pacific Ocean

pac_EPRscl = FluxBandScaled(areas,Pmask*(1-lsm[0]),bnd_ind,pac_EPRbnds)*factor
pac_evapscl = FluxBandScaled(areas,Pmask*(1-lsm[0]),bnd_ind,pac_evapbnds)*factor
pac_precscl = FluxBandScaled(areas,Pmask*(1-lsm[0]),bnd_ind,pac_precbnds)*factor
pac_runscl = FluxBandScaled(areas,Pmask*(1-lsm[0]),bnd_ind,pac_run10deg)*factor

# Indian Ocean

ind_EPRscl = FluxBandScaled(areas,Imask*(1-lsm[0]),bnd_ind,ind_EPRbnds)*factor
ind_evapscl = FluxBandScaled(areas,Imask*(1-lsm[0]),bnd_ind,ind_evapbnds)*factor
ind_precscl = FluxBandScaled(areas,Imask*(1-lsm[0]),bnd_ind,ind_precbnds)*factor
ind_runscl = FluxBandScaled(areas,Imask*(1-lsm[0]),bnd_ind,ind_run10deg)*factor

# Indian Ocean 30N-40N is nothing, set to nan:
ind_EPRscl[2] = pl.float32('nan'); ind_evapscl[2] = pl.float32('nan');
ind_precscl[2] = pl.float32('nan'); ind_runscl[2] = pl.float32('nan')

# Plot stuff as bar charts
fig, ax = pl.subplots(3,1,figsize=(12,12))
ind = pl.arange(bnd_ind.shape[0]-1)
width = 0.2 # width of each bar
bnd_lbls = ['50N-60N','40N-50N','30N-40N','20N-30N','10N-20N','0-10N','10S-0',
            '20S-10S','30S-20S']
ax1 = pl.subplot(3,1,1)
#for i in range(bnd_ind.shape[0]-1):
q1 = ax1.bar(ind,atl_EPRscl,width,color='k')
q2 = ax1.bar(ind+width,atl_evapscl,width,color='red')
q3 = ax1.bar(ind+2*width,atl_precscl,width,color='blue')
q4 = ax1.bar(ind+3*width,atl_runscl,width,color='green')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.arange(0.4,10.4,1))
ax1.xaxis.set_major_formatter(pl.NullFormatter())
ax1.annotate('(a)',(0,1.02),xycoords='axes fraction',size=25)
pl.ylim(-250,300)
pl.setp(ax1.get_yticklabels(),fontsize=16)
pl.grid(axis='y')
#pl.yticks(pl.linspace(-0.5,1,4))
#ax1.set_ylabel('Freshwater fluxes (cm/yr)',fontsize=16)
ax1.set_title('Atlantic Ocean',fontsize=24)

ax2 = pl.subplot(3,1,2)
w1 = ax2.bar(ind,pac_EPRscl,width,color='k')
w2 = ax2.bar(ind+width,pac_evapscl,width,color='red')
w3 = ax2.bar(ind+2*width,pac_precscl,width,color='blue')
w4 = ax2.bar(ind+3*width,pac_runscl,width,color='green')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.arange(0.4,10.4,1))
ax2.xaxis.set_major_formatter(pl.NullFormatter())
ax2.annotate('(b)',(0,1.02),xycoords='axes fraction',size=25)
pl.ylim(-250,300)
pl.grid(axis='y')
pl.setp(ax2.get_yticklabels(),fontsize=16)
#pl.yticks(pl.linspace(-0.5,1,4))
ax2.set_ylabel('Surface water fluxes (cm/yr)',fontsize=24)
ax2.set_title('Pacific Ocean',fontsize=24)

ax3 = pl.subplot(3,1,3)
s1 = ax3.bar(ind,ind_EPRscl,width,color='k')
s2 = ax3.bar(ind+width,ind_evapscl,width,color='red')
s3 = ax3.bar(ind+2*width,ind_precscl,width,color='blue')
s4 = ax3.bar(ind+3*width,ind_runscl,width,color='green')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.arange(0.4,10.4,1))
ax3.set_xticklabels(bnd_lbls,rotation=45,fontsize=20)
ax3.annotate('(c)',(0,1.02),xycoords='axes fraction',size=25)
#ax3.set_xlabel('Latitude bands',fontsize=15)
pl.ylim(-250,300)
pl.setp(ax3.get_yticklabels(),fontsize=16)
pl.grid(axis='y')
#pl.yticks(pl.linspace(-0.5,1,4))
#ax3.set_ylabel('Freshwater fluxes (cm/yr)',fontsize=16)
ax3.set_title('Indian Ocean',fontsize=24)
ax3.legend((s1[0],s2[0],s3[0],s4[0]),('$\\bar{e}-\\bar{p}-\\bar{r}$','$\\bar{e}$','$\\bar{p}$','$\\bar{r}$'),
           loc=2,fontsize=25)

#fig.suptitle('Ocean freshwater fluxes in 10$^\circ$ latitude bands scaled by area',
#             fontsize=17,x=0.55,y=0.95)
pl.subplots_adjust(left=0.10,right=0.98)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/fluxes_EPR_10degscl.png')


# Plot E-P-R bands for all oceans together:
fig, ax = pl.subplots(1,1,figsize=(9,8))
atloc = ax.bar(ind,atl_EPRscl,width,color='k')
pacoc = ax.bar(ind+width,pac_EPRscl,width,color='red')
indoc = ax.bar(ind+2*width,ind_EPRscl,width,color='blue')
pl.grid(axis='y')
pl.xticks(pl.arange(0.3,10.3,1))
ax.set_xticklabels(bnd_lbls,rotation=45)
ax.set_xlabel('Latitude bands',fontsize=16)
ax.set_ylabel('$E-P-R$  (cm/yr)',fontsize=16)
ax.legend((atloc[0],pacoc[0],indoc[0]),('Atlantic','Pacific','Indian'),loc=2)
#pl.title('Net freshwater flux in 10$^\circ$ latitude bands scaled by area')
pl.subplots_adjust(top=0.95,bottom=0.15)


#-----------CALCULATE THE DIFFERENCES BETWEEN ATLANTIC AND PACIFIC-------------

EPR_diff = pac_EPRscl - atl_EPRscl
evap_diff = pac_evapscl - atl_evapscl
prec_diff = pac_precscl - atl_precscl
run_diff = pac_runscl - atl_runscl

fig, ax = pl.subplots(figsize=(10,10))
ax.bar(ind,EPR_diff,width,color='k',label='$\\bar{e}-\\bar{p}-\\bar{r}$')
ax.bar(ind+width,evap_diff,width,color='red',label='$\\bar{e}$')
ax.bar(ind+2*width,prec_diff,width,color='b',label='$\\bar{p}$')
ax.bar(ind+3*width,run_diff,width,color='g',label='$\\bar{r}$')
pl.grid(axis='y')
pl.xticks(pl.arange(0.4,10.4,1))
ax.set_xticklabels(bnd_lbls,rotation=45,fontsize=20)
#ax.set_xlabel('Latitude bands',fontsize=16)
ax.set_ylabel('Surface water fluxes  (cm/yr)',fontsize=24,labelpad=-10)
pl.setp(ax.get_yticklabels(),fontsize=24)
ax.legend(loc=0,fontsize=25)
pl.title('Pacific minus Atlantic',fontsize=25)
pl.subplots_adjust(bottom=0.11)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/fluxdiffs_10deg.png')