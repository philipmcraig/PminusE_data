# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 15:19:39 2015

Code to calculate net E-P(-R) for 5 degree latitude bands for each basin and
plot as bar charts.

Last updated: 18/09/2015 5:23PM 1st commit

@author: np838619
"""

# import stuff up here
from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import pandas
from functions import *

# read in data here
# basin masks are in netcdf files
maskfile = Dataset('tcdq_basins_60.nc','r')
Amask = maskfile.variables['maskA'][:] # atlantic basin mask
Bmask = maskfile.variables['maskB'][:] # baltic basin mask
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

# runoff data:
runoff = pandas.read_csv('DT03_runoff.csv',sep=',',header=1,skiprows=0,
    names=['row','lat','Atl','Pac','Ind','Glo'])
# to get the values of a column type runoff.name e.g. runoff.Pac
run_atl = runoff.Atl
run_pac = runoff.Pac
run_ind = runoff.Ind

# Make a new mask for Atlantic, Mediterranean & Baltic:
ABmask = Amask + Bmask #+ Mmask

# Multiply each flux array by each mask to extract correct flux:
atl_EmP = EmP*ABmask
atl_evap = evap*ABmask
atl_prec = prec*ABmask
pac_EmP = EmP*Pmask
pac_evap = evap*Pmask
pac_prec = prec*Pmask
ind_EmP = EmP*Imask
ind_evap = evap*Imask
ind_prec = prec*Imask

# define arrays for 5 degree latitude bands for each basin here
bands = pl.linspace(-35,65,21)
# need to convert these to nearest index on ERA-Interim latitude grid:
bnd_ind = pl.zeros_like(bands)
for i in range(bands.shape[0]):
    bnd_ind[i] = NearestIndex(lat,bands[i])
bnd_ind = pl.flipud(bnd_ind) # flip this upside down so that it starts with the smallest lat index

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

# Use NetEmPBands functions and perhaps EPRBands & TotalRunoff to calculate net
# fluxes in each band:

# Atlantic Ocean:
atl_EmPbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(atl_EmP)*factor,rho)
atl_evapbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(atl_evap)*factor,rho)
atl_precbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(atl_prec)*factor,rho)

atl_run5deg = pl.zeros_like(atl_EmPbnds)
for i in range(len(bnd_ind)-1):
    latpair = (bnd_ind[i],bnd_ind[i+1])
    atl_run5deg[i] = TotalRunoff(runoff.lat,lat,latpair,run_atl)

atl_EPRbnds = EPRbands(atl_EmPbnds,bnd_ind,runoff.lat,lat,run_atl)

# Baltic Sea enters the North Sea in 55N-60N band, add the 60N-65N band onto
# the 55N-60N band and set the 60N-65N band to nan:
atl_EPRbnds[1] = atl_EPRbnds[0] + atl_EPRbnds[1]
atl_EPRbnds[0] = pl.float64('nan')
atl_evapbnds[1] = atl_evapbnds[0] + atl_evapbnds[1]
atl_evapbnds[0] = pl.float64('nan')
atl_precbnds[1] = atl_precbnds[0] + atl_precbnds[1]
atl_precbnds[0] = pl.float64('nan')
atl_run5deg[1] = atl_run5deg[1] + atl_run5deg[0]
atl_run5deg[0] = pl.float64('nan')

# Pacific Ocean:
pac_EmPbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(pac_EmP)*factor,rho)
pac_evapbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(pac_evap)*factor,rho)
pac_precbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(pac_prec)*factor,rho)

pac_run5deg = pl.zeros_like(pac_EmPbnds)
for i in range(len(bnd_ind)-1):
    latpair = (bnd_ind[i],bnd_ind[i+1])
    pac_run5deg[i] = TotalRunoff(runoff.lat,lat,latpair,run_pac)

pac_EPRbnds = EPRbands(pac_EmPbnds,bnd_ind,runoff.lat,lat,run_pac)

# No data for Pacific 35S-30S at the moment, will probably change this, but for
# now set the 35S-30S band to nan:
pac_EPRbnds[-1] = pl.float64('nan')
pac_evapbnds[-1] = pl.float64('nan')
pac_precbnds[-1] = pl.float64('nan')
pac_run5deg[-1] = pl.float64('nan')

# Indian Ocean:
ind_EmPbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(ind_EmP)*factor,rho)
ind_evapbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(ind_evap)*factor,rho)
ind_precbnds = NetEmPBands(bnd_ind,areas,pl.squeeze(ind_prec)*factor,rho)

ind_run5deg = pl.zeros_like(ind_EmPbnds)
for i in range(len(bnd_ind)-1):
    latpair = (bnd_ind[i],bnd_ind[i+1])
    ind_run5deg[i] = TotalRunoff(runoff.lat,lat,latpair,run_ind)

ind_EPRbnds = EPRbands(ind_EmPbnds,bnd_ind,runoff.lat,lat,run_ind)

# Indian Ocean stops at 30N, so set everything above 30N to nan:
ind_EPRbnds[:7] = pl.float64('nan')
ind_evapbnds[:7] = pl.float64('nan')
ind_precbnds[:7] = pl.float64('nan')
ind_run5deg[:7] = pl.float64('nan')


#---------------------------SCALE EVERYTHING BY AREA---------------------------

factor = (10**6)*100*86400*365 # factor to rescale from Sv/m^2 to cm/yr

# Atlantic Ocean

atl_EPRscl = FluxBandScaled(areas,ABmask,bnd_ind,atl_EPRbnds)*factor
atl_evapscl = FluxBandScaled(areas,ABmask,bnd_ind,atl_evapbnds)*factor
atl_precscl = FluxBandScaled(areas,ABmask,bnd_ind,atl_precbnds)*factor
atl_runscl = FluxBandScaled(areas,ABmask,bnd_ind,atl_run5deg)*factor

# Pacific Ocean

pac_EPRscl = FluxBandScaled(areas,Pmask,bnd_ind,pac_EPRbnds)*factor
pac_evapscl = FluxBandScaled(areas,Pmask,bnd_ind,pac_evapbnds)*factor
pac_precscl = FluxBandScaled(areas,Pmask,bnd_ind,pac_precbnds)*factor
pac_runscl = FluxBandScaled(areas,Pmask,bnd_ind,pac_run5deg)*factor

# Indian Ocean

ind_EPRscl = FluxBandScaled(areas,Imask,bnd_ind,ind_EPRbnds)*factor
ind_evapscl = FluxBandScaled(areas,Imask,bnd_ind,ind_evapbnds)*factor
ind_precscl = FluxBandScaled(areas,Imask,bnd_ind,ind_precbnds)*factor
ind_runscl = FluxBandScaled(areas,Imask,bnd_ind,ind_run5deg)*factor


# Plot stuff as bar charts
fig, ax = pl.subplots(3,1,figsize=(12,12))
ind = pl.arange(bnd_ind.shape[0]-1)
width = 0.2 # width of each bar
bnd_lbls = ['60N-65N','55N-60N','50N-55N','45N-50N','40N-45N','35N-40N','30N-35N',
            '25N-30N','20N-25N','15N-20N','10N-15N','5N-10N','0-5N','5S-0','10S-5S',
            '15S-10S','20S-15S','25S-20S','30S-25S','35S-30S']
ax1 = pl.subplot(3,1,1)
#for i in range(bnd_ind.shape[0]-1):
q1 = ax1.bar(ind,atl_EPRscl,width,color='k')
q2 = ax1.bar(ind+width,atl_evapscl,width,color='red')
q3 = ax1.bar(ind+2*width,atl_precscl,width,color='blue')
q4 = ax1.bar(ind+3*width,atl_runscl,width,color='green')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.linspace(2*width,19+2*width,20))
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(-250,300)
#pl.yticks(pl.linspace(-0.5,1,4))
ax1.set_ylabel('Freshwater fluxes (cm/yr)',fontsize=16)
ax1.set_title('Atlantic Ocean')

ax2 = pl.subplot(3,1,2)
w1 = ax2.bar(ind,pac_EPRscl,width,color='k')
w2 = ax2.bar(ind+width,pac_evapscl,width,color='red')
w3 = ax2.bar(ind+2*width,pac_precscl,width,color='blue')
w4 = ax2.bar(ind+3*width,pac_runscl,width,color='green')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.linspace(3*width/2,19+3*width/2,20))
ax2.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(-250,300)
#pl.yticks(pl.linspace(-0.5,1,4))
ax2.set_ylabel('Freshwater fluxes (cm/yr)',fontsize=16)
ax2.set_title('Pacific Ocean')

ax3 = pl.subplot(3,1,3)
s1 = ax3.bar(ind,ind_EPRscl,width,color='k')
s2 = ax3.bar(ind+width,ind_evapscl,width,color='red')
s3 = ax3.bar(ind+2*width,ind_precscl,width,color='blue')
s4 = ax3.bar(ind+3*width,ind_runscl,width,color='green')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.linspace(3*width/2,19+3*width/2,20))
ax3.set_xticklabels(bnd_lbls,rotation=45)
ax3.set_xlabel('Latitude bands',fontsize=15)
pl.ylim(-250,300)
#pl.yticks(pl.linspace(-0.5,1,4))
ax3.set_ylabel('Freshwater fluxes (cm/yr)',fontsize=16)
ax3.set_title('Indian Ocean')
ax3.legend((s1[0],s2[0],s3[0],s4[0]),('$E-P-R$','$E$','$P$','$R$'),loc=2)

fig.suptitle('Ocean freshwater fluxes in 5$^\circ$ latitude bands scaled by area',
             fontsize=17,x=0.55,y=0.95)
pl.subplots_adjust(left=0.08,right=0.98)