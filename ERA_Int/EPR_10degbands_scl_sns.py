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
lsm = maskfile.variables['LSM'][0] # ERA-Interim land-sea mask
maskfile.close()

# add the Atlantic & Baltic masks together:
ABMmask = Amask + Bmask + Mmask

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

evap_tot = pl.squeeze(evap_tot); prec_tot = pl.squeeze(prec_tot)


# define arrays for 10 degree latitude bands for each basin here
bands = pl.linspace(60,-30,10)
# need to convert these to nearest index on ERA-Interim latitude grid:
bnd_ind = pl.zeros_like(bands)
for i in range(bands.shape[0]):
    bnd_ind[i] = NearestIndex(lat,bands[i])


# RUNOFF DATA: 
# read the 4 csv file with the climatological monthly means for each basin
# DT02 runoff into the Atlantic Ocean in 1deg bands:
atl_run = pandas.read_csv('DT02_Atl_mnths.csv',sep=',',header=1,skiprows=0,
    names=['row','lat','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
           'Oct','Nov','Dec'])

# DT02 runoff into the Medditeranean Sea (one value per month)
med_run = pandas.read_csv('DT02_Med_mnths.csv',sep=',',header=1,skiprows=0,
    names=['row','lat','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
           'Oct','Nov','Dec'])

# DT02 runoff into the Pacific Ocean in 1deg bands:
pac_run = pandas.read_csv('DT02_Pac_mnths.csv',sep=',',header=1,skiprows=0,
    names=['row','lat','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
           'Oct','Nov','Dec'])

# DT02 runoff into the Indian Ocean in 1deg bands:
ind_run = pandas.read_csv('DT02_Ind_mnths.csv',sep=',',header=1,skiprows=0,
    names=['row','lat','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
           'Oct','Nov','Dec'])

# stick it all into arrays

Prun_mnths = pl.array([pac_run.Jan,pac_run.Feb,pac_run.Mar,pac_run.Apr,
                       pac_run.May,pac_run.Jun,pac_run.Jul,pac_run.Aug,
                       pac_run.Sep,pac_run.Oct,pac_run.Nov,pac_run.Dec])

Irun_mnths = pl.array([ind_run.Jan,ind_run.Feb,ind_run.Mar,ind_run.Apr,
                       ind_run.May,ind_run.Jun,ind_run.Jul,ind_run.Aug,
                       ind_run.Sep,ind_run.Oct,ind_run.Nov,ind_run.Dec])

Arun_DJF = pl.array([atl_run.Dec,atl_run.Jan,atl_run.Feb])
Arun_MAM = pl.array([atl_run.Mar,atl_run.Apr,atl_run.May])
Arun_JJA = pl.array([atl_run.Jun,atl_run.Jul,atl_run.Aug])
Arun_SON = pl.array([atl_run.Sep,atl_run.Oct,atl_run.Nov])

Prun_DJF = pl.array([pac_run.Dec,pac_run.Jan,pac_run.Feb])
Prun_MAM = pl.array([pac_run.Mar,pac_run.Apr,pac_run.May])
Prun_JJA = pl.array([pac_run.Jun,pac_run.Jul,pac_run.Aug])
Prun_SON = pl.array([pac_run.Sep,pac_run.Oct,pac_run.Nov])

Irun_DJF = pl.array([ind_run.Dec,ind_run.Jan,ind_run.Feb])
Irun_MAM = pl.array([ind_run.Mar,ind_run.Apr,ind_run.May])
Irun_JJA = pl.array([ind_run.Jun,ind_run.Jul,ind_run.Aug])
Irun_SON = pl.array([ind_run.Sep,ind_run.Oct,ind_run.Nov])

# season runoff means:
Arun_DJF = pl.mean(Arun_DJF,axis=0); Arun_MAM = pl.mean(Arun_MAM,axis=0)
Arun_JJA = pl.mean(Arun_JJA,axis=0); Arun_SON = pl.mean(Arun_SON,axis=0)
Arun_sns = pl.array([Arun_DJF,Arun_MAM,Arun_JJA,Arun_SON])

Prun_DJF = pl.mean(Prun_DJF,axis=0); Prun_MAM = pl.mean(Prun_MAM,axis=0)
Prun_JJA = pl.mean(Prun_JJA,axis=0); Prun_SON = pl.mean(Prun_SON,axis=0)
Prun_sns = pl.array([Prun_DJF,Prun_MAM,Prun_JJA,Prun_SON])

Irun_DJF = pl.mean(Irun_DJF,axis=0); Irun_MAM = pl.mean(Irun_MAM,axis=0)
Irun_JJA = pl.mean(Irun_JJA,axis=0); Irun_SON = pl.mean(Irun_SON,axis=0)
Irun_sns = pl.array([Irun_DJF,Irun_MAM,Irun_JJA,Irun_SON])

atl_run_tot = pl.zeros([4,bnd_ind.shape[0]-1])
pac_run_tot = pl.zeros([4,bnd_ind.shape[0]-1])
ind_run_tot = pl.zeros([4,bnd_ind.shape[0]-1])

# calculate total seasonal runoff for latitude bands:
for m in range(4):
    for b in range(atl_run_tot.shape[1]):
        atl_run_tot[m,b] = TotalRunoff(atl_run.lat,lat,(bnd_ind[b],bnd_ind[b+1]),Arun_sns[m])
        pac_run_tot[m,b] = TotalRunoff(pac_run.lat,lat,(bnd_ind[b],bnd_ind[b+1]),Prun_sns[m])
        ind_run_tot[m,b] = TotalRunoff(ind_run.lat,lat,(bnd_ind[b],bnd_ind[b+1]),Irun_sns[m])


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
#med_EmPnet = NetEmPCalc(areas,med_EmP*factor,rho)
#med_evapnet = NetEmPCalc(areas,med_evap*factor,rho)
#med_precnet = NetEmPCalc(areas,med_prec*factor,rho)
## Baltic Sea:
#bal_EmPnet = NetEmPCalc(areas,bal_EmP*factor,rho)
#bal_evapnet = NetEmPCalc(areas,bal_evap*factor,rho)
#bal_precnet = NetEmPCalc(areas,bal_prec*factor,rho)

# Calculate runoff in Med & Baltic Seas:
    # Difference between global runoff & sum of Atl, Pac & Ind runoff
#API = run_atl + run_pac + run_ind
#rundiff = runoff.Glo - API
#    # Find indices of max/min latitudes of Med & Baltic Seas:
#med_nth = NearestIndex(runoff.lat,47.3); med_sth = NearestIndex(runoff.lat,30.5)
#bal_nth = NearestIndex(runoff.lat,66.0); bal_sth = NearestIndex(runoff.lat,53.8)
#    # calculate runoff
#med_run = pl.sum(rundiff[med_nth:med_sth+1])
#bal_run = pl.sum(rundiff[bal_nth:bal_sth+1])
#
## Calculate E-P-R for Med and Baltic Seas:
#med_EPRnet = med_EmPnet - med_run
#bal_EPRnet = bal_EmPnet - bal_run


evap_sns = SeasonalMeans(evap_tot)
prec_sns = SeasonalMeans(prec_tot)
EmP_sns = SeasonalMeans(evap_tot-prec_tot)

atl_evap = evap_sns*ABMmask; atl_prec = prec_sns*ABMmask; atl_EmP = EmP_sns*ABMmask
pac_evap = evap_sns*Pmask; pac_prec = prec_sns*Pmask; pac_EmP = EmP_sns*Pmask
ind_evap = evap_sns*Imask; ind_prec = prec_sns*Imask; ind_EmP = EmP_sns*Imask

atl_evp_bnds = pl.zeros([4,bands.size-1]); atl_prc_bnds = pl.zeros_like(atl_evp_bnds)
pac_evp_bnds = pl.zeros_like(atl_evp_bnds); pac_prc_bnds = pl.zeros_like(atl_evp_bnds)
ind_evp_bnds = pl.zeros_like(atl_evp_bnds); ind_prc_bnds = pl.zeros_like(ind_evp_bnds)
atl_EmP_bnds = pl.zeros_like(pac_evp_bnds); pac_EmP_bnds = pl.zeros_like(atl_evp_bnds)
ind_EmP_bnds = pl.zeros_like(pac_evp_bnds)

for i in range(4):
    atl_evp_bnds[i] = NetEmPBands(bnd_ind,areas,atl_evap[i]*factor,rho)
    atl_prc_bnds[i] = NetEmPBands(bnd_ind,areas,atl_prec[i]*factor,rho)
    atl_EmP_bnds[i] = NetEmPBands(bnd_ind,areas,atl_EmP[i]*factor,rho)
    pac_evp_bnds[i] = NetEmPBands(bnd_ind,areas,pac_evap[i]*factor,rho)
    pac_prc_bnds[i] = NetEmPBands(bnd_ind,areas,pac_prec[i]*factor,rho)
    pac_EmP_bnds[i] = NetEmPBands(bnd_ind,areas,pac_EmP[i]*factor,rho)
    ind_evp_bnds[i] = NetEmPBands(bnd_ind,areas,ind_evap[i]*factor,rho)
    ind_prc_bnds[i] = NetEmPBands(bnd_ind,areas,ind_prec[i]*factor,rho)
    ind_EmP_bnds[i] = NetEmPBands(bnd_ind,areas,ind_EmP[i]*factor,rho)


atl_EPR_bnds = pl.zeros_like(atl_EmP_bnds); pac_EPR_bnds = pl.zeros_like(pac_EmP_bnds)
ind_EPR_bnds = pl.zeros_like(ind_EmP_bnds)
for i in range(4):
    atl_EPR_bnds[i] = EPRbands(atl_EmP_bnds[i],bnd_ind,atl_run.lat,lat,Arun_sns[i])
    pac_EPR_bnds[i] = EPRbands(pac_EmP_bnds[i],bnd_ind,pac_run.lat,lat,Prun_sns[i])
    ind_EPR_bnds[i] = EPRbands(ind_EmP_bnds[i],bnd_ind,ind_run.lat,lat,Irun_sns[i])


#---------------------------SCALE EVERYTHING BY AREA---------------------------

factor = (10**6)*100*86400*365 # factor to rescale from Sv/m^2 to cm/yr

# Atlantic Ocean

#atl_EPRscl = FluxBandScaled(areas,Amask*(1-lsm[0]),bnd_ind,atl_EPRbnds)*factor
# Recalculate 50N-60N band to include Baltic Sea:
#atl_EPRscl[0] = BandPlusSeaScaled(areas,atl_EPRbnds[0],bal_EPRnet,Amask,Bmask,
#                                        lsm[0],(bnd_ind[0],bnd_ind[1]))*factor
# Recalculate 30N-40N band to include Med Sea:
#atl_EPRscl[2] = BandPlusSeaScaled(areas,atl_EPRbnds[2],med_EPRnet,Amask,Mmask,
#                                        lsm[0],(bnd_ind[2],bnd_ind[3]))*factor

#atl_evapscl = FluxBandScaled(areas,Amask*(1-lsm[0]),bnd_ind,atl_evapbnds)*factor
# Recalculate 50N-60N band to include Baltic Sea:
#atl_evapscl[0] = BandPlusSeaScaled(areas,atl_evapbnds[0],bal_evapnet,Amask,Bmask,
#                                        lsm[0],(bnd_ind[0],bnd_ind[1]))*factor
# Recalculate 30N-40N band to include Med Sea:
#atl_evapscl[2] = BandPlusSeaScaled(areas,atl_evapbnds[2],med_evapnet,Amask,Mmask,
#                                       lsm[0],(bnd_ind[2],bnd_ind[3]))*factor

#atl_precscl = FluxBandScaled(areas,Amask*(1-lsm[0]),bnd_ind,atl_precbnds)*factor
# Recalculate 50N-60N band to include Baltic Sea:
#atl_precscl[0] = BandPlusSeaScaled(areas,atl_precbnds[0],bal_precnet,Amask,Bmask,
#                                        lsm[0],(bnd_ind[0],bnd_ind[1]))*factor
# Recalculate 30N-40N band to include Med Sea:
#atl_precscl[2] = BandPlusSeaScaled(areas,atl_precbnds[2],med_precnet,Amask,Mmask,
#                                        lsm[0],(bnd_ind[2],bnd_ind[3]))*factor

#atl_runscl = FluxBandScaled(areas,Amask*(1-lsm[0]),bnd_ind,atl_run10deg)*factor
# Recalculate 50N-60N band to include Baltic Sea:
#atl_runscl[0] = BandPlusSeaScaled(areas,atl_run10deg[0],bal_run,Amask,Bmask,
#                                        lsm[0],(bnd_ind[0],bnd_ind[1]))*factor
# Recalculate 30N-40N band to include Med Sea:
#atl_runscl[2] = BandPlusSeaScaled(areas,atl_run10deg[2],med_run,Amask,Mmask,
#                                        lsm[0],(bnd_ind[2],bnd_ind[3]))*factor

# Pacific Ocean

#pac_EPRscl = FluxBandScaled(areas,Pmask*(1-lsm[0]),bnd_ind,pac_EPRbnds)*factor
#pac_evapscl = FluxBandScaled(areas,Pmask*(1-lsm[0]),bnd_ind,pac_evapbnds)*factor
#pac_precscl = FluxBandScaled(areas,Pmask*(1-lsm[0]),bnd_ind,pac_precbnds)*factor
#pac_runscl = FluxBandScaled(areas,Pmask*(1-lsm[0]),bnd_ind,pac_run10deg)*factor

# Indian Ocean

#ind_EPRscl = FluxBandScaled(areas,Imask*(1-lsm[0]),bnd_ind,ind_EPRbnds)*factor
#ind_evapscl = FluxBandScaled(areas,Imask*(1-lsm[0]),bnd_ind,ind_evapbnds)*factor
#ind_precscl = FluxBandScaled(areas,Imask*(1-lsm[0]),bnd_ind,ind_precbnds)*factor
#ind_runscl = FluxBandScaled(areas,Imask*(1-lsm[0]),bnd_ind,ind_run10deg)*factor

# Indian Ocean 30N-40N is nothing, set to nan:
#ind_EPRscl[2] = pl.float32('nan'); ind_evapscl[2] = pl.float32('nan');
#ind_precscl[2] = pl.float32('nan'); ind_runscl[2] = pl.float32('nan')


atl_evp_scl = pl.zeros_like(atl_evp_bnds); atl_prc_scl = pl.zeros_like(atl_evp_bnds)
pac_evp_scl = pl.zeros_like(atl_evp_bnds); pac_prc_scl = pl.zeros_like(atl_evp_bnds)
ind_evp_scl = pl.zeros_like(atl_evp_bnds); ind_prc_scl = pl.zeros_like(atl_evp_bnds)
atl_EmP_scl = pl.zeros_like(atl_evp_bnds); pac_EmP_scl = pl.zeros_like(atl_evp_bnds)
ind_EmP_scl = pl.zeros_like(atl_evp_bnds)
atl_EPR_scl = pl.zeros_like(atl_evp_bnds); pac_EPR_scl = pl.zeros_like(pac_evp_bnds)
ind_EPR_scl = pl.zeros_like(ind_prc_bnds)
atl_run_scl = pl.zeros_like(atl_prc_bnds); pac_run_scl = pl.zeros_like(pac_EPR_bnds)
ind_run_scl = pl.zeros_like(ind_EmP_bnds)

for i in range(4):
    atl_evp_scl[i] = FluxBandScaled(areas,Amask*(1-lsm),bnd_ind,atl_evp_bnds[i])*factor
    atl_prc_scl[i] = FluxBandScaled(areas,Amask*(1-lsm),bnd_ind,atl_prc_bnds[i])*factor
    atl_EmP_scl[i] = FluxBandScaled(areas,Amask*(1-lsm),bnd_ind,atl_EmP_bnds[i])*factor
    atl_EPR_scl[i] = FluxBandScaled(areas,Amask*(1-lsm),bnd_ind,atl_EPR_bnds[i])*factor
    atl_run_scl[i] = FluxBandScaled(areas,Amask*(1-lsm),bnd_ind,atl_run_tot[i])*factor
    pac_evp_scl[i] = FluxBandScaled(areas,Pmask*(1-lsm),bnd_ind,pac_evp_bnds[i])*factor
    pac_prc_scl[i] = FluxBandScaled(areas,Pmask*(1-lsm),bnd_ind,pac_prc_bnds[i])*factor
    pac_EmP_scl[i] = FluxBandScaled(areas,Pmask*(1-lsm),bnd_ind,pac_EmP_bnds[i])*factor
    pac_EPR_scl[i] = FluxBandScaled(areas,Pmask*(1-lsm),bnd_ind,pac_EPR_bnds[i])*factor
    pac_run_scl[i] = FluxBandScaled(areas,Pmask*(1-lsm),bnd_ind,pac_run_tot[i])*factor
    ind_evp_scl[i] = FluxBandScaled(areas,Imask*(1-lsm),bnd_ind,ind_evp_bnds[i])*factor
    ind_prc_scl[i] = FluxBandScaled(areas,Imask*(1-lsm),bnd_ind,ind_prc_bnds[i])*factor
    ind_EmP_scl[i] = FluxBandScaled(areas,Imask*(1-lsm),bnd_ind,ind_EmP_bnds[i])*factor
    ind_EPR_scl[i] = FluxBandScaled(areas,Imask*(1-lsm),bnd_ind,ind_EPR_bnds[i])*factor
    ind_run_scl[i] = FluxBandScaled(areas,Imask*(1-lsm),bnd_ind,ind_run_tot[i])*factor

# Indian Ocean 30N-40N is nothing, set to nan:
ind_EPR_scl[:,2] = pl.float32('nan'); ind_evp_scl[:,2] = pl.float32('nan')
ind_prc_scl[:,2] = pl.float32('nan'); ind_run_scl[:,2] = pl.float32('nan')

# Plot stuff as bar charts
fig, ax = pl.subplots(2,2,figsize=(14,11),sharex=True)
ind = pl.arange(bnd_ind.shape[0]-1)
width = 0.2 # width of each bar
bnd_lbls = ['50N-60N','40N-50N','30N-40N','20N-30N','10N-20N','0-10N','10S-0',
            '20S-10S','30S-20S']
ax1 = pl.subplot(221)
q1 = ax1.bar(ind,ind_EPR_scl[0],width,color='k',label='$\\bar{e}-\\bar{p}-\\bar{r}$')
q2 = ax1.bar(ind+width,ind_evp_scl[0],width,color='red',label='$\\bar{e}$')
q3 = ax1.bar(ind+2*width,ind_prc_scl[0],width,color='blue',label='$\\bar{p}$')
q4 = ax1.bar(ind+3*width,ind_run_scl[0],width,color='green',label='$\\bar{r}$')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.arange(0.4,10.4,1))
ax1.xaxis.set_major_formatter(pl.NullFormatter())
ax1.annotate('(a)',(0,1.02),xycoords='axes fraction',size=25)
pl.ylim(-300,300)
pl.setp(ax1.get_yticklabels(),fontsize=16)
pl.grid(axis='y'); ax1.legend(loc=3,fontsize=20)
#pl.yticks(pl.linspace(-0.5,1,4))
ax1.set_ylabel('cm/yr',fontsize=18,labelpad=-2)
ax1.set_title('DJF',fontsize=24)

ax2 = pl.subplot(222)
w1 = ax2.bar(ind,ind_EPR_scl[1],width,color='k')
w2 = ax2.bar(ind+width,ind_evp_scl[1],width,color='red')
w3 = ax2.bar(ind+2*width,ind_prc_scl[1],width,color='blue')
w4 = ax2.bar(ind+3*width,ind_run_scl[1],width,color='green')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.arange(0.4,10.4,1))
ax2.xaxis.set_major_formatter(pl.NullFormatter())
ax2.annotate('(b)',(0,1.02),xycoords='axes fraction',size=25)
pl.ylim(-300,300)
pl.grid(axis='y')
pl.setp(ax2.get_yticklabels(),fontsize=16)
#pl.yticks(pl.linspace(-0.5,1,4))
#ax2.set_ylabel('Surface water fluxes (cm/yr)',fontsize=24)
ax2.set_title('MAM',fontsize=24)

ax3 = pl.subplot(223)
s1 = ax3.bar(ind,ind_EPR_scl[2],width,color='k',label='$\\bar{e}-\\bar{p}-\\bar{r}$')
s2 = ax3.bar(ind+width,ind_evp_scl[2],width,color='red',label='$\\bar{e}$')
s3 = ax3.bar(ind+2*width,ind_prc_scl[2],width,color='blue',label='$\\bar{p}$')
s4 = ax3.bar(ind+3*width,ind_run_scl[2],width,color='green',label='$\\bar{r}$')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.arange(0.4,10.4,1))
ax3.set_xticklabels(bnd_lbls,rotation=22,fontsize=14)
ax3.annotate('(c)',(0,1.02),xycoords='axes fraction',size=25)
#ax3.set_xlabel('Latitude bands',fontsize=15)
pl.ylim(-300,300)
pl.setp(ax3.get_yticklabels(),fontsize=16)
pl.grid(axis='y')
#pl.yticks(pl.linspace(-0.5,1,4))
ax3.set_ylabel('cm/yr',fontsize=18,labelpad=-2)
ax3.set_title('JJA',fontsize=24)


ax4 = pl.subplot(224)
s1 = ax4.bar(ind,ind_EPR_scl[3],width,color='k',label='$\\bar{e}-\\bar{p}-\\bar{r}$')
s2 = ax4.bar(ind+width,ind_evp_scl[3],width,color='red',label='$\\bar{e}$')
s3 = ax4.bar(ind+2*width,ind_prc_scl[3],width,color='blue',label='$\\bar{p}$')
s4 = ax4.bar(ind+3*width,ind_run_scl[3],width,color='green',label='$\\bar{r}$')
pl.axhline(y=0,color='lightgrey',ls='--',linewidth=0.5) # add line to show y=0
pl.xticks(pl.arange(0.4,10.4,1))
ax4.set_xticklabels(bnd_lbls,rotation=22,fontsize=14)
ax4.annotate('(d)',(0,1.02),xycoords='axes fraction',size=25)
#ax3.set_xlabel('Latitude bands',fontsize=15)
pl.ylim(-300,300)
pl.setp(ax4.get_yticklabels(),fontsize=16)
pl.grid(axis='y')
#ax3.set_ylabel('Freshwater fluxes (cm/yr)',fontsize=16)
ax4.set_title('SON',fontsize=24)

pl.suptitle('Indian',y=0.95,fontsize=25)
pl.subplots_adjust(left=0.10,right=0.98)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/ind_10degscl_sns.png')



#-----------CALCULATE THE DIFFERENCES BETWEEN ATLANTIC AND PACIFIC-------------

EPR_diff = pac_EPR_scl - atl_EPR_scl
evap_diff = pac_evp_scl - atl_evp_scl
prec_diff = pac_prc_scl - atl_prc_scl
run_diff = pac_run_scl - atl_run_scl

fig, ax = pl.subplots(2,2,figsize=(14,11))

ax1 = pl.subplot(221)
ax1.bar(ind,EPR_diff[0],width,color='k',label='$\\bar{e}-\\bar{p}-\\bar{r}$')
ax1.bar(ind+width,evap_diff[0],width,color='red',label='$\\bar{e}$')
ax1.bar(ind+2*width,prec_diff[0],width,color='b',label='$\\bar{p}$')
ax1.bar(ind+3*width,run_diff[0],width,color='g',label='$\\bar{r}$')
pl.grid(axis='y')
pl.xticks(pl.arange(0.4,10.4,1))
ax1.xaxis.set_major_formatter(pl.NullFormatter())
ax1.set_ylabel('cm/yr',fontsize=24,labelpad=-10)
pl.setp(ax1.get_yticklabels(),fontsize=16); pl.ylim(-200,200)
ax1.annotate('(a)',(0,1.02),xycoords='axes fraction',size=25)
ax1.set_title('DJF',fontsize=24)

ax2 = pl.subplot(222)
ax2.bar(ind,EPR_diff[1],width,color='k',label='$\\bar{e}-\\bar{p}-\\bar{r}$')
ax2.bar(ind+width,evap_diff[1],width,color='red',label='$\\bar{e}$')
ax2.bar(ind+2*width,prec_diff[1],width,color='b',label='$\\bar{p}$')
ax2.bar(ind+3*width,run_diff[1],width,color='g',label='$\\bar{r}$')
pl.grid(axis='y')
pl.xticks(pl.arange(0.4,10.4,1))
ax2.xaxis.set_major_formatter(pl.NullFormatter()); pl.ylim(-200,200)
pl.setp(ax2.get_yticklabels(),fontsize=16)
ax2.legend(loc=0,fontsize=17)
ax2.annotate('(b)',(0,1.02),xycoords='axes fraction',size=25)
ax2.set_title('MAM',fontsize=24)

ax3 = pl.subplot(223)
ax3.bar(ind,EPR_diff[2],width,color='k',label='$\\bar{e}-\\bar{p}-\\bar{r}$')
ax3.bar(ind+width,evap_diff[2],width,color='red',label='$\\bar{e}$')
ax3.bar(ind+2*width,prec_diff[2],width,color='b',label='$\\bar{p}$')
ax3.bar(ind+3*width,run_diff[2],width,color='g',label='$\\bar{r}$')
pl.grid(axis='y')
pl.xticks(pl.arange(0.4,10.4,1))
ax3.set_xticklabels(bnd_lbls,rotation=22,fontsize=14)
pl.setp(ax3.get_yticklabels(),fontsize=16); pl.ylim(-200,200)
ax3.set_ylabel('cm/yr',fontsize=24,labelpad=-10)
ax3.annotate('(c)',(0,1.02),xycoords='axes fraction',size=25)
ax3.set_title('JJA',fontsize=24)

ax4 = pl.subplot(224)
ax4.bar(ind,EPR_diff[3],width,color='k',label='$\\bar{e}-\\bar{p}-\\bar{r}$')
ax4.bar(ind+width,evap_diff[3],width,color='red',label='$\\bar{e}$')
ax4.bar(ind+2*width,prec_diff[3],width,color='b',label='$\\bar{p}$')
ax4.bar(ind+3*width,run_diff[3],width,color='g',label='$\\bar{r}$')
pl.grid(axis='y')
pl.xticks(pl.arange(0.4,10.4,1))
ax4.set_xticklabels(bnd_lbls,rotation=22,fontsize=14); pl.ylim(-200,200)
pl.setp(ax4.get_yticklabels(),fontsize=16)
ax4.annotate('(d)',(0,1.02),xycoords='axes fraction',size=25)
ax4.set_title('SON',fontsize=24)

pl.suptitle('Pacific minus Atlantic',y=0.975,fontsize=25)
#pl.subplots_adjust(bottom=0.11)
pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/fluxdiffs_10deg_sns.png')