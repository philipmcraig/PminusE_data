# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 10:19:36 2015

@author: np838619

Code to read in evaporation and total precipitation data for ERA-Interim
accumulated surface forecast data, and calculate the basin-integrated E & P for
each ocean basin for each year from 1979-2014 to show interannual variability.

Last updated: 2/03/2016 6:22PM 12th commit
"""

# import modules up here: pylab, netCDF4, functions
from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from scipy.stats import pearsonr
from functions import *
from matplotlib.lines import Line2D

pl.close('all')
# read in the basin masks from netCDF file
basinfile = Dataset('tcdq_basins_60.nc','r')
atlmask = basinfile.variables['maskA'][:] # mask for Atlantic Ocean
balmask = basinfile.variables['maskB'][:] # mask for Baltic Sea
medmask = basinfile.variables['maskM'][:] # mask for Med Sea
pacmask = basinfile.variables['maskP'][:] # mask for Pacific Ocean
indmask = basinfile.variables['maskI'][:] # mask for Indian Ocean
lon = basinfile.variables['lon'][:] # ERA-Interim longitude
lat = basinfile.variables['lat'][:] # ERA-Interim latitude
basinfile.close()

maskfile = Dataset('tcdq_ann_mean_ERAI.nc','r')
lsm = maskfile.variables['LSM'][:]
maskfile.close()

# add the Atlantic & Baltic masks together:
atlmask = (atlmask + balmask + medmask)*(1-lsm)
pacmask = pacmask*(1-lsm)
indmask = indmask*(1-lsm)

years = pl.linspace(1979,2014,36)
year_input = [str(int(i)) for i in years]


# Read in all the Evaporation and Precipitation data from hafs files:


# empty array for filenames
filenames = pl.zeros([len(years),96],dtype='S17') # Only need every 4th file!

# use PrintFiles function to read & store all the filenames
#loop over years:
for year in range(len(years)):
    #path = path to mm folder + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + str(year_input[year])
    #filenames[year] = PrintFiles(path,type)
    filenames[year] = PrintFiles(path,'hafs')

# get rid of all the unneccessary filenames
lessfiles = pl.zeros([filenames.shape[0],filenames.shape[1]/4],dtype='S17')
for year in range(filenames.shape[0]):
    filelist = []
    for name in range(filenames.shape[1]):
        if '12.nc' in filenames[year,name]:
            filelist.append(filenames[year,name])
    for i in range(len(filelist)):
        lessfiles[year,i] = filelist[i]
lessfiles = pl.sort(lessfiles)

# set up empty arrays for P & E
#empty array for evaporation
evap = pl.zeros([len(years),lessfiles.shape[1],1,1,256,512])
#empty array for precipitation
prec = pl.zeros_like(evap)

# read all the netCDF files to extract P & E
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


evap_tot = pl.zeros([evap.shape[0],evap.shape[1]/2,evap.shape[2],evap.shape[3],
                     evap.shape[4],evap.shape[5]])
prec_tot = pl.zeros_like(evap_tot)

# loop over number of years:
for year in range(evap.shape[0]):
    # loop over number of months:
    for month in range(int(evap.shape[1]/2)):
        evap_tot[year,month] = -1*(evap[year][2*month] + evap[year][2*month+1])
        prec_tot[year,month] = (prec[year][2*month] + prec[year][2*month+1])
   
     
# calculate the mean P & E for each year
evap_yrs_mn = pl.mean(evap_tot,axis=1)
prec_yrs_mn = pl.mean(prec_tot,axis=1)

# Get rid of the 1D axes:
evap_yrs_mn = pl.squeeze(evap_yrs_mn)
prec_yrs_mn = pl.squeeze(prec_yrs_mn)

# WRITE ANNUAL MEANS OF EVAP AND PRECIP TO A NETCDF FILE:
#newnc = Dataset('evapprec_36yrs.nc','w')
#
#lat_dim = newnc.createDimension('lat', 256)
#lon_dim = newnc.createDimension('lon', 512)
#lat_in = newnc.createVariable('lat', pl.float32, ('lat',))
#lat_in.units = 'degrees_north'
#lat_in.long_name = 'latitude'
#lon_in = newnc.createVariable('lon', pl.float32, ('lon',))
#lon_in.units = 'degrees_east'
#lon_in.long_name = 'longitude'
#
#time_dim = newnc.createDimension('time',36)
#time = newnc.createVariable('time', pl.float64, ('time',))
#time.units = 'years'
#time.long_name = 'time'
#
#EV = newnc.createVariable('evap',pl.float64,('time','lat','lon'))
#EV.units = 'm/day'
#EV.standard_name = 'Evaporation'
#
#PR = newnc.createVariable('prec',pl.float64,('time','lat','lon'))
#PR.units = 'm/day'
#PR.standard_name = 'Precipitation'
#
#lat_in[:] = lat # straight from ERA-Interim nc file
#lon_in[:] = lon # straight from ERA-Interim nc file
#EV[:,:,:] = evap_yrs_mn[:]
#PR[:,:,:] = prec_yrs_mn[:]
#
#newnc.close()

# Annual mean E & P
evap_mn = pl.mean(evap_yrs_mn,axis=0)
prec_mn = pl.mean(prec_yrs_mn,axis=0)


oc_evap = evap_yrs_mn*(1-lsm)
oc_prec = prec_yrs_mn*(1-lsm)
oc_evap_mn = evap_mn*(1-lsm)
oc_prec_mn = prec_mn*(1-lsm)

# extract the relevant ocean basins:

atl_evap = evap_yrs_mn*atlmask
atl_prec = prec_yrs_mn*atlmask
atl_evap_mn = evap_mn*atlmask
atl_prec_mn = prec_mn*atlmask

pac_evap = evap_yrs_mn*pacmask
pac_prec = prec_yrs_mn*pacmask
pac_evap_mn = evap_mn*pacmask
pac_prec_mn = prec_mn*pacmask

ind_evap = evap_yrs_mn*indmask
ind_prec = prec_yrs_mn*indmask
ind_evap_mn = evap_mn*indmask
ind_prec_mn = prec_mn*indmask


# Read in all the tcdq data from ggaw files:

# empty array for filenames
filenames = pl.zeros([len(years),12],dtype='S17') # Only need every 4th file!

# use PrintFiles function to read & store all the filenames
#loop over years:
for year in range(len(years)):
    #path = path to mm folder + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + str(year_input[year])
    #filenames[year] = PrintFiles(path,type)
    filenames[year] = PrintFiles(path,'ggaw')

# set up empty array for tcdq:
tcdq = pl.zeros([len(years),filenames.shape[1],1,1,256,512])

# read all the netCDF files to extract tcdq
#loop over years:
for year in range(len(years)):
    #loop over filenames:
    for name in range(filenames.shape[1]):
        #load ncfile
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                str(year_input[year]) + '/' + str(filenames[year,name]),'r')
        #extract E & TP data
        tcdq[year,name] = ncfile.variables['TCDQ'][:]
        ncfile.close()

# calculate the mean tcdq for each year
tcdq_yrs_mn = pl.mean(tcdq,axis=1)

# Get rid of the 1D axes:
tcdq_yrs_mn = pl.squeeze(tcdq_yrs_mn)

# Annual mean tcdq:
tcdq_mn = pl.mean(tcdq_yrs_mn,axis=0)

atl_tcdq = tcdq_yrs_mn*atlmask
atl_tcdq_mn = tcdq_mn*atlmask

pac_tcdq = tcdq_yrs_mn*pacmask
pac_tcdq_mn = tcdq_mn*pacmask

ind_tcdq = tcdq_yrs_mn*indmask
ind_tcdq_mn = tcdq_mn*indmask

#------------------------------------------------------------------------------

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



#-----------------------CALCULATE THE AREAS OF GRID CELLS----------------------

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
#------------------------------------------------------------------------------


# set up empty arrays, size 36, for basin-integrated P & E
oc_evapint = pl.zeros([evap_yrs_mn.shape[0]])
oc_precint = pl.zeros_like(oc_evapint)
atl_evapint = pl.zeros([evap_yrs_mn.shape[0]])
atl_precint = pl.zeros_like(atl_evapint)
atl_evapano = pl.zeros_like(atl_evapint)
atl_precano = pl.zeros_like(atl_evapint)
atl_tcdqint = pl.zeros_like(atl_evapint)
atl_tcdqano = pl.zeros_like(atl_evapint)
pac_evapint = pl.zeros_like(atl_evapint)
pac_precint = pl.zeros_like(atl_evapint)
pac_evapano = pl.zeros_like(atl_evapint)
pac_precano = pl.zeros_like(atl_evapint)
pac_tcdqint = pl.zeros_like(atl_evapint)
pac_tcdqano = pl.zeros_like(atl_evapint)
ind_evapint = pl.zeros_like(atl_evapint)
ind_precint = pl.zeros_like(atl_evapint)
ind_evapano = pl.zeros_like(atl_evapint)
ind_precano = pl.zeros_like(atl_evapint)
ind_tcdqint = pl.zeros_like(atl_evapint)
ind_tcdqano = pl.zeros_like(atl_evapint)
# with a rescaling factor, use NetEmPCalc to get basin-integrated E & P for each year
# No rescaling factor needed for tcdq
rho = 1000 # density of water
factor = rho/86400 # rescaling factor
for yr in range(len(years)):
    oc_evapint[yr] = NetEmPCalc(areas,oc_evap[yr]*factor,rho)
    oc_precint[yr] = NetEmPCalc(areas,oc_prec[yr]*factor,rho)
    atl_evapint[yr] = NetEmPCalc(areas*(1-lsm[0]),atl_evap[yr]*factor,rho)
    atl_precint[yr] = NetEmPCalc(areas*(1-lsm[0]),atl_prec[yr]*factor,rho)
    atl_evapano[yr] = NetEmPCalc(areas*(1-lsm[0]),(atl_evap[yr]-atl_evap_mn)*factor,rho)
    atl_precano[yr] = NetEmPCalc(areas*(1-lsm[0]),(atl_prec[yr]-atl_prec_mn)*factor,rho)
    atl_tcdqint[yr] = NetEmPCalc(areas*(1-lsm[0]),atl_tcdq[yr],rho)
    atl_tcdqano[yr] = NetEmPCalc(areas*(1-lsm[0]),atl_tcdq[yr]-atl_tcdq_mn,rho)
    pac_evapint[yr] = NetEmPCalc(areas*(1-lsm[0]),pac_evap[yr]*factor,rho)
    pac_precint[yr] = NetEmPCalc(areas*(1-lsm[0]),pac_prec[yr]*factor,rho)
    pac_evapano[yr] = NetEmPCalc(areas*(1-lsm[0]),(pac_evap[yr]-pac_evap_mn)*factor,rho)
    pac_precano[yr] = NetEmPCalc(areas*(1-lsm[0]),(pac_prec[yr]-pac_prec_mn)*factor,rho)
    pac_tcdqint[yr] = NetEmPCalc(areas*(1-lsm[0]),pac_tcdq[yr],rho)
    pac_tcdqano[yr] = NetEmPCalc(areas*(1-lsm[0]),pac_tcdq[yr]-pac_tcdq_mn,rho)
    ind_evapint[yr] = NetEmPCalc(areas*(1-lsm[0]),ind_evap[yr]*factor,rho)
    ind_precint[yr] = NetEmPCalc(areas*(1-lsm[0]),ind_prec[yr]*factor,rho)
    ind_evapano[yr] = NetEmPCalc(areas*(1-lsm[0]),(ind_evap[yr]-ind_evap_mn)*factor,rho)
    ind_precano[yr] = NetEmPCalc(areas*(1-lsm[0]),(ind_prec[yr]-ind_prec_mn)*factor,rho)
    ind_tcdqint[yr] = NetEmPCalc(areas*(1-lsm[0]),ind_tcdq[yr],rho)
    ind_tcdqano[yr] = NetEmPCalc(areas*(1-lsm[0]),ind_tcdq[yr]-ind_tcdq_mn,rho)

# E-P for each year
atl_EmPint = atl_evapint - atl_precint
pac_EmPint = pac_evapint - pac_precint
ind_EmPint = ind_evapint - ind_precint
# E-P anomalies:
atl_EmPano = atl_EmPint - pl.mean(atl_EmPint)
pac_EmPano = pac_EmPint - pl.mean(pac_EmPint)
ind_EmPano = ind_EmPint - pl.mean(ind_EmPint)

# calculate standard deviations
atl_evp_sd = pl.std(atl_evapint); print 'Atlantic E standard deviation: ', atl_evp_sd
atl_prc_sd = pl.std(atl_precint); print 'Atlantic P standard deviation: ', atl_prc_sd
atl_EmP_sd = pl.std(atl_EmPint); print 'Atlantic E-P standard deviation: ', atl_EmP_sd
pac_evp_sd = pl.std(pac_evapint); print 'Pacific E standard deviation: ', pac_evp_sd
pac_prc_sd = pl.std(pac_precint); print 'Pacific P standard deviation: ', pac_prc_sd
pac_EmP_sd = pl.std(pac_EmPint); print 'Pacific E-P standard deviation: ', pac_EmP_sd
ind_evp_sd = pl.std(ind_evapint); print 'Indian E standard deviation: ', ind_evp_sd
ind_prc_sd = pl.std(ind_precint); print 'Indian P standard deviation: ', ind_prc_sd
ind_EmP_sd = pl.std(ind_EmPint); print 'Indian E-P standard deviation: ', ind_EmP_sd

# Calculate correlations between E-P and E,P
atl_evp_cr = pearsonr(atl_evapint,atl_EmPint); print 'Atlantic correlation with E & E-P: ', atl_evp_cr[0]
atl_prc_cr = pearsonr(atl_precint,atl_EmPint); print 'Atlantic correlation with P & E-P: ', atl_prc_cr[0]
pac_evp_cr = pearsonr(pac_evapint,pac_EmPint); print 'Pacific correlation with E & E-P: ', pac_evp_cr[0]
pac_prc_cr = pearsonr(pac_precint,pac_EmPint); print 'Pacific correlation with P & E-P: ', pac_prc_cr[0]
ind_evp_cr = pearsonr(ind_evapint,ind_EmPint); print 'Indian correlation with E & E-P: ', ind_evp_cr[0]
ind_prc_cr = pearsonr(ind_precint,ind_EmPint); print 'Indian correlation with P & E-P: ', ind_prc_cr[0]

# plot time series
"""fig, ax  = pl.subplots(3,1,figsize=(10,10))

ax1 = pl.subplot(311)
pl.plot(pl.linspace(0,35,36),atl_evapint,color='r')
pl.plot(pl.linspace(0,35,36),atl_precint,color='b')
pl.plot(pl.linspace(0,35,36),atl_EmPint,color='k')
pl.plot(pl.linspace(0,35,36),atl_tcdqint,color='grey')
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(0,4)
ax1.set_yticks(pl.linspace(0,4,5))
ax1.set_ylabel('Sv',fontsize=16)
pl.title('Atlantic Ocean',fontsize=16)

ax2 = pl.subplot(312)
pl.plot(pl.linspace(0,35,36),pac_evapint,color='r')
pl.plot(pl.linspace(0,35,36),pac_precint,color='b')
pl.plot(pl.linspace(0,35,36),pac_EmPint,color='k')
pl.plot(pl.linspace(0,35,36),pac_tcdqint,color='grey')
pl.axhline(y=0,color='k',ls='--')
ax2.xaxis.set_major_formatter(pl.NullFormatter())
ax2.set_yticks(pl.linspace(-1,9,6))
ax2.set_ylabel('Sv',fontsize=16)
pl.title('Pacific Ocean',fontsize=16)

ax3 = pl.subplot(313)
pl.plot(pl.linspace(0,35,36),ind_evapint,color='r',label='$E$')
pl.plot(pl.linspace(0,35,36),ind_precint,color='b',label='$P$')
pl.plot(pl.linspace(0,35,36),ind_EmPint,color='k',label='$E-P$')
pl.plot(pl.linspace(0,35,36),ind_tcdqint,color='grey',label='$\\nabla\cdot\\rho q u$')
ax3.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=12)
ax3.set_xlabel('Years',fontsize=16)
pl.ylim(0,4)
ax3.set_yticks(pl.linspace(0,4,5))
ax3.set_ylabel('Sv',fontsize=16)
pl.title('Indian Ocean',fontsize=16)

ax3.legend(loc=2,ncol=4)
fig.suptitle('Ocean Surface Freshwater Fluxes',fontsize=18)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/E_TP_interann_basint.png')


# Plot anomalies:

fig, ax = pl.subplots(3,1,figsize=(10,10))

ax1 = pl.subplot(311)
pl.plot(pl.linspace(0,35,36),atl_evapano,color='r',label='$E$')
pl.plot(pl.linspace(0,35,36),-atl_precano,color='b',label='$P$')
pl.plot(pl.linspace(0,35,36),atl_EmPano,color='k',label='$E-P$')
#pl.plot(pl.linspace(0,35,36),atl_tcdqano,color='grey',label='$\\nabla\cdot\\rho q u$')
pl.axhline(y=0,color='k',ls='--')
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(-0.6,0.6)
#ax1.set_yticks(pl.linspace(0,4,5))
ax1.set_ylabel('Sv',fontsize=16)
pl.title('Atlantic Ocean',fontsize=16)
ax1.legend(loc=2,ncol=3)

ax2 = pl.subplot(312)
pl.plot(pl.linspace(0,35,36),pac_evapano,color='r')
pl.plot(pl.linspace(0,35,36),-pac_precano,color='b')
pl.plot(pl.linspace(0,35,36),pac_EmPano,color='k')
#pl.plot(pl.linspace(0,35,36),pac_tcdqano,color='grey')
pl.axhline(y=0,color='k',ls='--')
ax2.xaxis.set_major_formatter(pl.NullFormatter())
#ax2.set_yticks(pl.linspace(-1,9,6))
ax2.set_ylabel('Sv',fontsize=16)
pl.title('Pacific Ocean',fontsize=16)


ax3 = pl.subplot(313)
pl.plot(pl.linspace(0,35,36),ind_evapano,color='r')
pl.plot(pl.linspace(0,35,36),-ind_precano,color='b')
pl.plot(pl.linspace(0,35,36),ind_EmPano,color='k')
#pl.plot(pl.linspace(0,35,36),ind_tcdqano,color='grey')
pl.axhline(y=0,color='k',ls='--')
ax3.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=12)
ax3.set_xlabel('Years',fontsize=16)
pl.ylim(-0.6,0.6)
#ax3.set_yticks(pl.linspace(0,4,5))
ax3.set_ylabel('Sv',fontsize=16)
pl.title('Indian Ocean',fontsize=16)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/E_TP_anom_basint.png')"""


# Calculate residual time series:

atlres_int = atl_tcdqint - atl_EmPint # integrated residual
atlres_ano = atlres_int - pl.mean(atlres_int)

H = -atl_precint + atlres_int
H_ano = H  - pl.mean(H)

#-------------------------SCALE FLUXES BY BASIN AREA---------------------------

scl_fctr = (10**6)*100*86400*365 # factor to rescale from Sv/m^2 to cm/yr

# setup empty arrays for scaled fluxes:
oc_evapscl = pl.zeros_like(oc_evapint)
oc_precscl = pl.zeros_like(oc_precint)
atl_evapscl = pl.zeros_like(atl_evapint)
atl_precscl = pl.zeros_like(atl_evapscl)
atl_evpscl_an = pl.zeros_like(atl_evapscl)
atl_prcscl_an = pl.zeros_like(atl_evapscl)
atl_tcdqscl = pl.zeros_like(atl_evapscl)
atl_tcdqscl_an = pl.zeros_like(atl_evapscl)
pac_evapscl = pl.zeros_like(atl_evapscl)
pac_precscl = pl.zeros_like(atl_evapscl)
pac_evpscl_an = pl.zeros_like(atl_evapscl)
pac_prcscl_an = pl.zeros_like(atl_evapscl)
pac_tcdqscl = pl.zeros_like(pac_evapscl)
pac_tcdqscl_an = pl.zeros_like(pac_evapscl)
ind_evapscl = pl.zeros_like(atl_evapscl)
ind_precscl = pl.zeros_like(atl_evapscl)
ind_evpscl_an = pl.zeros_like(atl_evapscl)
ind_prcscl_an = pl.zeros_like(atl_evapscl)
ind_tcdqscl = pl.zeros_like(ind_evapscl)
ind_tcdqscl_an = pl.zeros_like(ind_evapscl)

atlres_scl = pl.zeros_like(oc_evapscl)
atlres_scl_an = pl.zeros_like(atlres_scl)
H_scl = pl.zeros_like(oc_evapscl)
H_scl_an = pl.zeros_like(oc_evapscl)

for yr in range(len(years)):
    oc_evapscl[yr] = FluxScaled(areas,1-lsm,oc_evapint[yr])*scl_fctr
    oc_precscl[yr] = FluxScaled(areas,1-lsm,oc_precint[yr])*scl_fctr
    atl_evapscl[yr] = FluxScaled(areas*(1-lsm[0]),atlmask,atl_evapint[yr])*scl_fctr
    atl_precscl[yr] = FluxScaled(areas*(1-lsm[0]),atlmask,atl_precint[yr])*scl_fctr
    atl_evpscl_an[yr] = FluxScaled(areas*(1-lsm[0]),atlmask,atl_evapano[yr])*scl_fctr
    atl_prcscl_an[yr] = FluxScaled(areas*(1-lsm[0]),atlmask,atl_precano[yr])*scl_fctr
    atl_tcdqscl[yr] = FluxScaled(areas*(1-lsm[0]),atlmask,atl_tcdqint[yr])*scl_fctr
    atl_tcdqscl_an[yr] = FluxScaled(areas*(1-lsm[0]),atlmask,atl_tcdqano[yr])*scl_fctr
    pac_evapscl[yr] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_evapint[yr])*scl_fctr
    pac_precscl[yr] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_precint[yr])*scl_fctr
    pac_evpscl_an[yr] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_evapano[yr])*scl_fctr
    pac_prcscl_an[yr] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_precano[yr])*scl_fctr
    pac_tcdqscl[yr] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_tcdqint[yr])*scl_fctr
    pac_tcdqscl_an[yr] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_tcdqano[yr])*scl_fctr
    ind_evapscl[yr] = FluxScaled(areas*(1-lsm[0]),indmask,ind_evapint[yr])*scl_fctr
    ind_precscl[yr] = FluxScaled(areas*(1-lsm[0]),indmask,ind_precint[yr])*scl_fctr
    ind_evpscl_an[yr] = FluxScaled(areas*(1-lsm[0]),indmask,ind_evapano[yr])*scl_fctr
    ind_prcscl_an[yr] = FluxScaled(areas*(1-lsm[0]),indmask,ind_precano[yr])*scl_fctr
    ind_tcdqscl[yr] = FluxScaled(areas*(1-lsm[0]),indmask,ind_tcdqint[yr])*scl_fctr
    ind_tcdqscl_an[yr] = FluxScaled(areas*(1-lsm[0]),indmask,ind_tcdqano[yr])*scl_fctr
    
    atlres_scl[yr] = FluxScaled(areas,atlmask,atlres_int[yr])*scl_fctr
    atlres_scl_an[yr] = FluxScaled(areas,atlmask,atlres_ano[yr])*scl_fctr
    H_scl[yr] = FluxScaled(areas,atlmask,H[yr])*scl_fctr
    H_scl_an[yr] = FluxScaled(areas,atlmask,H_ano[yr])*scl_fctr

# E-P scaled by area:
atl_EmPscl = atl_evapscl-atl_precscl
pac_EmPscl = pac_evapscl-pac_precscl
ind_EmPscl = ind_evapscl-ind_precscl
# Scaled E-P anomalies:
atl_EmPscl_an = atl_EmPscl - pl.mean(atl_EmPscl)
pac_EmPscl_an = pac_EmPscl - pl.mean(pac_EmPscl)
ind_EmPscl_an = ind_EmPscl - pl.mean(ind_EmPscl)




"""fig, ax = pl.subplots(3,1,figsize=(10,10))

ax1 = pl.subplot(311)
pl.plot(pl.linspace(0,35,36),atl_evapscl,color='r')
pl.plot(pl.linspace(0,35,36),atl_precscl,color='b')
pl.plot(pl.linspace(0,35,36),atl_EmPscl,color='k')
pl.plot(pl.linspace(0,35,36),atl_tcdqscl,color='grey')
pl.axhline(y=0,ls='--',color='k')
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(-10,170)
ax1.set_yticks(pl.linspace(-10,170,7))
ax1.set_ylabel('cm/yr',fontsize=16)
pl.title('Atlantic Ocean',fontsize=16)

ax2 = pl.subplot(312)
pl.plot(pl.linspace(0,35,36),pac_evapscl,color='r',label='$E$')
pl.plot(pl.linspace(0,35,36),pac_precscl,color='b',label='$P$')
pl.plot(pl.linspace(0,35,36),pac_EmPscl,color='k',label='$E-P$')
pl.plot(pl.linspace(0,35,36),pac_tcdqscl,color='grey',label='$\\nabla\cdot\\rho q u$')
pl.axhline(y=0,ls='--',color='k')
ax2.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(-10,170)
ax2.set_yticks(pl.linspace(-10,170,7))
ax2.set_ylabel('cm/yr',fontsize=16)
pl.title('Pacific Ocean',fontsize=16)

ax3 = pl.subplot(313)
pl.plot(pl.linspace(0,35,36),ind_evapscl,color='r',)
pl.plot(pl.linspace(0,35,36),ind_precscl,color='b')
pl.plot(pl.linspace(0,35,36),ind_EmPscl,color='k')
pl.plot(pl.linspace(0,35,36),ind_tcdqscl,color='grey')
pl.axhline(y=0,ls='--',color='k')
ax3.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=12)
ax3.set_xlabel('Years',fontsize=16)
pl.ylim(-10,170)
ax3.set_yticks(pl.linspace(-10,170,7))
ax3.set_ylabel('cm/yr',fontsize=16)
pl.title('Indian Ocean',fontsize=16)

ax2.legend(loc=6)
fig.suptitle('Ocean Surface Freshwater Fluxes scaled by area',fontsize=18)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/E_TP_interann_basscl.png')"""


# Read in GPCP precipitation annual means:

GPCP_means = pl.genfromtxt('/home/np838619/GPCP/GPCP_annmeans.txt',skip_header=3)
gpcp_atl = GPCP_means[:,0]
gpcp_pac = GPCP_means[:,1]
gpcp_ind = GPCP_means[:,2]

# Read in GPCP precipitation anomalies:

GPCP_anoms = pl.genfromtxt('/home/np838619/GPCP/GPCP_anoms.txt',skip_header=3)
gpcp_atl_an = GPCP_anoms[:,0]
gpcp_pac_an = GPCP_anoms[:,1]
gpcp_ind_an = GPCP_anoms[:,2]

#residuals = pl.genfromtxt('/home/np838619/PminusE_data/ERA_Int/res_basins_yrs.txt',skip_header=3)
#atlres = residuals[:,0]
#pacres = residuals[:,1]
#indres = residuals[:,2]

#atlres_scl = atl_tcdqscl - atl_EmPscl
pacres_scl = pac_tcdqscl - pac_EmPscl
indres_scl = ind_tcdqscl - ind_EmPscl

#z1 = atlres_scl - pl.mean(atlres_scl)
#z2 = -atl_prcscl_an + atlres_scl_an

v1 = pacres_scl - pl.mean(pacres_scl)
v2 = -pac_prcscl_an + v1

b1 = indres_scl - pl.mean(indres_scl)
b2 = -ind_prcscl_an + b1

# Plot anomalies:

fig, ax = pl.subplots(3,1,figsize=(10,9))

ax1 = pl.subplot(311)
pl.plot(pl.linspace(0,35,36),atl_evpscl_an,color='r',label='$\\bar{e}$')
pl.plot(pl.linspace(0,35,36),-atl_prcscl_an,color='b',label='$-\\bar{p}$')
pl.plot(pl.linspace(0,35,36),-gpcp_atl_an,color='b',ls='--',label='$-$GPCP')
pl.plot(pl.linspace(0,35,36),atl_EmPscl_an,color='k',label='$\\bar{e}-\\bar{p}$')
#pl.plot(pl.linspace(0,35,36),H_scl_an,color='g',marker='+',markersize=12.,ls='None',label='$-P+ $res')
pl.plot(pl.linspace(0,35,36),atl_tcdqscl_an,color='grey',label='div$\mathbf{Q}$',
        linestyle=':',linewidth=3)
#pl.axhline(y=0,ls='--',color='k')
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(-15,15); pl.grid(axis='y')
#ax1.set_yticks(pl.linspace(-10,170,7))
ax1.set_ylabel('cm/yr',fontsize=20,labelpad=-10)
#ax1.annotate('(a)',(0,1.02),xycoords='axes fraction',size=25)
pl.setp(ax1.get_yticklabels(),fontsize=16)
pl.title('(a) Atlantic',fontsize=20)

LINES = [Line2D([0],[0],color='r'),Line2D([0],[0],color='b'),
         Line2D([0],[0],color='b',ls='--'),Line2D([0],[0],color='k'),
         Line2D([0],[0],color='grey',ls=':')]
LABS = ['$\\bar{e}$','$\\bar{p}$','$-$GPCP','$\\bar{e}-\\bar{p}$','div$\mathbf{Q}$']
leg = ax1.legend(loc=2,ncol=3,framealpha=1.,frameon=False,fontsize=20,
           handletextpad=0.3,columnspacing=0.5)
for legobj in leg.legendHandles:
    legobj.set_linewidth(3.0)

ax2 = pl.subplot(312)
pl.plot(pl.linspace(0,35,36),pac_evpscl_an,color='r',label='$\\bar{e}$')
pl.plot(pl.linspace(0,35,36),-pac_prcscl_an,color='b',label='$-\\bar{p}$')
pl.plot(pl.linspace(0,35,36),-gpcp_pac_an,color='b',ls='--',label='$-$GPCP')
pl.plot(pl.linspace(0,35,36),pac_EmPscl_an,color='k',label='$\\bar{e}-\\bar{p}$')
#pl.plot(pl.linspace(0,35,36),v2,color='g',marker='+',markersize=12.,ls='None')
pl.plot(pl.linspace(0,35,36),pac_tcdqscl_an,color='grey',label='div$\mathbf{Q}$',
        linestyle=':',linewidth=3)
#pl.axhline(y=0,ls='--',color='k')
ax2.xaxis.set_major_formatter(pl.NullFormatter())
#ax2.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=18)
#ax2.set_xlabel('Years',fontsize=22)
#ax2.annotate('(b)',(0,1.02),xycoords='axes fraction',size=22)
pl.ylim(-15,15); pl.grid(axis='y')
#ax2.set_yticks(pl.linspace(-10,170,7))
ax2.set_ylabel('cm/yr',fontsize=20,labelpad=-10)
pl.setp(ax2.get_yticklabels(),fontsize=16)
#ax2.tick_params(axis='x', which='major', pad=10)
pl.title('(b) Pacific',fontsize=20)
#ax2.legend(loc=2,ncol=3,framealpha=1.,frameon=False,fontsize=20,
#           handletextpad=0.3,columnspacing=0.5)
#pl.subplots_adjust(top=0.94,bottom=0.15)

ax3 = pl.subplot(313)
pl.plot(pl.linspace(0,35,36),ind_evpscl_an,color='r')
pl.plot(pl.linspace(0,35,36),-ind_prcscl_an,color='b')
pl.plot(pl.linspace(0,35,36),-gpcp_ind_an,color='b',ls='--')
pl.plot(pl.linspace(0,35,36),ind_EmPscl_an,color='k')
#pl.plot(pl.linspace(0,35,36),b2,color='g',marker='+',markersize=12.,ls='None')
pl.plot(pl.linspace(0,35,36),ind_tcdqscl_an,color='grey',linestyle=':',linewidth=3)
#pl.axhline(y=0,ls='--',color='k')
ax3.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=12)
ax3.set_xlabel('Years',fontsize=20)
#ax3.annotate('(c)',(0,1.02),xycoords='axes fraction',size=25)
pl.setp(ax3.get_xticklabels(),fontsize=16)
#pl.rcParams['xtick.major.pad']='16'
pl.ylim(-15,15); pl.grid(axis='y')
#ax3.set_yticks(pl.linspace(-10,170,7))
ax3.set_ylabel('cm/yr',fontsize=20,labelpad=-10)
pl.setp(ax3.get_yticklabels(),fontsize=18)
ax3.tick_params(axis='x', which='major', pad=10)
pl.title('(c) Indian',fontsize=20)
pl.tight_layout()
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/E_TP_anom_basscl_GPCP.png')
#pl.close()

#ax, fig = pl.subplots()
#pl.plot(pl.linspace(0,35,36),atl_tcdqscl,color='grey',label='tcdq')
#pl.plot(pl.linspace(0,35,36),atl_EmPscl,color='k',label='$E-P$')
#pl.plot(pl.linspace(0,35,36),atl_precscl,color='b',label='$P$')
#pl.plot(pl.linspace(0,35,36),gpcp_atl,color='b',ls='--',label='GPCP')
#pl.plot(pl.linspace(0,35,36),atl_precscl-atlres_scl,color='g',label='$P-$ res')
#pl.title('Atlantic')
#pl.matplotlib.pyplot.legend(loc=6)
#
#ax, fig = pl.subplots()
#pl.plot(pl.linspace(0,35,36),pac_tcdqscl,color='grey',label='tcdq')
#pl.plot(pl.linspace(0,35,36),pac_EmPscl,color='k',label='$E-P$')
#pl.plot(pl.linspace(0,35,36),pac_precscl,color='b',label='$P$')
#pl.plot(pl.linspace(0,35,36),gpcp_pac,color='b',ls='--',label='GPCP')
#pl.plot(pl.linspace(0,35,36),pac_precscl-pacres_scl,color='g',label='$P-$ res')
#pl.title('Pacific')
#pl.matplotlib.pyplot.legend(loc=6)
#
#ax, fig = pl.subplots()
#pl.plot(pl.linspace(0,35,36),ind_tcdqscl,color='grey',label='tcdq')
#pl.plot(pl.linspace(0,35,36),ind_EmPscl,color='k',label='$E-P$')
#pl.plot(pl.linspace(0,35,36),ind_precscl,color='b',label='$P$')
#pl.plot(pl.linspace(0,35,36),gpcp_ind,color='b',ls='--',label='GPCP')
#pl.plot(pl.linspace(0,35,36),ind_precscl-indres_scl,color='g',label='$P-$ res')
#pl.title('Indian')
#pl.matplotlib.pyplot.legend(loc=6)