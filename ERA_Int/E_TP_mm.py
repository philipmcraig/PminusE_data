# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 17:46:12 2016

Code to calculate monthly means of 1979-2014 basin-integrated evaporation and 
precipitation from accumulated surface forecasts. Also want to calculate a
monthly moisture budget residual.

Last updated: 27/01/2016 2:06AM 2nd commit

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from scipy.stats import pearsonr, linregress
from scipy.signal import detrend
from functions import *

def ScatterPlot(gpcp,eraint):
    """
    """
    r = linregress(gpcp,eraint.flatten())
    x = pl.linspace(0,10,11); y = r[0]*x + r[1]
    pl.scatter(gpcp,eraint.flatten())
    pl.plot(x,y,label='regression line')
    pl.plot(x,x,ls='--',color='k',label='1:1 line')
    pl.xlim(0,5); pl.ylim(0,5)
    pl.xlabel('GPCP (mm/day)',fontsize=18)
    pl.text(3.75,0.25,'$r=$'+str(round(r[2],3)),fontsize=15)


def RMSE(obs,model):
    """
    """
    rmse = pl.sqrt(pl.mean((obs-model)**2))
    
    return rmse
    
    

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
atlmask = atlmask + balmask + medmask

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

# remove 1-d axes:
evap_tot = pl.squeeze(evap_tot)
prec_tot = pl.squeeze(prec_tot)

EmP_mnths = evap_tot - prec_tot

#-----------------------EXTRACT THE RELEVANT OCEAN BASINS----------------------

atl_evap = evap_tot*atlmask
atl_prec = prec_tot*atlmask
#atl_evap_mn = evap_mn*atlmask
#atl_prec_mn = prec_mn*atlmask

pac_evap = evap_tot*pacmask
pac_prec = prec_tot*pacmask
#pac_evap_mn = evap_mn*pacmask
#pac_prec_mn = prec_mn*pacmask

ind_evap = evap_tot*indmask
ind_prec = prec_tot*indmask
#ind_evap_mn = evap_mn*indmask
#ind_prec_mn = prec_mn*indmask
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
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



# Get rid of the 1D axes:
tcdq = pl.squeeze(tcdq)

atl_tcdq = tcdq*atlmask
pac_tcdq = tcdq*pacmask
ind_tcdq = tcdq*indmask

#------------------------------------------------------------------------------

res_mnths = tcdq - EmP_mnths*(1000/86400) # residual in kg/m2s

atlres = res_mnths*atlmask
pacres = res_mnths*pacmask
indres = res_mnths*indmask

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


# set up empty arrays, size 36x12, for basin-integrated P & E

atl_evapint = pl.zeros([evap_tot.shape[0],evap_tot.shape[1]])
atl_precint = pl.zeros_like(atl_evapint)

pac_evapint = pl.zeros_like(atl_evapint)
pac_precint = pl.zeros_like(atl_evapint)

ind_evapint = pl.zeros_like(atl_evapint)
ind_precint = pl.zeros_like(atl_evapint)


# with a rescaling factor, use NetEmPCalc to get basin-integrated E & P for each month
# No rescaling factor needed for tcdq
rho = 1000 # density of water
factor = rho/86400 # rescaling factor
for yr in range(evap_tot.shape[0]):
    for mt in range(evap_tot.shape[1]):
        atl_evapint[yr,mt] = NetEmPCalc(areas*atlmask*(1-lsm[0]),atl_evap[yr,mt]*factor,rho)
        atl_precint[yr,mt] = NetEmPCalc(areas*atlmask*(1-lsm[0]),atl_prec[yr,mt]*factor,rho)
        pac_evapint[yr,mt] = NetEmPCalc(areas*pacmask*(1-lsm[0]),pac_evap[yr,mt]*factor,rho)
        pac_precint[yr,mt] = NetEmPCalc(areas*pacmask*(1-lsm[0]),pac_prec[yr,mt]*factor,rho)
        ind_evapint[yr,mt] = NetEmPCalc(areas*indmask*(1-lsm[0]),ind_evap[yr,mt]*factor,rho)
        ind_precint[yr,mt] = NetEmPCalc(areas*indmask*(1-lsm[0]),ind_prec[yr,mt]*factor,rho)


#-------------------------SCALE FLUXES BY BASIN AREA---------------------------

scl_fctr = (10**6)*1000*86400 # factor to rescale from Sv/m^2 to mm/day

# setup empty arrays for scaled fluxes:
atl_evapscl = pl.zeros_like(atl_evapint)
atl_precscl = pl.zeros_like(atl_evapscl)

pac_evapscl = pl.zeros_like(atl_evapscl)
pac_precscl = pl.zeros_like(atl_evapscl)

ind_evapscl = pl.zeros_like(atl_evapscl)
ind_precscl = pl.zeros_like(atl_evapscl)


for yr in range(evap_tot.shape[0]):
    for mt in range(evap_tot.shape[1]):
        atl_evapscl[yr,mt] = FluxScaled(areas*(1-lsm[0]),atlmask,atl_evapint[yr,mt])*scl_fctr
        atl_precscl[yr,mt] = FluxScaled(areas*(1-lsm[0]),atlmask,atl_precint[yr,mt])*scl_fctr
        pac_evapscl[yr,mt] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_evapint[yr,mt])*scl_fctr
        pac_precscl[yr,mt] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_precint[yr,mt])*scl_fctr
        ind_evapscl[yr,mt] = FluxScaled(areas*(1-lsm[0]),indmask,ind_evapint[yr,mt])*scl_fctr
        ind_precscl[yr,mt] = FluxScaled(areas*(1-lsm[0]),indmask,ind_precint[yr,mt])*scl_fctr
#------------------------------------------------------------------------------


# Read the GPCP precip from .txt file:

GPCP = pl.genfromtxt('/home/np838619/GPCP/GPCP_mm.txt')
gpcp_atl = GPCP[:,0]
gpcp_pac = GPCP[:,1]
gpcp_ind = GPCP[:,2]


fig, ax  = pl.subplots(1,3,figsize=(15,5))
ax1 = pl.subplot(131)
ScatterPlot(gpcp_atl,atl_precscl)
ax1.set_ylabel('ERA-Interim $P$ (mm/day)',fontsize=18)
pl.title('Atlantic',fontsize=22)
ax1.legend(loc=0)
ax2 = pl.subplot(132)
ScatterPlot(gpcp_pac,pac_precscl)
pl.title('Pacific',fontsize=22)
ax3 = pl.subplot(133)
ScatterPlot(gpcp_ind,ind_precscl)
pl.title('Indian',fontsize=22)
pl.subplots_adjust(left=0.04,right=0.96,bottom=0.11)


a = linregress(detrend(gpcp_atl),detrend(atl_precscl.flatten()))
b = linregress(detrend(gpcp_pac),detrend(pac_precscl.flatten()))
c = linregress(detrend(gpcp_ind),detrend(ind_precscl.flatten()))

f = RMSE(gpcp_atl,atl_precscl.flatten())
g = RMSE(gpcp_pac,pac_precscl.flatten())
h = RMSE(gpcp_ind,ind_precscl.flatten())


#-------------------------CALCULATE MONTHLY RESIDUAL---------------------------

atl_resint = pl.zeros([atlres.shape[0],atlres.shape[1]])
pac_resint = pl.zeros_like(atl_resint)
ind_resint = pl.zeros_like(atl_resint)

rho = 1000 # density of water
# no rescaling factor needed
for yr in range(res_mnths.shape[0]):
    for mt in range(res_mnths.shape[1]):
        atl_resint[yr,mt] = NetEmPCalc(areas*atlmask*(1-lsm[0]),atlres[yr,mt],rho) # results in Sv
        pac_resint[yr,mt] = NetEmPCalc(areas*pacmask*(1-lsm[0]),pacres[yr,mt],rho)
        ind_resint[yr,mt] = NetEmPCalc(areas*indmask*(1-lsm[0]),indres[yr,mt],rho)


# Scale by area:

scl_fctr = (10**6)*100*86400*365 # factor to rescale from Sv/m^2 to cm/yr

atl_resscl = pl.zeros_like(atl_resint)
pac_resscl = pl.zeros_like(atl_resscl)
ind_resscl = pl.zeros_like(atl_resscl)

for yr in range(res_mnths.shape[0]):
    for mt in range(res_mnths.shape[1]):
        atl_resscl[yr,mt] = FluxScaled(areas*(1-lsm[0]),atlmask,atl_resint[yr,mt])*scl_fctr
        pac_resscl[yr,mt] = FluxScaled(areas*(1-lsm[0]),pacmask,pac_resint[yr,mt])*scl_fctr
        ind_resscl[yr,mt] = FluxScaled(areas*(1-lsm[0]),indmask,ind_resint[yr,mt])*scl_fctr