# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 11:49:25 2015

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import pandas
from scipy.stats import pearsonr, spearmanr
from functions import *

# read in the basin masks from netCDF file
basinfile = Dataset('tcdq_basins_60.nc','r')
atlmask = basinfile.variables['maskA'][:] # mask for Atlantic Ocean
balmask = basinfile.variables['maskB'][:] # mask for Baltic Sea
pacmask = basinfile.variables['maskP'][:] # mask for Pacific Ocean
indmask = basinfile.variables['maskI'][:] # mask for Indian Ocean
lon = basinfile.variables['lon'][:] # ERA-Interim longitude
lat = basinfile.variables['lat'][:] # ERA-Interim latitude
basinfile.close()

# add the Atlantic & Baltic masks together:
atlmask = atlmask + balmask

# list of years
years = pl.linspace(1979,2014,36)
# need next part to get rid of decimal point
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


#romve 1-dimensional axes:
evap_tot = pl.squeeze(evap_tot)
prec_tot = pl.squeeze(prec_tot)

atl_evap = evap_tot*atlmask
atl_prec = prec_tot*atlmask
pac_evap = evap_tot*pacmask
pac_prec = prec_tot*pacmask


#------------------------------------------------------------------------------

# which latitudes will be required for splitting the basins into subregions?
atl_lats = pl.array([60,35,15,-15,-35])
pac_lats = pl.array([65.5,35,10,-10,-30])
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


atl_evapint = pl.zeros([atl_evap.shape[0],atl_evap.shape[1]])
atl_precint = pl.zeros_like(atl_evapint)
pac_evapint = pl.zeros_like(atl_evapint)
pac_precint = pl.zeros_like(pac_evapint)

# with a rescaling factor, use NetEmPCalc to get basin-integrated E & P for each year
# No rescaling factor needed for tcdq
rho = 1000 # density of water
factor = rho/86400 # rescaling factor
for yr in range(len(years)):
    for mt in range(0,12):
        atl_evapint[yr,mt] = NetEmPCalc(areas,atl_evap[yr,mt]*factor,rho)
        atl_precint[yr,mt] = NetEmPCalc(areas,atl_prec[yr,mt]*factor,rho)
        pac_evapint[yr,mt] = NetEmPCalc(areas,pac_evap[yr,mt]*factor,rho)
        pac_precint[yr,mt] = NetEmPCalc(areas,pac_prec[yr,mt]*factor,rho)


atl_evp_bnds = pl.zeros([36,12,pac_ind.shape[0]-1])
atl_prc_bnds = pl.zeros_like(atl_evp_bnds)
pac_evp_bnds = pl.zeros_like(atl_evp_bnds)
pac_prc_bnds = pl.zeros_like(atl_evp_bnds)


for yr in range(pac_evp_bnds.shape[0]):
    for mt in range(pac_evp_bnds.shape[1]):
        atl_evp_bnds[yr,mt] = NetEmPBands(atl_ind,areas[:,:],atl_evap[yr,mt,:,:],rho)*factor
        atl_prc_bnds[yr,mt] = NetEmPBands(atl_ind,areas[:,:],atl_prec[yr,mt,:,:],rho)*factor
        pac_evp_bnds[yr,mt] = NetEmPBands(pac_ind,areas[:,:],pac_evap[yr,mt,:,:],rho)*factor
        pac_prc_bnds[yr,mt] = NetEmPBands(pac_ind,areas[:,:],pac_prec[yr,mt,:,:],rho)*factor


#--------------------READ IN CSV FILES FOR CLIMATE INDICES---------------------

# OCEANIC NINO INDEX
    # Nino 3.4 region is longitude indices 270-340, P correlation is good inside there

enso = pandas.read_csv('oceanic_nino_index.csv',sep=',',header=0,skiprows=0,
    names=year_input)

z = pl.zeros([36,12])
for yr in range(len(year_input)):
    p = amo[[year_input[yr]]]
    for mt in range(12):
        z[yr,mt] = p.loc[mt]
pac_evap_enso = pearsonr(z.flatten(),pac_evapint.flatten())
pac_prec_enso = pearsonr(z.flatten(),pac_precint.flatten())

v = pac_evp_bnds[:,:,2]
b = pac_prc_bnds[:,:,2]

vm = pl.mean(v); bm = pl.mean(b)

va = v-vm; ba = b-bm
print pearsonr(z.flatten(),va.flatten())
print pearsonr(z.flatten(),ba.flatten())
print pearsonr(z.flatten(),(v-b).flatten())


x = pl.zeros_like(z)