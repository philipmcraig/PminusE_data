# -*- coding: utf-8 -*-
"""
Created on Wed Sep 09 4:27:00 2015

Code for calculating the areas of ERA-Interim grid cells on the half grid and
multiplying cell areas by the E-P value to determine the total E-P, in Sverdrups,
for ocean basins (& latitude bands) and scaling everything by area.

@author: np838619

Last updated: 10/9/2015 4:26PM 1st commit
"""


from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import pandas
from functions import *

def FluxScaled(areas,mask,flux):
    basin_areas = areas*mask
    flux_scaled = flux/pl.sum(basin_areas)
    
    return flux_scaled

def FluxBandScaled(areas,mask,indices,flux_bands):
    basin_areas = areas*mask
    flux_scaled = pl.zeros_like(flux_bands)
    for i in range(flux_bands.shape[0]):
        indpair = (indices[i],indices[i+1])
        band_area = pl.sum(basin_areas[indpair[0]:indpair[1]+1])
        flux_scaled[i] = flux_bands[i]/band_area
        
    return flux_scaled

#----open up netcdf files for basins here & extract lon, lat & tcdq arrays-----
ncfile = Dataset('/home/np838619/PminusE_data/ERA_Int/tcdq_basins_60.nc','r') # Gaussian grid
lon = ncfile.variables['lon'][:] # extract lon
lat = ncfile.variables['lat'][:] # extract lat
tcdq_atl = ncfile.variables['tcdq_Atl'][:]  # extract Atlantic tcdq
atlmask = ncfile.variables['maskA'][:]
tcdq_bal = ncfile.variables['tcdq_Bal'][:] # extract Baltic tcdq
balmask = ncfile.variables['maskB'][:]
#tcdq_med = ncfile.variables['tcdq_Med'][:]
medmask = ncfile.variables['maskM'][:]
tcdq_pac = ncfile.variables['tcdq_Pac'][:]
pacmask = ncfile.variables['maskP'][:]
tcdq_ind = ncfile.variables['tcdq_Ind'][:]
indmask = ncfile.variables['maskI'][:]
ncfile.close() # close nc file

nc2 = Dataset('/home/np838619/PminusE_data/ERA_Int/tcdq_36years_means.nc','r')
tcdq = nc2.variables['tcdq'][:]
nc2.close()

tcdq = pl.mean(tcdq[:],axis=0)

# Need to covert nan to 0 so that the arrays add together properly
tcdq_atl = NaN2Zero(tcdq_atl)
tcdq_bal = NaN2Zero(tcdq_bal)
tcdq_med = tcdq*medmask
tcdq_AB = tcdq_atl + tcdq_bal + tcdq_med

tcdq_pac = NaN2Zero(tcdq_pac)
tcdq_ind = NaN2Zero(tcdq_ind)

# add the Atlantic & Baltic masks together:
ABmask = atlmask + balmask + medmask
tcdq_A = tcdq*ABmask


# -----------------Read in Dai & Trenberth runoff csv file here----------------
runoff = pandas.read_csv('DT03_runoff.csv',sep=',',header=1,skiprows=0,
    names=['row','lat','Atl','Pac','Ind','Glo'])
# to get the values of a column type runoff.name e.g. runoff.Pac
run_atl = runoff.Atl
run_pac = runoff.Pac
run_ind = runoff.Ind

# which latitudes will be required for splitting the basins into subregions?
atl_lats = pl.array([60,45,24,-16,-35])
pac_lats = pl.array([65,47,24,-30])
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


#--------------calculate cell areas here, use function from above--------------
# set up empty array for area, size lat_half X lon_half
areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
radius = 6.37*(10**6)
# loop over latitude and longitude
for i in range(lat_half.shape[0]-1): # loops over 256
    latpair = (lat_half[i],lat_half[i+1])
    for j in range(lon.shape[0]): # loops over 512
        #lonpair = (lon_half[i],lon_half[i+1])
        areas[i,j] = AreaFullGaussianGrid(radius,delta_lambda,latpair)
        #areas[i,j] = AreaCalc2(radius,delta_lambda,lat_half[i],delta_phi)
#print "%.4g" % pl.sum(areas)
        

#-calculate the total E-P for each cell down here and add everything up to get- 
#-------------------total E-P for the basin in Sverdrups-----------------------
#--------------------see NetEmPCalc function for details-----------------------
rho = 10**3 # density of water in kg/m^3
# calculate total E-P for each basin:
A_total = NetEmPCalc(areas,tcdq_A,rho) # Atlantic (35S-45N) + Med/Black seas
M_total = NetEmPCalc(areas,tcdq_med,rho)
P_total = NetEmPCalc(areas,tcdq_pac,rho) # Pacififc (30S-47N)
I_total = NetEmPCalc(areas,tcdq_ind,rho) # Indian (>35S)
#print AM_total, P_total, I_total
EmP = pl.array([A_total,P_total,I_total])

M_EmP_scl = FluxScaled(areas,medmask,M_total)

#------------------------------ADD RUNOFF-------------------------------------
# Calculate the total runoff for each ocean basin from the Dai & Trenberth dataset
#---------------------see TotalRunoff function for details---------------------

# use the indices pf the northern & southern boundaries for the relevant region of each ocean basin 
atl_lim = (atl_ind[0],atl_ind[-1])
pac_lim = (pac_ind[0],pac_ind[-1])
ind_lim = (ind_ind[0],ind_ind[-1])
A_EPR = A_total - TotalRunoff(runoff.lat,lat,atl_lim,run_atl) # Atlantic (35S-45N) + Med/Black seas
P_EPR = P_total - TotalRunoff(runoff.lat,lat,pac_lim,run_pac) # Pacififc (30S-47N)
I_EPR = I_total - TotalRunoff(runoff.lat,lat,ind_lim,run_ind) # Indian (>35S)
#print Atl_EPR, P_EPR, I_EPR
EPR = pl.array([A_EPR,P_EPR,I_EPR])

# Scale E-P-R by area:

A_EPR_scl = FluxScaled(areas,ABmask,A_EPR)
P_EPR_scl = FluxScaled(areas,pacmask,P_EPR)
I_EPR_scl = FluxScaled(areas,indmask,I_EPR)

# rescale to change from Sv/m^2 to cm/yr:
factor = (10**6)*100*86400*365
print 'Atlantic Ocean = ', A_EPR_scl*factor, ' cm/yr'
print 'Pacific Ocean = ', P_EPR_scl*factor, ' cm/yr'
print 'Indian Ocean = ', I_EPR_scl*factor, ' cm/yr'


#------------------------CALCULATE E-P FOR OCEAN SUBREGIONS--------------------
#----------------------see NetEmPBands function for details--------------------
# From literature the bands are:
    # Atlantic bands: 24N-45N, 16S-24N, 35S-16S
    # Pacific bands: 24N-47N, 30S-24N
    # Indian bands: >8S, 20S-8S, 35S-20S

A_bands = NetEmPBands(atl_ind,areas,tcdq_AB,rho) 
print A_bands

P_bands = NetEmPBands(pac_ind,areas,tcdq_pac,rho)
print P_bands

I_bands = NetEmPBands(ind_ind,areas,tcdq_ind,rho)
print I_bands


#----------------------CALCULATE RUNOFF FOR SUBREGIONS-------------------------
# see EPRbands function for details

A_bands_EPR = EPRbands(A_bands,atl_ind,runoff.lat,lat,run_atl)
print A_bands_EPR

P_bands_EPR = EPRbands(P_bands,pac_ind,runoff.lat,lat,run_pac)
print P_bands_EPR

I_bands_EPR = EPRbands(I_bands,ind_ind,runoff.lat,lat,run_ind)
print I_bands_EPR

# Scale E-P-R by area:
Abnds_scl = FluxBandScaled(areas,ABmask,atl_ind,A_bands_EPR)
Pbnds_scl = FluxBandScaled(areas,pacmask,pac_ind,P_bands_EPR)
Ibnds_scl = FluxBandScaled(areas,indmask,ind_ind,I_bands_EPR)

# rescale to change from Sv/m^2 to cm/yr:

print 'Atlantic Ocean bands = ' ,Abnds_scl*factor, ' cm/yr'
print 'Pacific Ocean bands = ', Pbnds_scl*factor, ' cm/yr'
print 'Indian Ocean bands = ', Ibnds_scl*factor, ' cm/yr'

#---------------------------WRITE RESULTS TO FILE------------------------------

"""f = open('ERAI_FWflux60.txt','w')
f.write('Results from ERA-Interim vertically integrated moisture flux divergence' +
' and Dai & Trenberth (2002) continental discharge. Presented here are net E-P-R' +
' for ocean basins and latitude bands, and net E-P for ocean basins and latitude bands. \n\n' + 
'Net E-P-R for ocean basins (Sv): \nAtlantic (35S-60N)  Pacific (30S-BS)   Indian (>35S) \n')
EPR.tofile(f,sep="   ",format='%.3f')
f.write('\n \nNet E-P-R for latitude bands (Sv): \n \nAtlantic: 45N-60N     24N-45N     16S-24N     35S-16S \n')
A_bands_EPR.tofile(f,sep="   ",format='%.3f')
f.write('\n\nPacific: 47N-BS    24N-47N     30S-24N \n')
P_bands_EPR.tofile(f,sep="   ",format='%.3f')
f.write('\n\nIndian: >8S    20S-8S  35S-20S \n')
I_bands_EPR.tofile(f,sep="  ",format='%.3f')
f.write('\n \nNet E-P for ocean basins (Sv): \nAtlantic (35S-60N)  Pacific (30S-47N)   Indian (>35S) \n')
EmP.tofile(f,sep="   ",format='%.3f')
f.write('\n\nNet E-P for latitude bands (Sv): \n \nAtlantic: 45N-60N      24N-45N     16S-24N     35S-16S \n')
A_bands.tofile(f,sep="   ",format='%.3f')
f.write('\n\nPacific: 47N-BS   24N-47N   30S-24N \n')
P_bands.tofile(f,sep="   ",format='%.3f')
f.write('\n\nIndian: >8S    20S-8S  35S-20S \n')
I_bands.tofile(f,sep=" ",format='%.3f')
f.close()"""