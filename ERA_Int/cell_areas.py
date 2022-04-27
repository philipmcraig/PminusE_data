# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 11:29:37 2015

Code for calculating the areas of ERA-Interim grid cells on the half grid and
multiplying cell areas by the E-P value to determine the total E-P, in Sverdrups,
for ocean basins (& latitude bands)

@author: np838619
"""

#---------------------------import modules here--------------------------------
# pylab, netCDF4 should be all I need!
from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import pandas


def AreaFullGaussianGrid(radius,delta_lambda,latpair):
    """ Function for calculating the area of a grid cell on the surface of a sphere
    where the grid is a full Gaussian Grid, e.g.:
    https://badc.nerc.ac.uk/help/coordinates/cell-surf-area.html
    
    Args:
        radius (float): predefined radius of Earth, 6.37*(10**6) meters
        delta_lambda (float): longitudinal distance between grid points in radians,
                            for a full Gaussian grid delta_lambda is constant
        latpair (tuple): pair of latitude points in radians for calculating 
                        delta_mu, the difference of the sines of each latitude
        
    Returns:
        area (float): area of a grid cell between 2 latitudes in m**2
    
    Example:
        >>> for i in range(latitude):
        ...    for j in range(longitude):
        ...        areas[i,j] = AreaFullGaussianGrid(radius,deta_lambda,latpair[i])
    """

    delta_mu = pl.sin(latpair[0])-pl.sin(latpair[1])
    
    area = (radius**2)*delta_lambda*delta_mu

    return area
    

def AreaCalc2(radius,delta_lambda,phi,delta_phi):
    # this function is for a regular grid
    area = (radius**2)*pl.cos(phi)*delta_lambda*delta_phi
    
    return area

#-----open up netcdf file for basins here & extract lon, lat & tcdq arrays-----
ncfile = Dataset('/home/np838619/PminusE_data/ERA_Int/tcdq_basins2.nc','r') # Gaussian grid
lon = ncfile.variables['lon'][:] # extract lon
lat = ncfile.variables['lat'][:] # extract lat
tcdq_atl = ncfile.variables['tcdq_Atl'][:]  # extract Atlantic tcdq
tcdq_med = ncfile.variables['tcdq_MedBS'][:] # extract Med+BS tcdq
ncfile.close() # close nc file

ncfile2 = Dataset('/home/np838619/PminusE_data/ERA_Int/tcdq_basins3.nc','r') 
tcdq_pac = ncfile2.variables['tcdq_Pac'][:]
tcdq_ind = ncfile2.variables['tcdq_Ind'][:]
ncfile2.close()


# Need to covert nan to 0 so that the arrays add together properly
atl_nan = pl.isnan(tcdq_atl)
tcdq_atl[atl_nan] = 0.

med_nan = pl.isnan(tcdq_med)
tcdq_med[med_nan] = 0.

# add atlantic & Med+BS arrays so that total flux can be inspected
tcdq_AM = tcdq_atl + tcdq_med


pac_nan = pl.isnan(tcdq_pac)
tcdq_pac[pac_nan] = 0.

ind_nan = pl.isnan(tcdq_ind)
tcdq_ind[ind_nan] = 0.

# Convert lat & lon arrays to radians
lat[:] = pl.radians(lat[:])
lon[:] = pl.radians(lon[:])

#--------------------CREATE LATITUDE HALF-GRID ARRAY HERE----------------------
#---------------------------need to loop over arrays---------------------------
# set up empty array, size one more than lat
lat_half = pl.zeros([lat.shape[0]+1])

# set first & last elements of lat_half seperately as pi/2 & -pi/2
lat_half[0] = (pl.pi)/2
lat_half[-1] = -(pl.pi)/2
# loop over lat_half from index 1 to index -2:
for i in range(1,lat_half.shape[0]-1): # loops over 256
    lat_half[i] = 0.5*(lat[i]+lat[i-1])


#------------------------CALCULATE DELTA PHI ARRAY HERE------------------------
# set up empty array for delta phi, same size as lat
# loop over delta phi:
    # delta_phi = lat_half[i] - lat_half[i+1]


#--------------------------CALCULATE DELTA LAMBDA HERE-------------------------
# delta_lambda = 2*pi/nlon
nlon = lon.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon




#--------------calculate cell areas here, use function from above--------------
# set up empty array for area, size lat_half X lon_half
areas = pl.zeros([lat.shape[0],lon.shape[0]])
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
rho = 10**3 # density of water in kg/m^3
EmP_basin = areas*(tcdq_pac/rho) # divide divergence by density to convert from kg/m^2s to m/s
EmP_total = pl.sum(EmP_basin)
print EmP_total/(10**6) # 1000 times too big :(