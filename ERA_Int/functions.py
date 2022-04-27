# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 10:13:40 2015

@author: np838619

Functions for doing stuff with ERA-Interim netCDF files and calculating net
E-P and E-P-R for ocean basins and latitude bands; creating a half-grid for
a latitude array and calculating grid cell areas & many more.

Last updated: 24/01/16 5:59PM 8th commit
"""

from __future__ import division
import pylab as pl
import os
from matplotlib.colors import Normalize

class MidpointNormalize(Normalize):
    """Found this on the internet. Use to centre colour bar at zero.
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        a, b = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return pl.ma.masked_array(pl.interp(value, a, b))

def PrintFiles(path_to_files,start_of_file):
    """Function to isolate files that which start with a particular string
    within a particular directory then arranges them in the correct order.
    
    Args:
        path_to_files (str): string indicating the directory where files are located
        start_of_file (str): string of the beginning of the filename required
    
    Returns:
        filelist (list): list of required file names in a directory
    """
    list_of_files = os.listdir(path_to_files)
    filelist = []
    for each_file in list_of_files:
        if each_file.startswith(start_of_file):
            filelist.append(each_file)

    return filelist

def NearestIndex(array_in,point_in):
    """Function to the the nearest index to a specified geographic co-ordinate from an
    array of latitude or longitude co-ordinates
    
    Args:
        array_in (array): longitude or latitude array
        point_in (float): longitude or latitude co-ordinate
    
    Returns:
        index (int): index of array which has value closest to point_in
    """
    index = pl.absolute(array_in-point_in).argmin()
    
    return index

def MapCoords(coords,basinmap):
    """Function to convert from the matplotlib figure co-ordinates to the actual
    longitude-latitude co-ordinates in degrees.
    
    Args:
        coords (list): x,y co-ordinates of clicked points from matplotlib figure
        basinmap (Basemap object): map used for clicking
    
    Returns:
        lonlat (array): longitude-latitude co-ordinates of clicked points
    """
    # need empty array, len(coords) X 2:
    boundary = pl.zeros([len(coords),len(coords[0])]) 
    boundary[:] = coords[:] # make every row of boundary a pair of co-ordinates from coords
    # Transform co-ordinates from figure x,y to lon-lat:
    ilon, ilat = basinmap(boundary[:,0],boundary[:,1],inverse=True)
    lonlat=pl.array([ilon,ilat]).T # make array of lon, lat values
    
    return lonlat


def GridPoints(lonlat,lon,lat):
    """Function to convert the co-ordinates of all the points clicked on from
    Basemap longitude-latitudes to their nearest indices in the ERA-Interim
    grid.
    
    Args:
        lonlat (array): Basemap longitude-latitude co-ordinates of clicked points
        lon (array): ERA-Interim longitude array
        lat (array): ERA-Interim latitude array
    
    Returns:
        grid (array): co-ordinates of clicked points in their nearest indices on
                      the ERA-Interim grid, shifted a bit if the value is -ve
    """
    # set up empty array to contain the grid points closest to mouse clicks
    grid = pl.zeros_like(lonlat)
    for i in range(lonlat.shape[0]):
        # convert all the longitude/latitude co-ordinates to their nearest index
        # on the ERA-Interim grid
        grid[i,0] = NearestIndex(lon,lonlat[i,0])
        grid[i,1] = NearestIndex(lat,lonlat[i,1])
    # shift all the longitude points down by 256 (Basemap longitude goes from -180 to 180):
    grid[:,0] = grid[:,0]-len(lat)+1
    for i in range((len(grid[:,0]))):
        # all grid points must be +ve integers, if the longitude points are -ve
        # then add 511 to them:
        if grid[i,0] < 0:
            grid[i,0] = grid[i,0] + len(lon) - 1
    
    return grid


def SplitGrid(grid,lonlat):
    """Function to split indices array of clicked points in half into an array
    representing the western boundary and another array representing the eastern
    boundary.
    
    Args:
        grid (array): co-ordinates of clicked points in terms of their nearest
                      longitude-latitude indices on the ERA-Interim grid
    
    Returns:
        grid_west (array): co-ordinates of clicked points in terms of their
                           nearest longitude-latitude indices on the ERA-Interim
                           grid for the western boundary of an ocean basin
        grid_east (array): co-ordinates of clicked points in terms of their
                           nearest longitude-latitude indices on the ERA-Interim
                           grid for the eastern boundary of an ocean basin
    """
    # top half of grid array represents western boundary,
    # bottom half represents eastern boundary
    # if western & eastern boundaries are both north down then split grid array 
    # in half without flipping
    if grid[0,1]!=grid[grid.shape[0]-1,1]:
        grid_west = grid[:grid.shape[0]/2]
        grid_east = grid[lonlat.shape[0]/2:]
    # if western boundary is north down & eastern boundary is south up then 
    # flip the eastern boundary
    else:
        grid_west = grid[:grid.shape[0]/2]
        grid_east = pl.flipud(grid[grid.shape[0]/2:])
    
    return grid_west, grid_east
   
   
def MaskOutsideGrid(grid_west,grid_east,lon,lat):
    """Function to create an array the same size as the ERA-Interim grid with
    grid points = 1 if inside mask & = 0 if outside mask.
    
    Args:
        grid_west (array): co-ordinates of clicked points in terms of their
                           nearest longitude-latitude indices on the ERA-Interim
                           grid for the western boundary of an ocean basin
       grid_east (array): co-ordinates of clicked points in terms of their
                           nearest longitude-latitude indices on the ERA-Interim
                           grid for the eastern boundary of an ocean basin
       lon (array): ERA-Interim longitude array
       lat (array): ERA-Interim latitude array
   
    Returns:
        eragrid (array): array same size as ERA-Interim grid with grid points = 1
                         inside the mask & = 0 outside the mask
    """
    # need an empty grid the same size as era-interim grid, lons x lats (as this is what matplotlib does)
    eragrid = pl.zeros([len(lon),len(lat)])
    for i in range(len(grid_west[:,0])): # loop along latitudes required for mask
    # every point in eragrid array between longitude & latitudes indices specified in grid_west & grid_east = 1
    # if eastern boundary longitude point is east of 0 meridian then eragrid must be set to 1 seperately
    # set eragrid=1 from western longitude point to 0 meridian, then from 0 meridian to eastern longitude point
        if grid_east[i,0] < grid_west[i,0]:
            eragrid[grid_west[i,0]:,grid_east[i,1]]=1.0
            eragrid[:grid_east[i,0]+1,grid_east[i,1]]=1.0
    # if the eastern boundary point is west of the 0 meridian then set eragrid=1
    # from the western to eastern longitude point
        else:
            eragrid[grid_west[i,0]:grid_east[i,0],grid_east[i,1]]=1.0
    
    return eragrid


def Zero2Nan(array_in):
    """Function to create an array with all points inside the ocean basin = 1 
    & everywhere else = nan. This is useful for when saving in netCDF as nan
    doesn't interfere with the colourbar (shows up white in ncview). Should be
    done before plotting as shiftgrid messes things up.
    
    Args:
        array_in (array): input array where all points inside ocean basin = 0 & 
        points inside = 1.
    
    Returns:
        array (array): output array where all points inside the ocean basin = 1 
    & everywhere else = nan. Same dimensions as input.
    """
    newarray = pl.zeros_like(array_in)
    for i in range(newarray.shape[0]):
        for j in range(newarray.shape[1]):
            if array_in[i,j] == 0.:
                newarray[i,j] = pl.float32('nan')
            else:
                newarray[i,j] = array_in[i,j]
    
    return newarray


def BasinMask(newarray):
    """Function to make the final mask for an ocean basin that will be saved in
    a netCDF file. Everything outside the mask = 0 & everything inside = 1.
    
    Args:
        newarray (array): array for whole planet, same size as ERA-Interim grid,
                          inside basin = 1 & outside basin = NaN

    Returns:
        basinmask (array): final mask for an ocean basin, points inside = 1 & 
                           points outside  = 0
    """
    basinmask = pl.ones_like(newarray) # array of ones
    nanarray = pl.isnan(newarray) # find where the NaNs are in newarray
    basinmask[nanarray] = 0. # set the same locations to zero in basin mask
    
    return basinmask
    

def HalfGrid(lat):
    """Function to create latitude half-grid
    
    Args:
        lat (array): latitude array in radians
    
    Returns:
        lath (array): half-grid latitude array in radians
    """
    # set up empty array, size one more than lat
    lath = pl.zeros([lat.shape[0]+1])
    # set first & last elements of lat_half seperately as pi/2 & -pi/2:
    lath[0] = (pl.pi)/2; lath[-1] = -(pl.pi)/2
    # loop over lat_half from index 1 to index -2:
    for i in range(1,lath.shape[0]-1): # loops over 256
        lath[i] = 0.5*(lat[i]+lat[i-1])
    
    return lath


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


def NaN2Zero(tcdq):
    """Function to convert all the nans in tcdq array into zeros
    
    Args:
        tcdq (array): moisture flux divergence array in kg/m^2s with masked
                    locations assigned as nan
    
    Returns:
        tcdq (array): moisture flux divergence array in kg/m^2 with masked
                    locations assigned as zero
    """
    # find where the NaNs are in the array:
    nanarray = pl.isnan(tcdq)
    # change the NaNs to zero:
    tcdq[nanarray] = 0.
    
    return tcdq


def NetEmPCalc(area_array,tcdq,rho):
    """Function to calculate the net evaporation minus precipitation for an 
    ocean basin.
    
    Args:
        area_array (array): array of grid cell areas in m^2
        tcdq (array): array of moisture flux divergence in kg/m^2s
        rho (float): predefined density of water, usually 10^3 kg/m^3
    
    Returns:
        total (float): net surface water flux of ocean in Sverdrups 
                        (Sv, 10^6 m^3/s)
    """
    
    # multiply the area by tcdq/rho to get E-P values in m^3/s
    # dividing tcdq by density coverts the units from kg/m^2 into m/s
    basin = area_array*(tcdq/rho)
    # add up all the E-P across the basin and divide by 10^6 to convert units
    # to Sverdrups (1Sv = 10^6 m^3/s)
    total = pl.nansum(basin)/(10**6)
    
    return total


def NetEmPBands(indices,areas,tcdq,rho):
    """Function to calculate net E-P for latitude bands
    
    Args:
        indices (array): indices from ERA-Interim latitude grid which have the
                        value closest to the physical latitudes used in literature
                        to split up the oceans
        areas (array): areas of grid cells in m^2
        tcdq (array): moisture flux divergence in kg/m^2s for particular ocean
        
    Returns:
        bands (array): E-P values for latitude bands in Sverdrups (Sv, 10^6 m^3/s)
    """
    # set up empty array for E-P values, one element smaller than indices:
    bands = pl.zeros([indices.shape[0]-1])
    for i in range(indices.shape[0]-1):
        indpair = (indices[i],indices[i+1]) # select the restricting latitudes
        areas_band = areas[indpair[0]:indpair[1]] # restrict the areas array
        tcdq_band = tcdq[indpair[0]:indpair[1]] # restrict the divergence array
        bands[i] = NetEmPCalc(areas_band,tcdq_band,rho) # calculate net E-P
        
    return bands


def EmPyears(tcdq_years,basin_mask,eramask,areas,rho):
    """Function to calculate net E-P for each year over an ocean basin.
    
    Args:
        tcdq_years (array): yearly means of moisture flux divergence
        basin_mask (array): ones indicate grid points required for a particular
                            basin
        eramask (array): land-sea mask from ERA-Interim, values range from 0-1
                        where 0 is a grid point with no land coverage and 1 is
                        is a grid point with complete land coverage
        areas (array): area of each grid cell
        rho (int): density (kg/m^3) of water, divide by rho to change from
                    kg/m^3s to m/s
    
    Returns:
        EmP_years (array): net E-P (Sverdrups, Sv) for each year over an ocean
                            basin
    """
    basin_years = tcdq_years*basin_mask*(1-eramask)
    EmP_years = pl.zeros([tcdq_years.shape[0]])
    for i in range(tcdq_years.shape[0]):
        EmP_years[i] = NetEmPCalc(areas,basin_years[i],rho)
    
    return EmP_years
 
   
def TotalRunoff(runlats,eralats,indices,basin_runoff):
    """Function to calculate the total runoff for a given ocean basin between
    two latitude points
    
    Args:
        runlats (array): latitudes from runoff csv file
        eralats (array): latitudes from ERA-Interim grid
        indices (tuple): indices of the latitudes we want to restrict the runoff
                        between
        basin_runoff (array): runoff in 1deg latitude bands for specific ocean basin
    
    Returns:
        total (float): total runoff between latitudes specificed
    """
    
    # find the nearest indices in runlats to the ERA-Interim latitude at each of
    # the indices closest to the physical limits from literature
    run_ind0 = NearestIndex(runlats,eralats[indices[0]])
    run_ind1 = NearestIndex(runlats,eralats[indices[1]])
    # add up all the runoff values between the 2 runlat indices
    total = pl.sum(basin_runoff[run_ind0:run_ind1+1])
    
    return total


def EPRbands(EmP_bands,indices,runlats,eralats,basin_runoff):
    """Function to calculate net E-P-R for latitude bands
    
    Args:
        EmP_bands (array): Net E-P values (Sv) for the relevant latitude bands
        indices (array): indices for ERA-Interim latitude grid which have the
                        value closest to the physical latitudes used in literature
                        to split up the oceans
        runlats (array): latitudes from runoff csv file
        eralats (array): latitudes from ERA-Interim grid
        basin_runoff (array): runoff in 1deg latitude bands for specific ocean basin
    
    Returns:
        EPR_bands (array): net E-P-R values for latitude bands in Sverdrups
    """
    # set up empty array for E-P-R, same size as E-P array
    EPR_bands = pl.zeros_like(EmP_bands)
    for i in range(indices.shape[0]-1):
        indpair = (indices[i],indices[i+1]) # select the restricting latitudes
        # calculate E-P-R for each band, use TotalRunoff function
        EPR_bands[i] = EmP_bands[i] - TotalRunoff(runlats,eralats,indpair,basin_runoff)
    
    return EPR_bands


def dtcwv_dt(tcwv_months):
    """Function to calculate the d(TCWV)/dt term in the moisture budget:
                d(TCWV)/dt + TCDQ = E-P,
    where TCWV is total column water vapour & TCDQ is vertically integrated
    moisture flux divergence.
    
    Args:
        tcwv_months (array): monthly mean total column water vapour from the
                             ERA-Interim reanalysis (kg/m^2)
    
    Returns:
        tcwv_months_dt (array): the time derivative term in the moisture budget
                                for each month (kg/m^2s)
    """
    time_step = 86400*(365/12) # no. of seconds in a day X no. of days in a month
    # empty array for time derivatives, same size as tcwv_months:
    tcwv_months_dt = pl.zeros_like(tcwv_months)
    # loop from -1 (Dec) to shape-1 (Nov) to include all months in loop
    for month in range(-1,tcwv_months.shape[0]-1):
        # centred difference i.e. (Mar-Jan)/2*time_step
        tcwv_months_dt[month] = (tcwv_months[month+1] - tcwv_months[month-1])/(2*time_step)
    
    return tcwv_months_dt
    

def SeasonalMeans(variable_array):
    """Function to calculate the climatological seasonal means (DJF,MAM,JJA,SON)
    of a variable.
    
    Args:
        variable_array (array): variable of which the climatological seasonal 
                                means are required
    
    Returns:
        variable_seasons (array): climatological seasonal means of the input 
                                  variable
    """
    # Calculate the mean of each trio of months for each year, missing out the 
    # first year as there is no data for December the year before the data starts
    MAM = pl.mean(variable_array[:,2:5],axis=1) # March, April, May
    JJA = pl.mean(variable_array[:,5:8],axis=1) # June, July, August
    SON = pl.mean(variable_array[:,8:11],axis=1) # September, October, November
    DJF = pl.zeros_like(MAM) # December, January, February
    DJF[0] = (variable_array[-1,-1]+pl.sum(variable_array[0,:2],axis=0))/3
    for i in range(1,DJF.shape[0]):
        DJF[i] = (variable_array[i-1,-1]+pl.sum(variable_array[i,:2],axis=0))/3
    
    # Calculate the climatological mean of each season:  
    MAM_mn = pl.mean(MAM,axis=0)
    JJA_mn = pl.mean(JJA,axis=0)
    SON_mn = pl.mean(SON,axis=0)
    DJF_mn = pl.mean(DJF,axis=0)
    
    # Stick all seasons in one array before outputting:
    variable_seasons = pl.array([DJF_mn,MAM_mn,JJA_mn,SON_mn])
    
    return variable_seasons
    
def FluxScaled(areas,basin_mask,flux):
    """Function to scale the basin-integrated freshwater flux of an ocean by 
    its area.
    
    Args:
        areas (array): area of each grid cell
        basin_mask (array): mask for a particular ocean basin, 1 represents a 
                            point inside the basin & 0 a point outside it
        flux (float): value of net freshwater flux (Sv) for an ocean basin
    
    Returns:
        flux_scaled (float): value of net freshwater flux for an ocean basin
                             scaled by the area of the basin (Sv/m^2)
    """
    basin_areas = areas*basin_mask # total area of ocean basin
    flux_scaled = flux/pl.sum(basin_areas) # net freshwater flux scaled by area
    
    return flux_scaled

def FluxBandScaled(areas,basin_mask,indices,flux_bands):
    """Function to scale the net freshwater flux of a latitude band in an ocean
    basin by its area.
    
    Args:
        areas (array): area of each grid cell
        basin_mask (array): mask for a particular ocean basin, 1 represents a 
                            point inside the basin & 0 a point outside it
        indices (array): indices for ERA-Interim latitude grid which have the
                         value closest to the physical latitudes used in literature
                         to split up the oceans
        flux_bands (array): net freshwater flux (Sv) for each latitude band in
                            an ocean basin
    
    Returns:
        flux_scaled (array): net freshwater flux (Sv/m^2) for each latitude band
                             scaled by the band's area
    """
    basin_areas = areas*basin_mask # total area of ocean basin
    # need empty array for scaled flux, same size as flux_bands
    flux_scaled = pl.zeros_like(flux_bands)
    for i in range(flux_bands.shape[0]):
        indpair = (indices[i],indices[i+1]) # bounding latitude indices for a band
        band_area = pl.sum(basin_areas[indpair[0]:indpair[1]+1]) # area of latitude band
        flux_scaled[i] = flux_bands[i]/band_area # net freshwater flux scaled by area
        
    return flux_scaled