# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 10:16:39 2015

@author: np838619

Code to calculate the yearly mean E-P from ERA-Interim reanalysis using monthly
mean TCDQ (vertically integrated moisture flux divergence). The time derivative
term of the moisture budget is not included as it is equivalent to zero on the
time scales considered here. Time series of the inter-annual variability for
each ocean basin and latitude bands therein are produced and standard deviations
representing inter-annual variabilty for each basin and latitude band are also
presented.

Last updated: 24/01/16 6:37PM 5th commit
"""

# import stuff up here: pylab, netcdf4, functions
from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from functions import * 


# open up netcdf files here:
    # tcdq_36years_means.nc has the mean tcdq for each year
    # tcdq_basins2.nc has tcdq for Atlantic (35S-45N) & Med/BS
    # tcdq_basins3.nc has tcdq for Pacific (30S-47N) and Indian (>35S)

yearfile = Dataset('tcdq_36years_means.nc','r') # file containing mean tcdq for each year
tcdq_years = yearfile.variables['tcdq'][:]
lat = yearfile.variables['lat'][:] # need lat & lon for areas and restricting to basins
lon = yearfile.variables['lon'][:]
yearfile.close()

basinfile = Dataset('tcdq_basins_60.nc','r') # file containing masks for Atlantic & Med/BS
atlmask = basinfile.variables['maskA'][:] # Atlantic Ocean mask
balmask = basinfile.variables['maskB'][:] # Baltic Sea mask
medmask = basinfile.variables['maskM'][:] # Med + Black Seas mask
pacmask = basinfile.variables['maskP'][:] # Pacific Ocean mask
indmask = basinfile.variables['maskI'][:] # Indian Ocean mask
basinfile.close()

maskfile = Dataset('tcdq_ann_mean_ERAI.nc','r')
eramask = maskfile.variables['LSM'][:]
maskfile.close()


# add the Atlantic and Med/BS masks together:
Amask = atlmask + balmask + medmask


# calculate halfgrid here:
lat_rad = pl.radians(lat[:])
lon_rad = pl.radians(lon[:])
lat_half = HalfGrid(lat_rad)

nlon = lon_rad.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon

#----------------------calculate the areas of each grid cell-------------------
# set up empty array for area, size lat_half X lon_half
areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
radius = 6.37*(10**6) # radius of Earth, in metres
for i in range(lat_half.shape[0]-1):
    latpair = (lat_half[i],lat_half[i+1])
    for j in range(lon.shape[0]):
        areas[i,j] = AreaFullGaussianGrid(radius,delta_lambda,latpair)


#-------------------calculate net E-P for Atlantic for each year---------------
rho = 10**3
atlEmP_years = EmPyears(tcdq_years,Amask,eramask,areas,rho)
pacEmP_years = EmPyears(tcdq_years,pacmask,eramask,areas,rho)
indEmP_years = EmPyears(tcdq_years,indmask,eramask,areas,rho)



#----------------------CALCULATE STANDARD DEVIATIONS---------------------------
atl_sd = pl.std(atlEmP_years)
pac_sd = pl.std(pacEmP_years)
ind_sd = pl.std(indEmP_years)

#-------------------------DO THE SAME FOR LATITUDE BANDS-----------------------
 # Arrays for latitude boundaries:

atl_lats = pl.array([60,45,24,-16,-35])
pac_lats = pl.array([65.5,47,24,-17,-30])
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
    
atl_tcdq_yrs = Amask*tcdq_years
pac_tcdq_yrs = pacmask*tcdq_years
ind_tcdq_yrs = indmask*tcdq_years
    
# Empty arrays for lat bands for entire ERA-I period:
atlbands_yrs = pl.zeros([tcdq_years.shape[0],atl_ind.shape[0]-1])
pacbands_yrs = pl.zeros([tcdq_years.shape[0],pac_ind.shape[0]-1])
indbands_yrs = pl.zeros([tcdq_years.shape[0],ind_ind.shape[0]-1])
# Loop over 36 years:
for i in range(tcdq_years.shape[0]):
    atlbands_yrs[i] = NetEmPBands(atl_ind,areas,atl_tcdq_yrs[i],rho)
    pacbands_yrs[i] = NetEmPBands(pac_ind,areas,pac_tcdq_yrs[i],rho)
    indbands_yrs[i] = NetEmPBands(ind_ind,areas,ind_tcdq_yrs[i],rho)


atlbands_sd = pl.std(atlbands_yrs,axis=0)
pacbands_sd = pl.std(pacbands_yrs,axis=0)
indbands_sd = pl.std(indbands_yrs,axis=0)


#------------------PLOT TIME SERIES FOR EACH BASIN------------------------------
# plot the time series of yearly means for all basins & each basin's latitude
# bands in the same subplot
clrs = ['b','g','r','k']
fill = ['#B2E6FF','#D1FFA3','#FFC2C2','gainsboro']
ax, fig = pl.subplots(2,2,figsize=(12,10))

years = pl.linspace(0,35,36)
#-all basins-------------------------------------------------------------------
ax1 = pl.subplot(2,2,1)
pl.plot(years,atlEmP_years,label="Atlantic")
pl.fill_between(years,atlEmP_years-atl_sd,atlEmP_years+atl_sd,color='#B2E6FF')
pl.plot(years,pacEmP_years,label="Pacific")
pl.fill_between(years,pacEmP_years-pac_sd,pacEmP_years+pac_sd,color='#D1FFA3')
pl.plot(years,indEmP_years,label="Indian")
pl.fill_between(years,indEmP_years-ind_sd,indEmP_years+ind_sd,color='#FFC2C2')
pl.axhline(linewidth=1,ls='--',color='k')
ax1.set_ylabel('$E-P$ (Sv)',fontsize=16)
pl.yticks(pl.linspace(-0.2,1.4,5))
ax1.set_xticklabels(pl.linspace(1979,2014,8).astype(int)) # set the xtick labels to  be years
ax1.xaxis.set_major_formatter(pl.NullFormatter()) # don't want labels on this subplot
pl.title('Inter-annual $E-P$  variability',fontsize=12)
ax1.legend(loc=2,ncol=3,fontsize=12,labelspacing=0.5,columnspacing=0.47,
           mode="expand")

#-Atlantic---------------------------------------------------------------------
atl_labels = ['47$^\circ$N-60$^\circ$N','24$^\circ$N-45$^\circ$N',
              '16$^\circ$S-24$^\circ$N','35$^\circ$S-16$^\circ$S']
ax2 = pl.subplot(2,2,2)
for i in range(atl_ind.shape[0]-1):
    pl.plot(years,atlbands_yrs[:,i],label=atl_labels[i],color=clrs[i])
    pl.fill_between(years,atlbands_yrs[:,i]-atlbands_sd[i],atlbands_yrs[:,i]+atlbands_sd[i],
                    color=fill[i])
    #pl.fill_between(years,atlbands_yrs[:,i]-atlbands_sd[i],
    #                atlbands_yrs[:,i]+atlbands_sd[i],color='#B2E6FF')
pl.axhline(linewidth=1,ls='--',color='k')
ax2.set_xticklabels(pl.linspace(1979,2014,8).astype(int))
ax2.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylim(-0.25,0.75) # change y-axis to fit legend
#pl.yticks(pl.linspace(-0.15,0.65,5))
pl.matplotlib.ticker.FormatStrFormatter('%.2f')
pl.title('Atlantic latitude bands',fontsize=12)
ax2.legend(loc=(0.015,0.275),ncol=2,fontsize=12,columnspacing=0.5,labelspacing=0.4)

#-Pacific----------------------------------------------------------------------
pac_labels = ['47$^\circ$- BS','24$^\circ$N-47$^\circ$N',
              '17$^\circ$S-24$^\circ$N','30$^\circ$S-17$^\circ$S']
ax3 = pl.subplot(2,2,3)
for i in range(pac_ind.shape[0]-1):
    pl.plot(years,pacbands_yrs[:,i],label=pac_labels[i],color=clrs[i])
    pl.fill_between(years,pacbands_yrs[:,i]-pacbands_sd[i],pacbands_yrs[:,i]+pacbands_sd[i],
                    color=fill[i])
pl.axhline(linewidth=1,ls='--',color='k')
# reduce fontsize of a-axis tick labels and rotate them 45deg so they fit
ax3.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=10,rotation=45)
ax3.set_xlabel('Years',fontsize=12)
pl.ylim(-0.45,0.75)
#pl.yticks(pl.linspace(-0.2,0.6,5))
ax3.set_ylabel('$E-P$ (Sv)',fontsize=16)
pl.title('Pacific latitude bands',fontsize=12)
ax3.legend(loc=2,ncol=3,fontsize=12,columnspacing=0.375,handletextpad=0.4,
           mode="expand")

#-Indian-----------------------------------------------------------------------
ind_labels = ['>8$^\circ$S','20$^\circ$S-8$^\circ$S','35$^\circ$S-20$^\circ$S']
ax4 = pl.subplot(2,2,4)
for i in range(ind_ind.shape[0]-1):
    pl.plot(years,indbands_yrs[:,i],label=ind_labels[i],color=clrs[i])
    pl.fill_between(years,indbands_yrs[:,i]-indbands_sd[i],indbands_yrs[:,i]+indbands_sd[i],
                    color=fill[i])
pl.axhline(linewidth=1,ls='--',color='k')
ax4.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=10,rotation=45)
ax4.set_xlabel('Years',fontsize=12)
#pl.yticks(pl.linspace(-0.3,0.9,5))
pl.ylim(-0.25,0.75)
pl.title('Indian latitude bands',fontsize=12)
ax4.legend(loc=2,ncol=3,fontsize=12,labelspacing=0.5,columnspacing=0.5)

pl.subplots_adjust(wspace=0.2,hspace=0.2,top=0.94,bottom=0.1,left=0.11,right=0.93)
#------------------------------------------------------------------------------


#---------------------WRITE STANDARD DEVIATIONS TO FILE------------------------

#f = open('annmean_EmP_stdev60.txt','w')
#f.write('Standard deviations of annual mean E-P for the Atlantic (35S-45N), ' + 
#'Pacific (30S-47N) and Indian (>35S) Oceans. Calculated from the E-P for each year' + 
#'in the ERA-Interim reanalysis (1979-2014). These standard deviations represent ' + 
#'the interannual variability of E-P for the reanalysis period. Also presented are ' +
#'standard deviations for various latitude bands within each ocean. \n\n' + 
#'Atlantic (35S-60N) Ocean: ')
#atl_sd.tofile(f,sep=" ",format="%.3f")
#f.write('\n\nPacific (30S-BS) Ocean: ')
#pac_sd.tofile(f,sep=" ",format="%.3f")
#f.write('\n\nIndian (>35S) Ocean: ')
#ind_sd.tofile(f,sep=" ",format="%.3f")
#f.write('\n\nAtlantic bands: \n45N-60N     24N-45N     16S-24N     35S-16S\n')
#atlbands_sd.tofile(f,sep="   ",format="%.3f")
#f.write('\n\nPacific bands: \n47N-BS    24N-47N   17S-24N   30S-17S \n')
#pacbands_sd.tofile(f,sep="   ",format="%.3f")
#f.write('\n\nIndian bands: \n>8S    20S-8S  35S-20S\n')
#indbands_sd.tofile(f,sep="   ",format="%.3f")
#f.close()