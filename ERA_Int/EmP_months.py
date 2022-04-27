# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 11:19:56 2015

@author: np838619

Code to calculate annual mean monthly E-P and E-P-R from ERA-Interim reanalysis
& Dai & Trenberth (2002) continental discharge dataset. E-P is calculated from
the moisture budget using monthly means of TCWV (total column water vapour) &
TCDQ (vertically integrated moisture flux divergence). Time series for each basin
and latitude bands therein are produced along with standard deviations (representing
intra-annual variability) for each basin and latitude band.

Last updated: 2/09/15 16:25PM 6th commit
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import pandas
from functions import *


def dtcwv_dt(tcwv_months):
    """
    """
    time_step = 86400*30 # no. of seconds in a day X no. of days in a month
    tcwv_months_dt = pl.zeros_like(tcwv_months)
    for month in range(-1,tcwv_months.shape[0]-1):
        tcwv_months_dt[month] = (tcwv_months[month+1] - tcwv_months[month-1])/(2*time_step)
    
    return tcwv_months_dt



# read in nc files with basin masks: tcdq_basins2 for Atl & Med, tcdq_basins3 for Pac & Ind
# extract lat & lon
# also need ERA-Interim mask
basinfile = Dataset('tcdq_basins_60.nc','r')
atlmask = basinfile.variables['maskA'][:] # mask for Atlantic Ocean
balmask = basinfile.variables['maskB'][:] # mask for Baltic Sea
pacmask = basinfile.variables['maskP'][:] # mask for Pacific Ocean
indmask = basinfile.variables['maskI'][:] # mask for Indian Ocean
lon = basinfile.variables['lon'][:] # ERA-Interim longitude
lat = basinfile.variables['lat'][:] # ERA-Interim latitude
basinfile.close()

#IPbasinfile = Dataset('tcdq_basins3.nc','r')
#pacmask = IPbasinfile.variables['maskP'][:] # mask for Pacific Ocean
#indmask = IPbasinfile.variables['maskI'][:] # mask for Indian Ocean
#IPbasinfile.close()

maskfile = Dataset('/panfs/jasmin/era/era-in/netc/ggis/1989/jan1989/ggis198901010000.nc','r')
eramask = maskfile.variables['LSM'][:] # ERA-Interim land-sea mask
maskfile.close()


# read in csv files for monthly runoff: Atlantic, Pacific, Indian & Med:

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
           

# STICK ALL THE MONTHLY RUNOFF ARRAYS INTO ONE ARRAY PER BASIN:
Arun_mnths = pl.array([atl_run.Jan,atl_run.Feb,atl_run.Mar,atl_run.Apr,
                       atl_run.May,atl_run.Jun,atl_run.Jul,atl_run.Aug,
                       atl_run.Sep,atl_run.Oct,atl_run.Nov,atl_run.Dec])

#Mrun_mnths = pl.array([med_run.Jan,med_run.Feb,med_run.Mar,med_run.Apr,
#                       med_run.May,med_run.Jun,med_run.Jul,med_run.Aug,
#                       med_run.Sep,med_run.Oct,med_run.Nov,med_run.Dec])

Prun_mnths = pl.array([pac_run.Jan,pac_run.Feb,pac_run.Mar,pac_run.Apr,
                       pac_run.May,pac_run.Jun,pac_run.Jul,pac_run.Aug,
                       pac_run.Sep,pac_run.Oct,pac_run.Nov,pac_run.Dec])

Irun_mnths = pl.array([ind_run.Jan,ind_run.Feb,ind_run.Mar,ind_run.Apr,
                       ind_run.May,ind_run.Jun,ind_run.Jul,ind_run.Aug,
                       ind_run.Sep,ind_run.Oct,ind_run.Nov,ind_run.Dec])

# WANT TO ADD THE MED RUNOFF ONTO ATLANTIC:
    # GIBRALTAR AT ~36N SO FIND NEAREST INDEX IN RUNOFF LATITUDE ARRAY
    # THEN ADD THE MED RUNOFF TO THE ATLANTIC RUNOFF AT THAT INDEX
#Gind = NearestIndex(atl_run.lat,36)
#Arun_mnths[:,Gind] = Arun_mnths[:,Gind] #+ Mrun_mnths[:,0]

# READ IN ALL THE GGAW NETCDF FILES HERE
# USE PrintFiles FUNCTION 
# NEED AN ARRAY OF THE YEARS
years = pl.linspace(1979,2014,36).astype(int)
year_input = [str(i) for i in years] # convert to strings fro use in path
# NEED TO SET UP EMPTY ARRAY FOR FILE NAMES
filenames = pl.zeros([len(years),12],dtype='S17')
# LOOP OVER YEARS TO GET ALL FILE NAMES USING PrintFiles
for year in range(len(year_input)):
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + year_input[year]
    filenames[year] = PrintFiles(path,'ggaw')
# NEED ARRAY FOR TCDQ
tcdq = pl.zeros([len(year_input),len(filenames[0]),1,1,256,512])
tcwv = pl.zeros_like(tcdq)
# LOOP OVER YEARS AND EXTRACT TCDQ FROM NC FILES USING FILENAMES ARRAY
for year in range(len(year_input)):
    for name in range(len(filenames[year])):
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                    str(year_input[year]) + '/' + str(filenames[year,name]),'r')
        tcdq[year,name] = ncfile.variables['TCDQ'][:]
        tcwv[year,name] = ncfile.variables['TCWV'][:]
        ncfile.close()


#TRY TO REMOVE THE LEVEL AXIS (& THERE MIGHT BE ANOTHER USELESS AXIS)
tcdq = pl.squeeze(tcdq,axis=(2,3))
tcwv = pl.squeeze(tcwv,axis=(2,3))
eramask = pl.squeeze(eramask,axis=(0,1))

# TAKE MEAN OF TCDQ ALONG 0 AXIS TO GET MEAN OF EACH MONTH FOR 36 YEARS
tcdq_mnths = pl.mean(tcdq,axis=0)
tcwv_mnths = pl.mean(tcwv,axis=0)

# for monthly means use the moisture budget:
tcwv_mnths_dt = dtcwv_dt(tcwv_mnths)
EmP_mm = tcwv_mnths_dt + tcdq_mnths


# which latitudes will be required for splitting the basins into subregions?
atl_lats = pl.array([60,45,24,-16,-35])
pac_lats = pl.array([65.5,47,24,-30])
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


# add the Atlantic and Baltic masks together:
Amask = atlmask + balmask

# calculate net E-P for each month using EmPyears (works the same for months)
rho = 10**3 # density of water
atlEmP_mnths = EmPyears(EmP_mm,Amask,eramask,areas,rho)
pacEmP_mnths = EmPyears(EmP_mm,pacmask,eramask,areas,rho)
indEmP_mnths = EmPyears(EmP_mm,indmask,eramask,areas,rho)

# INTRA-ANNUAL E-P VARIABILITY FOR EACH BASIN:
AEmP_mm_sd = pl.std(atlEmP_mnths)
PEmP_mm_sd = pl.std(pacEmP_mnths)
IEmP_mm_sd = pl.std(indEmP_mnths)


#---------------------CALCULATE NET E-P-R FOR EACH MONTH-----------------------

# use the indices of the northern & southern boundaries for the relevant region of each ocean basin 
atl_lim = (atl_ind[0],atl_ind[-1])
pac_lim = (pac_ind[0],pac_ind[-1])
ind_lim = (ind_ind[0],ind_ind[-1])
atl_EPR_mnths = pl.zeros_like(atlEmP_mnths)
pac_EPR_mnths = pl.zeros_like(pacEmP_mnths)
ind_EPR_mnths = pl.zeros_like(indEmP_mnths)
for i in range(atlEmP_mnths.shape[0]): # loop over 12 months
    atl_EPR_mnths[i] = atlEmP_mnths[i] - TotalRunoff(atl_run.lat,lat,atl_lim,Arun_mnths[i])
    pac_EPR_mnths[i] = pacEmP_mnths[i] - TotalRunoff(pac_run.lat,lat,pac_lim,Prun_mnths[i])
    ind_EPR_mnths[i] = indEmP_mnths[i] - TotalRunoff(ind_run.lat,lat,ind_lim,Irun_mnths[i])


# CALCULATE THE INTRA-ANNUAL E-P-R VARIABILITY FOR EACH BASIN:
A_EPR_sd = pl.std(atl_EPR_mnths)
P_EPR_sd = pl.std(pac_EPR_mnths)
I_EPR_sd = pl.std(ind_EPR_mnths)


#------------------TIME SERIES OF THE INTRA-ANNUAL VARIABILITY-----------------
# list of months for plotting
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

#put E-P & E-P-R on same plot
#------------------------------------------------------------------------------
mx = pl.subplots(2,1)
# E-P
ax1 = pl.subplot(2,1,1)
pl.plot(pl.linspace(0,11,12),atlEmP_mnths,label="Atlantic")
pl.fill_between(pl.linspace(0,11,12),atlEmP_mnths-AEmP_mm_sd,
                atlEmP_mnths+AEmP_mm_sd,color='#B2E6FF')
pl.plot(pl.linspace(0,11,12),pacEmP_mnths,label="Pacific")
pl.fill_between(pl.linspace(0,11,12),pacEmP_mnths-PEmP_mm_sd,
                pacEmP_mnths+PEmP_mm_sd,color='#D1FFA3',alpha=0.8)
pl.plot(pl.linspace(0,11,12),indEmP_mnths,label="Indian")
pl.fill_between(pl.linspace(0,11,12),indEmP_mnths-IEmP_mm_sd,
                indEmP_mnths+IEmP_mm_sd,color='#FFC2C2',alpha=0.65)
pl.axhline(linewidth=1, color='black',ls='--') # black dashed line at E=P+R
ax1.set_ylabel('$E-P$ (Sv)',fontsize=16)
pl.xticks(pl.linspace(0,11,12)) # set to 12 ticks on x-axis
pl.ylim(-0.2,1.6)
pl.yticks(pl.linspace(-0.2,1.6,4))
ax1.xaxis.set_major_formatter(pl.NullFormatter()) # get rid of labels for upper plot x-axis
pl.title('Intra-annual $E-P$ variability')
# only have legend in upper plot, mess around to make it look nice
ax1.legend(loc=(0.275,0.8),ncol=3,columnspacing=0.5,handletextpad=0.6,fontsize=13)

# E-P-R
ax2 = pl.subplot(2,1,2)
pl.plot(pl.linspace(0,11,12),atl_EPR_mnths,label="Atlantic")
pl.fill_between(pl.linspace(0,11,12),atl_EPR_mnths-A_EPR_sd,
                atl_EPR_mnths+A_EPR_sd,color='#B2E6FF')
pl.plot(pl.linspace(0,11,12),pac_EPR_mnths,label="Pacific")
pl.fill_between(pl.linspace(0,11,12),pac_EPR_mnths-P_EPR_sd,
                pac_EPR_mnths+P_EPR_sd,color='#D1FFA3',alpha=0.8)
pl.plot(pl.linspace(0,11,12),ind_EPR_mnths,label="Indian")
pl.fill_between(pl.linspace(0,11,12),ind_EPR_mnths-I_EPR_sd,
                ind_EPR_mnths+I_EPR_sd,color='#FFC2C2',alpha=0.65)
pl.axhline(linewidth=1, color='black',ls='--') # black dashed line at E=P+R
ax2.set_ylabel('$E-P-R$ (Sv)',fontsize=16)
pl.xticks(pl.linspace(0,11,12))
pl.ylim(-0.7,1.1)
pl.yticks(pl.linspace(-0.7,1.1,4))
ax2.set_xticklabels(months)
ax2.set_xlabel('Months',fontsize=16) # only label bottom plot x-axis
pl.title('Intra-annual $E-P-R$ variability')
pl.subplots_adjust(top=0.92)
#------------------------------------------------------------------------------


#-----------------------DO THE SAME FOR LATITUDE BANDS-------------------------
# set up empty arrays, size 12X3 (2 for Pacific), for each basin:
atl_bnds_mnths = pl.zeros([EmP_mm.shape[0],atl_ind.shape[0]-1])
pac_bnds_mnths = pl.zeros([EmP_mm.shape[0],pac_ind.shape[0]-1])
ind_bnds_mnths = pl.zeros([EmP_mm.shape[0],ind_ind.shape[0]-1])

# Get monthly mean E-P for each basin by multiplying by the basin masks &
# 1 - ERA-Interim mask as I can't be bothered writing a new function
atl_EmP_mm = Amask*EmP_mm*(1-eramask)
pac_EmP_mm = pacmask*EmP_mm*(1-eramask)
ind_EmP_mm = indmask*EmP_mm*(1-eramask)

for i in range(EmP_mm.shape[0]): # loop over 12 months
    atl_bnds_mnths[i] = NetEmPBands(atl_ind,areas,atl_EmP_mm[i],rho)
    pac_bnds_mnths[i] = NetEmPBands(pac_ind,areas,pac_EmP_mm[i],rho)
    ind_bnds_mnths[i] = NetEmPBands(ind_ind,areas,ind_EmP_mm[i],rho)

# Standard deviations for intra-annual E-P latitude bands:
AEmP_mb_sd = pl.std(atl_bnds_mnths,axis=0)
PEmP_mb_sd = pl.std(pac_bnds_mnths,axis=0)
IEmP_mb_sd = pl.std(ind_bnds_mnths,axis=0)

# Monthly E-P-R for latitude bands:
# set up empty arrays for each basin:
atl_EPR_mb = pl.zeros_like(atl_bnds_mnths)
pac_EPR_mb = pl.zeros_like(pac_bnds_mnths)
ind_EPR_mb = pl.zeros_like(ind_bnds_mnths)

for i in range(EmP_mm.shape[0]): # loop over 12 months
    atl_EPR_mb[i] = EPRbands(atl_bnds_mnths[i],atl_ind,atl_run.lat,lat,Arun_mnths[i])
    pac_EPR_mb[i] = EPRbands(pac_bnds_mnths[i],pac_ind,pac_run.lat,lat,Prun_mnths[i])
    ind_EPR_mb[i] = EPRbands(ind_bnds_mnths[i],ind_ind,ind_run.lat,lat,Irun_mnths[i])

# Standard deviations for intra-annual E-P-R latitude bands:
AEPR_mb_sd = pl.std(atl_EPR_mb,axis=0)
PEPR_mb_sd = pl.std(pac_EPR_mb,axis=0)
IEPR_mb_sd = pl.std(ind_EPR_mb,axis=0)

# Plot time series of E-P & E-P-R for latitude bands:
clrs=['b','g','r','k']
fill = ['#B2E6FF','#D1FFA3','#FFC2C2','gainsboro']
trans = [None,0.83,0.675,0.59]
#-Atlantic---------------------------------------------------------------------
atl_labels = ['45$^\circ$N-60$^\circ$N','24$^\circ$N-45$^\circ$N',
              '16$^\circ$S-24$^\circ$N','35$^\circ$S-16$^\circ$S']
fig, ax = pl.subplots(2,1)
# E-P
ax1 = pl.subplot(2,1,1)
for i in range(atl_ind.shape[0]-1):
    pl.plot(pl.linspace(0,11,12),atl_bnds_mnths[:,i],label=atl_labels[i],color=clrs[i])
    pl.fill_between(pl.linspace(0,11,12),atl_bnds_mnths[:,i]-AEmP_mb_sd[i],
                    atl_bnds_mnths[:,i]+AEmP_mb_sd[i],color=fill[i],alpha=trans[i])
pl.axhline(linewidth=1, color='black',ls='--') # add dashed black line at E=P
pl.xticks(pl.linspace(0,11,12))
ax1.xaxis.set_major_formatter(pl.NullFormatter()) # get rid of labels for upper plot x-axis
pl.yticks(pl.linspace(-0.2,1,4))
ax1.set_ylabel('$E-P$ (Sv)',fontsize=16)
pl.title('Atlantic intra-annual $E-P$  variability')


# E-P-R
ax2 = pl.subplot(2,1,2)
for i in range(atl_ind.shape[0]-1):
    pl.plot(pl.linspace(0,11,12),atl_EPR_mb[:,i],label=atl_labels[i],color=clrs[i])
    pl.fill_between(pl.linspace(0,11,12),atl_EPR_mb[:,i]-AEPR_mb_sd[i],
                    atl_EPR_mb[:,i]+AEPR_mb_sd[i],color=fill[i],alpha=trans[i])
pl.axhline(linewidth=1, color='black',ls='--') # add dashed black line at E=P+R
pl.xticks(pl.linspace(0,11,12))
ax2.set_xticklabels(months)
ax2.set_xlabel('Months',fontsize=16)
pl.yticks(pl.linspace(-0.4,0.8,4))
ax2.set_ylabel('$E-P-R$ (Sv)',fontsize=16)
pl.title('Atlantic intra-annual $E-P-R$  variability')
ax2.legend(loc=9,ncol=2,columnspacing=0.35,handletextpad=0.6,fontsize=13)
pl.subplots_adjust(top=0.92)
#------------------------------------------------------------------------------


#-Pacific----------------------------------------------------------------------
pac_labels = ['47$^\circ$N-BS','24$^\circ$N-47$^\circ$N','30$^\circ$S-24$^\circ$N']
fig, ax = pl.subplots(2,1)
# E-P
ax1 = pl.subplot(2,1,1)
for i in range(pac_ind.shape[0]-1):
    pl.plot(pl.linspace(0,11,12),pac_bnds_mnths[:,i],label=pac_labels[i])
    pl.fill_between(pl.linspace(0,11,12),pac_bnds_mnths[:,i]-PEmP_mb_sd[i],
                    pac_bnds_mnths[:,i]+PEmP_mb_sd[i],color=fill[i],alpha=trans[i])
pl.axhline(linewidth=1, color='black',ls='--') # add dashed black line at E=P
ax1.set_ylabel('$E-P$ (Sv)',fontsize=16)
pl.xticks(pl.linspace(0,11,12))
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.yticks(pl.linspace(-0.6,0.8,5))
pl.title('Pacific intra-annual $E-P$  variability')


# E-P-R
ax2 = pl.subplot(2,1,2)
for i in range(pac_ind.shape[0]-1):
    pl.plot(pl.linspace(0,11,12),pac_EPR_mb[:,i],label=pac_labels[i])
    pl.fill_between(pl.linspace(0,11,12),pac_EPR_mb[:,i]-PEPR_mb_sd[i],
                    pac_EPR_mb[:,i]+PEPR_mb_sd[i],color=fill[i],alpha=trans[i])
pl.axhline(linewidth=1, color='black',ls='--')
pl.xticks(pl.linspace(0,11,12))
ax2.set_xticklabels(months)
ax2.set_xlabel('Months',fontsize=16)
ax2.set_ylabel('$E-P-R$ (Sv)',fontsize=16)
pl.ylim(-0.6,0.8)
pl.yticks(pl.linspace(-0.6,0.8,5))
pl.title('Pacific intra-annual $E-P-R$  variability')
ax2.legend(loc=2,ncol=3,fontsize=12,handletextpad=0.5,columnspacing=0.5)
pl.subplots_adjust(top=0.92)
#------------------------------------------------------------------------------


#-Indian-----------------------------------------------------------------------
ind_labels = ['>8$^\circ$S','20$^\circ$S-8$^\circ$S','35$^\circ$S-20$^\circ$S']
fig, ax = pl.subplots(2,1)
# E-P
ax1 = pl.subplot(2,1,1)
for i in range(ind_ind.shape[0]-1):
    pl.plot(pl.linspace(0,11,12),ind_bnds_mnths[:,i],label=ind_labels[i])
    pl.fill_between(pl.linspace(0,11,12),ind_bnds_mnths[:,i]-IEmP_mb_sd[i],
                    ind_bnds_mnths[:,i]+IEmP_mb_sd[i],color=fill[i],alpha=trans[i])
pl.axhline(linewidth=1, color='black',ls='--')
ax1.set_ylabel('$E-P$ (Sv)',fontsize=16)
pl.xticks(pl.linspace(0,11,12))
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.yticks(pl.linspace(-0.8,0.8,5))
pl.title('Indian intra-annual $E-P$  variability')
ax1.legend(loc=3,ncol=3,fontsize=12,handletextpad=0.5,columnspacing=0.5)

# E-P-R
ax2 = pl.subplot(2,1,2)
for i in range(ind_ind.shape[0]-1):
    pl.plot(pl.linspace(0,11,12),ind_EPR_mb[:,i],label=ind_labels[i])
    pl.fill_between(pl.linspace(0,11,12),ind_EPR_mb[:,i]-IEPR_mb_sd[i],
                    ind_EPR_mb[:,i]+IEPR_mb_sd[i],color=fill[i],alpha=trans[i])
pl.axhline(linewidth=1, color='black',ls='--')
ax2.set_ylabel('$E-P-R$ (Sv)',fontsize=16)
pl.xticks(pl.linspace(0,11,12))
ax2.set_xticklabels(months)
ax2.set_xlabel('Months',fontsize=16)
pl.yticks(pl.linspace(-0.8,0.8,5))
pl.title('Indian intra-annual $E-P-R$  variability')
pl.subplots_adjust(top=0.92)
#------------------------------------------------------------------------------



# WRITE THE STANDARD DEVIATIONS TO A FILE

"""f = open('intra-ann_EmP_stdev60.txt','w')
f.write('Standard deviations of monthly means (intra-annual variability) of ' + 
'E-P and E-P-R for the Atlantic (35S-60N), Pacific (30S- Bering Strait) and Indian (>35S)' +
' Oceans. E-P is from ERA-Interim vertically integrated moisture flux divergence' +
' and runoff is from the Dai & Trenberth (2002) dataset. Also presented are ' + 
' standard deviations for various latitude bands within each basin.\n\nAtlantic' + 
' E-P (35S-45N) standard deviation =')
AEmP_mm_sd.tofile(f,sep=" ",format="%.3f")
f.write('\n\nPacific E-P (30S-BS) standard deviation = ')
PEmP_mm_sd.tofile(f,sep=" ",format="%.3f")
f.write('\n\nIndian E-P (>35S) standard deviation = ')
IEmP_mm_sd.tofile(f,sep=" ",format="%.3f")
f.write('\n\nAtlantic E-P bands: \n45N-60N     24N-45N    16S-24N      35S-16S\n')
AEmP_mb_sd.tofile(f,sep="      ",format="%.3f")
f.write('\n\nPacific E-P bands: \n47N-BS     24N-47N      30S-24N\n')
PEmP_mb_sd.tofile(f,sep="      ",format="%.3f")
f.write('\n\nIndian E-P bands: \n>8S    20S-8S    35S-20S\n')
IEmP_mb_sd.tofile(f,sep="   ",format="%.3f")
f.write('\n\n\nAtlantic E-P-R (35S-60N) standard deviation = ')
A_EPR_sd.tofile(f,sep=" ",format="%.3f")
f.write('\n\nPacific E-P-R (30S-BS) standard deviation = ')
P_EPR_sd.tofile(f,sep=" ",format="%.3f")
f.write('\n\nIndian E-P-R (>35S) standard deviation = ')
I_EPR_sd.tofile(f,sep=" ",format="%.3f")
f.write('\n\nAtlantic E-P-R bands: \n45N-60N     24N-45N    16S-24N      35S-16S\n')
AEPR_mb_sd.tofile(f,sep="      ",format="%.3f")
f.write('\n\nPacific E-P-R bands: \n47N-BS     24N-47N    30S-24N\n')
PEPR_mb_sd.tofile(f,sep="      ",format="%.3f")
f.write('\n\nIndian E-P-R bands: \n>8S    20S-8S    35S-20S\n')
IEPR_mb_sd.tofile(f,sep="   ",format="%.3f")
f.close()"""
