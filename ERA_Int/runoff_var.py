# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 09:57:05 2015

@author: np838619

Code to calculate the intra-annual variability of runoff and produce time series
of the seasonal cycles for each basin.

Last updated: 29/09/2015 4:58PM 2nd commit
"""

# import stuff up here: pylab, functions, division, pandas
from __future__ import division
import pylab as pl
import pandas
from netCDF4 import Dataset
from functions import *

# need the ERA-Interim latitude array:
ncfile = Dataset('tcdq_mean.nc','r')
lat = ncfile.variables['lat'][:]
ncfile.close()

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

# specify latitudes at which to restrict runoff:
atl_lats = pl.array([60,45,24,-16,-35])
pac_lats = pl.array([65.5,47,24,-30])
ind_lats = pl.array([73.5,-8,-20,-35])

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

#set up empty arrays for monthly runoff totals:
atl_totals = pl.zeros([Arun_mnths.shape[0]])
pac_totals = pl.zeros_like(atl_totals)
ind_totals = pl.zeros_like(atl_totals)
# calculate the total runoff for each basin for each month:
for m in range(Arun_mnths.shape[0]):
    atl_totals[m] = TotalRunoff(atl_run.lat,lat,(atl_ind[0],atl_ind[-1]),Arun_mnths[m])
    pac_totals[m] = TotalRunoff(pac_run.lat,lat,(pac_ind[0],pac_ind[-1]),Prun_mnths[m])
    ind_totals[m] = TotalRunoff(ind_run.lat,lat,(ind_ind[0],ind_ind[-1]),Irun_mnths[m])

# calculate standard deviations
atl_sd = pl.std(atl_totals)
pac_sd = pl.std(pac_totals)
ind_sd = pl.std(ind_totals)

# plot as single figure
mnth_lbls = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
fig, ax = pl.subplots()
pl.plot(pl.linspace(0,11,12),atl_totals,color='blue',label='Atlantic')
pl.plot(pl.linspace(0,11,12),pac_totals,color='green',label='Pacific')
pl.plot(pl.linspace(0,11,12),ind_totals,color='red',label='Indian')
pl.xticks(pl.linspace(0,11,12))
ax.set_xticklabels(mnth_lbls)
ax.set_xlabel('months',fontsize=17)
pl.ylim=(0,1.0)
ax.set_ylabel('$R$  (Sv)',fontsize=17)
pl.title('Intra-annual variability of runoff',fontsize=20)
ax.legend(loc=0)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/runoff_basins_intra.png')


# calculate total monthly runoff for latitude bands:

# set up empty arrays for each basin:
atl_bnds_tot = pl.zeros([Arun_mnths.shape[0],atl_ind.shape[0]-1])
pac_bnds_tot = pl.zeros([Prun_mnths.shape[0],pac_ind.shape[0]-1])
ind_bnds_tot = pl.zeros([Irun_mnths.shape[0],ind_ind.shape[0]-1])

for m in range(Arun_mnths.shape[0]):
    for b in range(atl_bnds_tot.shape[1]):
        atl_bnds_tot[m,b] = TotalRunoff(atl_run.lat,lat,(atl_ind[b],atl_ind[b+1]),Arun_mnths[m])
    for b in range(pac_bnds_tot.shape[1]):
        pac_bnds_tot[m,b] = TotalRunoff(pac_run.lat,lat,(pac_ind[b],pac_ind[b+1]),Prun_mnths[m])
        ind_bnds_tot[m,b] = TotalRunoff(ind_run.lat,lat,(ind_ind[b],ind_ind[b+1]),Irun_mnths[m])
        
# Plot as 3x1 subplot:
clrs = ['blue','green','red','k']
fig, ax = pl.subplots(3,1,figsize=(10,10))
ax1 = pl.subplot(311)
atl_lbls = ['45$^\circ$N-60$^\circ$N','24$^\circ$N-45$^\circ$N','16$^\circ$S-24$^\circ$N',
            '35$^\circ$S-16$^\circ$S']
for i in range(atl_bnds_tot.shape[1]):
    pl.plot(pl.linspace(0,11,12),atl_bnds_tot[:,i],color=clrs[i],label=atl_lbls[i])
pl.xticks(pl.linspace(0,11,12))
ax1.xaxis.set_major_formatter(pl.NullFormatter())
ax1.set_ylabel('$R$  (Sv)',fontsize=17)
ax1.legend(loc=0)
pl.title('Atlantic Ocean',fontsize=18)

ax2 = pl.subplot(312)
pac_lbls = ['47$^\circ$-BS','24$^\circ$-47$^\circ$','30$^\circ$S-24$^\circ$N']
for i in range(pac_bnds_tot.shape[1]):
    pl.plot(pl.linspace(0,11,12),pac_bnds_tot[:,i],color=clrs[i],label=pac_lbls[i])
pl.xticks(pl.linspace(0,11,12))
ax2.xaxis.set_major_formatter(pl.NullFormatter())
ax2.set_ylim([0,0.6])
ax2.set_ylabel('$R$  (Sv)',fontsize=17)
ax2.legend(loc=0)
pl.title('Pacific Ocean',fontsize=18)

ax3 = pl.subplot(313)
ind_lbls = ['>8$^\circ$S','20$^\circ$S-8$^\circ$S','35$^\circ$-20$^\circ$S']
for i in range(ind_bnds_tot.shape[1]):
    pl.plot(pl.linspace(0,11,12),ind_bnds_tot[:,i],color=clrs[i],label=ind_lbls[i])
pl.xticks(pl.linspace(0,11,12))
ax3.set_xticklabels(mnth_lbls)
ax3.set_xlabel('months',fontsize=17)
ax3.set_ylim([0.,0.6])
ax3.set_ylabel('$R$  (Sv)',fontsize=17)
ax3.legend(loc=0)
pl.title('Indian Ocean',fontsize=18)

fig.suptitle('Runoff seperated into latitude bands',fontsize=22)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/runoff_bnds_intra.png')