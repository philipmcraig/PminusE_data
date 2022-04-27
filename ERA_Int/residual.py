# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 12:12:33 2015

@author: np838619
Code to calculate residual between moisture flux divergence and E-P from ERA-Interim

Last updated: 3/3/2016 1:23PM Probably 2nd commit

"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.colors import Normalize
from matplotlib import ticker
from functions import MidpointNormalize
import scipy.ndimage.filters as scifilt
from functions import *

def WtdStd(weights,data):
    """
    """
    wtd_av = pl.average(data,weights=weights) # weighted average
    a = pl.sum(weights*(data-wtd_av)**2)
    b = pl.sum(weights)
    zeros = pl.where(weights==0)
    nonzero = weights.size - zeros[0].shape[0]
    std = pl.sqrt(a/(((nonzero-1)/nonzero)*b))
    
    return std


def GlobalMap(lon,lat,inarray,levels,norm,cmap,ticks,label,extend):
    """
    """
    
    m = Basemap(projection='cyl',resolution='l',llcrnrlat=-80.,urcrnrlat=80.,
                llcrnrlon=-180.,urcrnrlon=180.,lon_0=lon.mean(),lat_0=lat.mean(),
                lat_ts=20)
    inarray, lons = shiftgrid(180.0,inarray,lon,start=False)
    m.drawcoastlines()
    lons, lats = pl.meshgrid(lons,lat)
    X,Y = m(lons,lats)
    
    cf = m.contourf(X,Y,pl.squeeze(inarray),levels=levels,norm=norm,cmap=cmap,extend=extend)
    cbar = m.colorbar(cf,location='bottom',pad='10%',ticks=ticks)
    cbar.set_label(label,fontsize=16)
    
    m.drawparallels([-60,-30,0,30,60],labels=[1,0,0,1],linewidth=0,ax=None)
    m.drawmeridians([-90,0,90],labels=[1,0,0,1],linewidth=0,ax=None)
    
    return m


# read nc files
ncfile = Dataset('/home/np838619/PminusE_data/ERA_Int/flux7914.nc','r') # file with annual mean moisture budget
# extract lon & lat from one of the files:
lon = ncfile.variables['lon'][:]
lat = ncfile.variables['lat'][:]
tcdq = ncfile.variables['tcdq'][:]
EmP = ncfile.variables['E-P'][:]
# close the files:
ncfile.close()

nc2 = Dataset('/home/np838619/PminusE_data/ERA_Int/tcdq_36years_means.nc','r')
tcdq36 = nc2.variables['tcdq'][:]
nc2.close()

nc3 = Dataset('/home/np838619/PminusE_data/ERA_Int/evapprec_36yrs.nc','r')
evap36 = nc3.variables['evap'][:]
prec36 = nc3.variables['prec'][:]
nc3.close()

nc4 = Dataset('/home/np838619/PminusE_data/ERA_Int/TCWV.nc','r')
tcwv = nc4.variables['tcwv'][:]
nc4.close()

basins = Dataset('/home/np838619/PminusE_data/ERA_Int/tcdq_basins_60.nc','r')
atlmask = basins.variables['maskA'][:]
balmask = basins.variables['maskB'][:]
medmask = basins.variables['maskM'][:]
pacmask = basins.variables['maskP'][:]
indmask = basins.variables['maskI'][:]
basins.close()

maskfile = Dataset('tcdq_ann_mean_ERAI.nc','r')
eramask = maskfile.variables['LSM'][:]
maskfile.close()

atlmask = atlmask + balmask + medmask # add Atlantic & Med masks together
for i in range(atlmask.shape[0]):
    for j in range(atlmask.shape[0]):
        if atlmask[i,j] > 1.:
            atlmask[i,j] = 1.

#EmPfile = Dataset('direct_EmP.nc','r') # file with annual mean E-P
#EmP_globe = EmPfile.variables['E-P'][:]
#EmP_ocean = EmPfile.variables['E-P ocean'][:]
#mask = EmPfile.variables['LSM'][:]
#EmPfile.close()

# smooth the divergence with a gaussian filter:
tcdq_smth = scifilt.gaussian_filter(tcdq,sigma=(2,1),order=0)
EmP_smth = scifilt.gaussian_filter(EmP,sigma=(2,1),order=0)

zonal_smth=pl.mean(tcdq_smth,axis=1)
#pl.plot((pl.squeeze(gf_zonal)/(10**3))*86400*365)
#pl.figure(2)
#pl.imshow(pl.squeeze(gf))

#gf_masked = pl.zeros_like(gf)
#for i in range(0,gf.shape[1]):
#    for j in range(0,gf.shape[2]):
#        gf_masked[0,i,j] = pl.ma.masked_where(mask[0,i,j]>0., gf[0,i,j])
#        # any tcdq value with mask value >0 is NaN, non-ocean values!
#gf_ocean = pl.ma.masked_invalid(gf_masked)


# divergence - E-P
residual = tcdq*86400 - EmP*1000
#residual=tcdq_ocean
res_smth = tcdq_smth*86400 - EmP_smth*1000

pl.figure(1)
GlobalMap(lon,lat,res_smth,[-4,-2,-1,-0.5,-0.25,-0.125,0.125,0.25,0.5,1,2,4],
          pl.Normalize(-5,5),'seismic',[-4,-1,-0.25,0,1,0.25,4],'mm/day','both')
#pl.title('Moisture flux divergence minus $E-P$ (1979-2014)',fontsize=20)
#pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/residual.png')
pl.close()


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

#res_wtd = residual*areas # weight each residual value by grid cell area
#res_gav = pl.sum(res_wtd)/pl.sum(areas) # globally averaged residual
res_gav = pl.average(residual,weights=areas)
print res_gav

# Calculate the same weights that Paul Berrisford used: d(sin(lat))*d(lon)
Q = pl.zeros_like(residual)
for i in range(Q.shape[0]):
    for j in range(Q.shape[1]):
        Q[i,j] = pl.absolute(pl.sin(lat_half[i+1]-pl.sin(lat_half[i])))*delta_lambda

res_std = WtdStd(Q,residual)
print res_std

# Residual for each basin:

atlres = atlmask*residual*(1-eramask[0]); atl_area = atlmask*areas*(1-eramask[0])
pacres = pacmask*residual*(1-eramask[0]); pac_area = pacmask*areas*(1-eramask[0])
indres = indmask*residual*(1-eramask[0]); ind_area = indmask*areas*(1-eramask[0])

atl_av = pl.average(atlres,weights=atl_area)
pac_av = pl.average(pacres,weights=pac_area)
ind_av = pl.average(indres,weights=ind_area)

# tcwv for each basin:

atlwv = atlmask*tcwv*(1-eramask[0])
pacwv = pacmask*tcwv*(1-eramask[0])
indwv = indmask*tcwv*(1-eramask[0])

atlwv_av = pl.zeros([36]); pacwv_av = pl.zeros([36]); indwv_av = pl.zeros([36])

for yr in range(0,36):
    atlwv_av[yr] = pl.average(atlwv[yr],weights=atl_area)
    pacwv_av[yr] = pl.average(pacwv[yr],weights=pac_area)
    indwv_av[yr] = pl.average(indwv[yr],weights=ind_area)

rho = 10**3
# Annual mean residuals in Sv:
atl_net = NetEmPCalc(atl_area,atlres*(rho/(86400*1000)),rho); print 'Atlantic residual = ', atl_net, ' Sv'
pac_net = NetEmPCalc(pac_area,pacres*(rho/(86400*1000)),rho); print 'Pacific residual = ', pac_net, ' Sv'
ind_net = NetEmPCalc(ind_area,indres*(rho/(86400*1000)),rho); print 'Indian residual = ', ind_net, ' Sv'


# Calculate residual for each year:

EmP36 = evap36 - prec36

res36 = tcdq36 - EmP36*(1000/86400)
tcdq36_smth = scifilt.gaussian_filter(tcdq36,sigma=(0,2,1),order=0)
EmP36_smth = scifilt.gaussian_filter(EmP36,sigma=(0,2,1),order=0)
res36_smth = tcdq36_smth*86400 - EmP36_smth*1000

res36_std = pl.std(res36,axis=0) # temporal standard deviation of each point

#pl.figure(2)
#tx = [0,0.001,0.005,0.1,0.5,1]
#GlobalMap(lon,lat,res36_std,tx,
#              pl.Normalize(0,1.5),'coolwarm',tx,'mm/day','max')
#pl.close()

years = pl.linspace(1979,2014,36)


# Calculate globally averaged residual for each year:
res_yrs = pl.zeros_like(years) # empty array for each year's residual
for i in range(len(res_yrs)):
    res_yrs[i] = pl.average(res36[i],weights=areas)
    #print str(int(years[i])) + ' residual = ', res_yrs[i]

# Plot each res36 and save to create gif later:
#years = pl.linspace(1979,2014,36)
#pl.figure(3)
#for year in range(0,36):    
#    GlobalMap(lon,lat,res36_smth[i],[-4,-2,-1,-0.5,-0.25,-0.125,0.125,0.25,0.5,1,2,4],
#              pl.Normalize(-5,5),'seismic',[-4,-1,-0.25,0,1,0.25,4],'mm/day','both')
#    pl.title(str(int(years[year])),fontsize=20)
#    print 'Plotting ' + str(int(years[year]))
#    #pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/residual'+str(int(years[year]))+'.png')
#    if year != 35:
#        pl.clf()


# Calculate basin residuals for each year:

atl36 = pl.zeros_like(res_yrs)
pac36 = pl.zeros_like(atl36)
ind36 = pl.zeros_like(atl36)


for i in range(len(atl36)):
    atl36[i] = NetEmPCalc(atl_area,res36[i]*atlmask*(1-eramask[0]),rho)
    pac36[i] = NetEmPCalc(pac_area,res36[i]*pacmask*(1-eramask[0]),rho)
    ind36[i] = NetEmPCalc(ind_area,res36[i]*indmask*(1-eramask[0]),rho)
    
atl_sd = pl.std(atl36)
pac_sd = pl.std(pac36)
ind_sd = pl.std(ind36)


atl_lats = pl.array([45,-35])#pl.array([60,45,24,-16,-35])#
pac_lats = pl.array([47,-30])#pl.array([65.5,47,24,-17,-30])#
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

atl_bnds36 = pl.zeros([len(years),len(atl_ind)-1])
pac_bnds36 = pl.zeros([len(years),len(pac_ind)-1])
ind_bnds36 = pl.zeros([len(years),len(ind_ind)-1])

for i in range(atl_bnds36.shape[0]):
    atl_bnds36[i] = NetEmPBands(atl_ind,atl_area,res36[i]*atlmask*(1-eramask[0]),rho)
atl_bnds_sd = pl.std(atl_bnds36,axis=0)

for i in range(pac_bnds36.shape[0]):
    pac_bnds36[i] = NetEmPBands(pac_ind,pac_area,res36[i]*pacmask*(1-eramask[0]),rho)
pac_bnds_sd = pl.std(pac_bnds36,axis=0)

for i in range(ind_bnds36.shape[0]):
    ind_bnds36[i] = NetEmPBands(ind_ind,ind_area,res36[i]*indmask*(1-eramask[0]),rho)
ind_bnds_sd = pl.std(ind_bnds36,axis=0)


# Per unit area flux, cm/yr
scl_fctr = (10**6)*(100*365*86400)
atl36_scl = FluxScaled(atl_area,atlmask*(1-eramask[0]),atl36)*scl_fctr
pac36_scl = FluxScaled(pac_area,pacmask*(1-eramask[0]),pac36)*scl_fctr
ind36_scl = FluxScaled(ind_area,indmask*(1-eramask[0]),ind36)*scl_fctr
# multiply by 10/365 for mm/day

fig, ax = pl.subplots()
pl.plot(pl.linspace(0,35,36),atl36_scl,label='Atlantic')
pl.plot(pl.linspace(0,35,36),pac36_scl,label='Pacific')
pl.plot(pl.linspace(0,35,36),ind36_scl,label='Indian')
pl.axhline(y=0,ls='--',color='k')
ax.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=12)
ax.set_xlabel('Years',fontsize=16)
ax.set_ylabel('Residual (cm/yr)',fontsize=16)
ax.legend()
#pl.close()

resbas = zip(atl36,pac36,ind36)

# WRITE AREA-AVERAGED BASIN RESIDUALS TO FILE:
#f = open('/home/np838619/PminusE_data/ERA_Int/res_basins_yrs.txt','w')
#
#f.write('Area-averaged annual mean (1979-2014) ERA-Interim moisture budget residuals'
#+ '(divQ-(E-P)) for each ocean basin.\n\n')
#f.write('Atlantic      Pacific   Indian\n')
#pl.savetxt(f,resbas,fmt='%9.3f')
#f.close()