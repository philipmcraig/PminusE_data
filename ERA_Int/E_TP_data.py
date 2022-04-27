# -*- coding: utf-8 -*-
"""
Created on Fri May  1 15:14:14 2015

@author: np838619

Last updated: 17/09/15 2:38PM 7th commit
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset, MFDataset
import os
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.colors import Normalize
from matplotlib import ticker
import scipy.ndimage as scifilt
from functions import PrintFiles, MidpointNormalize, SeasonalMeans

# code for reading E and P ERA-Interim and calculating E-P
# what folders are E & P in? monthly_means folder
# hafs files 3 hourly, 3 per month, 96 per year :(
# hafs files gives accumulations from 00 or 12, only 0012 & 1212 are needed
# evaporation and total precipitation
pl.close('all')
# list of years
years = pl.linspace(1979,2014,36)
# need next part to get rid of decimal point
year_input = [str(int(i)) for i in years]

# empty array for filenames
filenames = pl.zeros([len(years),96],dtype='S17') # Only need every 4th file!

#loop over years:
for year in range(len(years)):
    #path = path to mm folder + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + str(year_input[year])
    #filenames[year] = PrintFiles(path,type)
    filenames[year] = PrintFiles(path,'hafs')

lessfiles = pl.zeros([filenames.shape[0],filenames.shape[1]/4],dtype='S17')
for year in range(filenames.shape[0]):
    filelist = []
    for name in range(filenames.shape[1]):
        if '12.nc' in filenames[year,name]:
            filelist.append(filenames[year,name])
    for i in range(len(filelist)):
        lessfiles[year,i] = filelist[i]

#empty array for evaporation
evap = pl.zeros([len(years),lessfiles.shape[1],1,1,256,512])
#empty array for precipitation
prec = pl.zeros_like(evap)

lessfiles = pl.sort(lessfiles)

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
     
# why are there -ve values of TP? Nothing to worry about.

evap_tot = pl.zeros([evap.shape[0],evap.shape[1]/2,evap.shape[2],evap.shape[3],
                     evap.shape[4],evap.shape[5]])
prec_tot = pl.zeros_like(evap_tot)
# loop over number of years:
for year in range(evap.shape[0]):
    # loop over number of months:
    for month in range(int(evap.shape[1]/2)):
        evap_tot[year,month] = -1*(evap[year][2*month] + evap[year][2*month+1])
        prec_tot[year,month] = (prec[year][2*month] + prec[year][2*month+1])

evap_years_mean = pl.mean(evap_tot,axis=1)
prec_years_mean = pl.mean(prec_tot,axis=1)

# calculate the climatological monthly means:
evap_mnths_mean = pl.mean(evap_tot,axis=0)
prec_mnths_mean = pl.mean(prec_tot,axis=0)

# calculate total E & P for each year, along axis=1, total metres of water for each year:
evap_years_sum = pl.sum(30*evap,axis=1)
prec_years_sum = pl.sum(30*prec,axis=1)

evap_mean = pl.mean(evap_years_mean,axis=0) # in m/day
prec_mean = pl.mean(prec_years_mean,axis=0) # in m/day

#calculate mean E & P for each year
#evap_years_mean = pl.mean(evap,axis=1)
#prec_years_mean = pl.mean(prec,axis=1)        

#calculate annual mean E & P:
evap_ann_mean = pl.mean(evap_years_sum,axis=0)
prec_ann_mean = pl.mean(prec_years_sum,axis=0)

#http://www.ecmwf.int/en/how-can-evaporation-have-both-positive-and-negative-values
#Evaporation is -ve values in evap* arrays, there fore E-P = -evap-prec
EmP_ann_mean = (evap_mean - prec_mean)#*365



# Calculate seasonal means:
EmP_mnths = evap_tot - prec_tot
EmP_sns = SeasonalMeans(EmP_mnths)
evap_sns = SeasonalMeans(evap_tot)
prec_sns = SeasonalMeans(prec_tot)

#evap_sns_filt = scifilt.gaussian_filter(evap_sns,sigma=(0,0,0,2.5,1),order=0)
#prec_sns_filt = scifilt.gaussian_filter(prec_sns,sigma=(0,0,0,2.5,1),order=0)
#EmP_sns_filt = scifilt.gaussian_filter(EmP_sns,sigma=(0,0,0,2.5,1),order=0)

#not sure about these results
#plot looks correct but units & scale are confusing

# Land-sea mask:
maskfile = Dataset('/panfs/jasmin/era/era-in/netc/ggis/1989/jan1989/ggis198901010000.nc','r')
mask = maskfile.variables['LSM'][:]
maskfile.close()

evap_masked = pl.zeros_like(prec_mean[0])
prec_masked = pl.zeros_like(evap_mean[0])
EmP_masked = pl.zeros_like(EmP_ann_mean[0])
for i in range(0,EmP_ann_mean[0].shape[1]):
    for j in range(0,EmP_ann_mean[0].shape[2]):
        EmP_masked[0,i,j] = pl.ma.masked_where(mask[0,0,i,j]>0.5, EmP_ann_mean[0,0,i,j])
        # any grid point with mask value >0.5 is NaN, non-ocean values!
        evap_masked[0,i,j] = pl.ma.masked_where(mask[0,0,i,j]>0.5, evap_mean[0,0,i,j])
        prec_masked[0,i,j] = pl.ma.masked_where(mask[0,0,i,j]>0.5, prec_mean[0,0,i,j])
EmP_ocean = pl.ma.masked_invalid(EmP_masked) #changes NaNs to '-'
evap_ocean = pl.ma.masked_invalid(evap_masked)
prec_ocean = pl.ma.masked_invalid(prec_masked)

latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                                    year_input[0] + '/' + filenames[0,0],'r')
lon = latlon.variables['longitude'][:]
lat = latlon.variables['latitude'][:]
latlon.close()

#------------------------WRITE DATA TO NETCDF FILE-----------------------------

#newnc = Dataset('fluxes_ann_mean.nc',mode='w',format='NETCDF4')
#lat_dim = newnc.createDimension('lat', 256)
#lon_dim = newnc.createDimension('lon', 512)
#time_dim = newnc.createDimension('time', None)
#lat_in = newnc.createVariable('lat', pl.float32, ('lat',))
#lat_in.units = 'degrees_north'
#lat_in.long_name = 'latitude'
#lon_in = newnc.createVariable('lon', pl.float32, ('lon',))
#lon_in.units = 'degrees_east'
#lon_in.long_name = 'longitude'
#time = newnc.createVariable('time', pl.float64, ('time',))
##time.units = '' #what should time units be?
#time.long_name = 'time'
#
#EmP = newnc.createVariable('E-P',pl.float64,('time','lat','lon'))
#EmP.units = 'm/day'
#EmP.standard_name = 'evaporation minus precipitation'
#
#lat_in[:] = lat # straight from ERA-Interim nc file
#lon_in[:] = lon # straight from ERA-Interim nc file
#EmP[:,:,:] = EmP_ann_mean
#
#evap_in = newnc.createVariable('E',pl.float64,('time','lat','lon'))
#evap_in.units = 'm/day'
#evap_in.standard_name = 'evaporation'
#evap_in[:,:,:] = evap_mean
#
#prec_in = newnc.createVariable('P',pl.float64,('time','lat','lon'))
#prec_in.units = 'm/day'
#prec_in.standard_name = 'precipitation'
#prec_in[:,:,:] = prec_mean
#
#newnc.close()


#--------------------------------ZONAL AVERAGES--------------------------------


# Over entire globe:
EmP_zonal = pl.mean(EmP_ann_mean*365,axis=-1)
evap_zonal = pl.mean(evap_mean*365,axis=-1)
prec_zonal = pl.mean(prec_mean*365,axis=-1)
pl.figure(1)
pl.plot(lat,pl.squeeze(pl.fliplr(EmP_zonal)),color='k',label='$ \overline{E}- \overline{P}$')
pl.plot(lat,pl.squeeze(pl.fliplr(evap_zonal)),color='red',label='$\overline{E}$')
pl.plot(lat,pl.squeeze(pl.fliplr(prec_zonal)),color='blue',label='$\overline{P}$')
pl.axhline(y=0,color='k',ls='--')
pl.xlim(-90,90)
pl.xticks([-90,-60,-30,0,30,60,90])
pl.legend(loc=0)
pl.ylabel('m/yr')
pl.xlabel('Latitude ($^\circ$)')
pl.title('Zonally averaged freshwater fluxes')

# Over global oceans:

EmP_ocz = pl.nanmean(EmP_ocean*365,axis=-1)
evap_ocz = pl.nanmean(evap_ocean*365,axis=-1)
prec_ocz = pl.nanmean(prec_ocean*365,axis=-1)
pl.figure(2)
pl.plot(lat,pl.squeeze(EmP_ocz),color='k',label='$ \overline{E}- \overline{P}$')
pl.plot(lat,pl.squeeze(evap_ocz),color='red',label='$\overline{E}$')
pl.plot(lat,pl.squeeze(prec_ocz),color='blue',label='$\overline{P}$')
pl.axhline(y=0,color='k',ls='--')
pl.ylim(-1.5,3)
pl.xlim(-90,90)
pl.xticks([-90,-60,-30,0,30,60,90])
pl.legend(loc=0)
pl.ylabel('m/yr')
pl.xlabel('Latitude ($^\circ$)')
pl.title('Zonally averaged oceanic freshwater fluxes')

#----------------------------PLOT ANNUAL MEAN E-P------------------------------

#latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
#                                    year_input[0] + '/' + filenames[0,0],'r')
#lon = latlon.variables['longitude'][:]
#lat = latlon.variables['latitude'][:]

# can't remember why these are here, need to check
lon_0 = lon.mean()
lat_0 = lat.mean()
pl.figure(3)
m = Basemap(projection='cyl',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
        llcrnrlon=-180,urcrnrlon=180,lat_ts=20)
EmP_ann_mean, lon = shiftgrid(180.0, EmP_ann_mean*365, lon, start=False)
m.drawcoastlines()
lons, lats = pl.meshgrid(lon,lat)
X, Y = m(lons,lats)

#norm = MidpointNormalize(midpoint=0)

cmap = pl.get_cmap('seismic')
#cmap.set_over('navy')
#cs = m.pcolormesh(X, Y, pl.squeeze(EmP_ann_mean), cmap=cmap, norm=pl.Normalize(-4,4,clip=False))
cf = m.contourf(X,Y,pl.squeeze(EmP_ann_mean),cmap=cmap,norm=pl.Normalize(-4,4,clip=False),
                levels=[-4,-2,-1,-0.5,-0.25,-0.125,0.125,0.25,0.5,1,2,4],extend="min")
                

cbar = m.colorbar(cf,location='bottom', pad = "10%",ticks=[-4,-1,-0.25,0,0.25,1,4],extend="min")
cbar.set_label('m yr$^{-1}$',fontsize=16)
#tick_locator = ticker.MaxNLocator(nbins=5)
#cbar.locator = tick_locator
#cbar.update_ticks()
#m.fillcontinents(color='gainsboro')
m.drawparallels([-60,-30,0,30,60],labels=[1,0,0,1],linewidth=0.5,ax=None,
                fontsize=10,dashes=[1,2])
m.drawmeridians([-90,0,90],labels=[1,0,0,1],linewidth=0,ax=None,fontsize=10,
                dashes=[1,2])
pl.title('ERA-Interim Annual Mean $E-P$  1979-2014',fontsize=16)
#pl.savefig('EmP_direct_1979-2014.png')


#---------------------------PLOT SEASONAL MEAN---------------------------------

sn_names = ['DJF','MAM','JJA','SON']
panels = ['(a)','(b)','(c)','(d)']

fig, ax = pl.subplots(2,2,figsize=(24,13))
for i in range(EmP_sns.shape[0]):
    latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                                year_input[0] + '/' + filenames[0,0],'r')
    lon = latlon.variables['longitude'][:]
    lat = latlon.variables['latitude'][:]
    latlon.close()
    lon_0 = lon.mean()
    lat_0 = lat.mean()
    pl.subplot(2,2,i+1)#,sharex=fig,sharey=fig) # subplot starts at 1 hence i+1 needed to index subplots
    m = Basemap(projection='cyl',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
                    llcrnrlon=-180,urcrnrlon=180,lat_ts=20)
    prec_sns[i], newlon = shiftgrid(180.0, prec_sns[i]*1000, lon, start=False)
    m.drawcoastlines(linewidth=0.5)
    lons, lats = pl.meshgrid(newlon,lat)
    X, Y = m(lons, lats)
    #norm = MidpointNormalize(midpoint=0) # calls a class that centres the colourmap at zero
    cmap = pl.get_cmap('YlGnBu')#YlOrRd
    #cmap.set_bad('lightgrey')
    # constrain all data within the same bounds to use use same colour map:
    #pcm=pl.pcolormesh(X,Y,pl.squeeze(EmP_sns[i]),cmap=cmap,norm=pl.Normalize(-10,10,clip=False))
    cf = m.contourf(X,Y, pl.squeeze(prec_sns[i]),cmap=cmap,norm=pl.Normalize(0,10,clip=False),
                levels=[0,2,4,6,8,10],extend="max")
    m.fillcontinents(color='gainsboro')
    m.drawparallels([-60,-30,0,30,60],labels=[0,0,0,0],linewidth=0.5,
                        ax=None,fontsize=22,dashes=[1,2])
    m.drawmeridians([-90,0,90],labels=[0,0,0,0],linewidth=0.5,ax=None,
                        fontsize=22,dashes=[1,2])
    if i in (0,2):
        m.drawparallels([-60,-30,0,30,60],labels=[1,0,0,1],linewidth=0.0,
                        ax=None,fontsize=22)
    if i > 1:
        m.drawmeridians([-90,0,90],labels=[1,0,0,1],linewidth=0.0,ax=None,
                        fontsize=22)
    pl.title(panels[i]+'  '+sn_names[i],fontsize=25)
    

f = pl.gcf()
colax = f.add_axes([0.9,0.25,0.015,0.5]) # set position of colour bar
#pcmbounds = pl.linspace(-20,20,11)
#cmap.set_over(color='maroon')
#cmap.set_under(color='#000047')
clb = pl.colorbar(cf, cax=colax,extend='both')#,boundaries=pcmbounds)
clb.set_label('$e-p$  (mm day$^{-1}$)',fontsize=30)
clb.ax.tick_params(labelsize=20) 
pl.subplots_adjust(wspace=0.04,hspace=-0.02,left=0.12,right=0.88,top=0.90,bottom=0.1)
pl.savefig('/home/np838619/PminusE_data/ERA_Int/plots/prec_sns_cont.png')