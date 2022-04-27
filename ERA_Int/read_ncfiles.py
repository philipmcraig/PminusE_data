"""Code to read a load of netcdf files and do stuff with vertically integrated 
moisture flux divergence and total column water vapour. Plots loads of stuff 
too.

Last updated: 17/09/2015 2:08PM 19th commit
"""


from __future__ import division
import pylab as pl
from netCDF4 import Dataset, MFDataset
import os
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.colors import Normalize
from matplotlib import ticker
import scipy.ndimage as scifilt
from functions import *

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
    """
    list_of_files = os.listdir(path_to_files)
    filelist = []
    for each_file in list_of_files:
        if each_file.startswith(start_of_file):
            filelist.append(each_file)

    return filelist

def TryDirectories(path_to_folders,years):
    """
    """
    list_of_dirs = os.listdir(path_to_folders)
    dirlist = []
    for each_folder in list_of_dirs:
        #if name of folder is in years
        if each_folder in years:
            dirlist.append(each_folder)
    
    # now how to loop over all files?
    # filelist_'each_folder' possible? NO
    #loop over dirlist
    #for each_folder in dirlist:                    

    # Perhaps this would work best as  2 functions?

    return dirlist

def Looping(filetype,main_directory):
    """This function wasn't used. It may have potential though.
    path to each file broken down into each string, possibily giving the chance
    to loop over years, months and times.
    """
    path = main_directory + filetype + '/' + year_input + '/' + month + '/' + time + '/' + '.nc'
    
def LatLon(filepath):
        """
        """
        ncfile =  Dataset(filepath,'r')
        lat = ncfile.variables['longitude'][:]
        lat = ncfile.variables['latitude'][:]
        
        lon_0 = lon.mean()
        lat_0 = lat.mean()
        
        return lat, lon, lat_0, lon_0
        

#def main():
   
"""# call PrintFiles function with directory path and beginning of file names
filenames = PrintFiles('/panfs/jasmin/era/era-in/netc/monthly_means/2014/','hgaw')

# need to take latitude and longitude from ONE FILE ONLY:
# could probably be done in the for loop using an if statement as well
latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/2014/' + filenames[0],'r')
lon = latlon.variables['longitude'][:]
lat = latlon.variables['latitude'][:]

# can't remember why these are here, need to check
lon_0 = lon.mean()
lat_0 = lat.mean()

#tcdq_array = pl.zeros([len(filenames)])
# define multi-dimensional array to hold data
# 48 files, then other dimensions are the dimensions of the variable
# would probably be better to get these straight from the file
tcdq_array = pl.ndarray(shape=(48, 1, 1, 256, 512))
#textfile = open('list_of_ncfiles_2014.txt','w') #for writing file names to
for i in range(len(filenames)):
    # read netcdf file, add each file name to directory path
    ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/2014/'+ filenames[i],'r')
    # stick the variable array into the multi-dimensional array:
    tcdq_array[i] = ncfile.variables['TCDQ'][:]
    ncfile.close() # closing the netcdf file saves memory!
    #textfile.write(i+'\n') #writes file names to text file
    #print i
    #textfile.close()

# Calculate total of variable along 0 axis:        
total = pl.sum(tcdq_array,axis=0)
# Calculate mean of variable:
mean = total/48
    
m = Basemap(projection='merc',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
        llcrnrlon=-180,urcrnrlon=180,lat_ts=20)
mean, lon = shiftgrid(180.0, mean, lon, start=False)
m.drawcoastlines()
lon, lat = pl.meshgrid(lon,lat)
X, Y = m(lon,lat)

norm = MidpointNormalize(midpoint=0)

cmap = pl.get_cmap('seismic')
cs = m.pcolormesh(X, Y, pl.squeeze(mean), cmap=cmap, norm=norm)

cbar = m.colorbar(cs,location='bottom', pad = "10%")
cbar.set_label('kg m$^{-2}$ s$^{-1}$')
tick_locator = ticker.MaxNLocator(nbins=3)
cbar.locator = tick_locator
cbar.update_ticks()
pl.title('Mean vertical integral of divergence of moisture flux 2014')
#pl.savefig('tcdq_mean_2014.png')
pl.show()"""

#dirs_and_files = TryDirectories('/panfs/jasmin/era/era-in/netc/monthly_means/')
#for j in dirs_and_files:
 #   print j

years = pl.linspace(1979,2014,36)
year_input = [str(int(i)) for i in years]

#months = pl.linspace(1,12,12)
#month_input = [str(int(j)) for j in months]

#times = pl.linspace(0,18,4)
#time_input = [str(int(k)).zfill(2) for k in times]


# need empty ndarray, shape: len(year_input)
#filenames = pl.ndarray(shape=(len(year_input),48))
filenames = pl.zeros([len(year_input),48],dtype='S17')
# dimensions of filenames are numbers of years and number of files
# loop over year_input:
"""for year in range(len(year_input)):
    #path to montly_means + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + year_input[year]
    ## call PrintFiles function
    filenames[year] = PrintFiles(path,'hgaw')
    
# now netcdf files must be read:
# make an empty array for all the tcdq data:
tcdq_years = pl.ndarray(shape=(len(year_input),len(filenames[0]),1,1,256,512))
tcwv_years = pl.zeros_like(tcdq_years)
#tcdq_years = pl.zeros([len(year_input),len(filenames[0]),1,1,256,512])
# loop over years:
for year in range(len(year_input)):
    #loop over rows in filenames:
    #for row in range(len(filenames[:,1])):
    # loop over filenames:
    for name in range(len(filenames[year])):
        # read netcdf file, path + year + filename
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                    str(year_input[year]) + '/' + str(filenames[year,name]),'r')
        # stick variable array into the multi-dimensional array
        tcdq_years[year,name] = ncfile.variables['TCDQ'][:]
        tcwv_years[year,name] = ncfile.variables['TCWV'][:]
        #close file
        ncfile.close()"""
        
# that loop goes over too many files, turns out the monthly means were already there as ggaw files
# making the code significantly more efficient:

filenames = pl.zeros([len(year_input),12],dtype='S13')
for year in range(len(year_input)):
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + year_input[year]
    filenames[year] = PrintFiles(path,'ggaw')
    
tcdq_years = pl.zeros([len(year_input),len(filenames[0]),1,1,256,512])
tcwv_years = pl.zeros_like(tcdq_years)
for year in range(len(year_input)):
    for name in range(len(filenames[year])):
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                    str(year_input[year]) + '/' + str(filenames[year,name]),'r')
        tcdq_years[year,name] = ncfile.variables['TCDQ'][:]
        tcwv_years[year,name] = ncfile.variables['TCWV'][:]
        ncfile.close()
        
tcdq_years = pl.squeeze(tcdq_years)

tcdq97 = tcdq_years[18,5:]; tcdq98 = tcdq_years[19,:5]
EN = pl.zeros([12,256,512])
EN[:7] = tcdq97[:]; EN[7:] = tcdq98[:]
ENmean = pl.mean(EN,axis=0)
pl.imshow(ENmean,cmap='seismic',norm=pl.Normalize(-0.0002,0.0002))

tcdq88 = tcdq_years[9,5:]; tcdq89 = tcdq_years[10,:5]
LN = pl.zeros([12,256,512])
LN[:7] = tcdq88[:]; LN[7:] = tcdq89[:]
LNmean = pl.mean(LN,axis=0)
pl.imshow(LNmean,cmap='seismic',norm=pl.Normalize(-0.0002,0.0002))


# may want to smooth data at this point with:
# years_smoothed = scipy.ndimage.gaussian_filter(tcdq_years, sigma=1.0, order=0)
# can continue with calculations as normal after, but plotting doesn't seem to work
# shiftgrid doesn't seem to work on smoothed data

tcdq_years_sum = pl.sum(30*tcdq_years,axis=1)
tcwv_years_sum = pl.sum(30*tcwv_years,axis=1)

# WRITE EACH YEAR TO netCDF file
"""all_years = Dataset('TCWV.nc','w')

lat_dim = all_years.createDimension('lat', 256)
lon_dim = all_years.createDimension('lon', 512)
lat_in = all_years.createVariable('lat', pl.float32, ('lat',))
lat_in.units = 'degrees_north'
lat_in.long_name = 'latitude'
lon_in = all_years.createVariable('lon', pl.float32, ('lon',))
lon_in.units = 'degrees_east'
lon_in.long_name = 'longitude'

#time_dim = all_years.createDimension('time',36)
#time = all_years.createVariable('time', pl.float64, ('time',))
#time.units = 'years'
#time.long_name = 'time'

TCWV = all_years.createVariable('tcwv',pl.float64,('time,'lat','lon'))
TCWV.units = 'kg m**-2'
TCWV.standard_name = 'total column water vapour'
lat_in[:] = lat # straight from ERA-Interim nc file
lon_in[:] = lon # straight from ERA-Interim nc file
TCWV[:,:,:] = tcwv_years_mean

all_years.close()"""

# calculate mean of variable for each year
tcdq_years_mean = pl.mean(tcdq_years,axis=1) # axis 1 is months axis
tcwv_years_mean = pl.mean(tcwv_years,axis=1)
# now each array along the 0 axis of years_mean is mean of variable for 1 year
# now the mean of variable for any year could be plotted

tcdq_mean = pl.mean(tcdq_years_mean,axis=0) # tcdq annual mean directly from monthly means
#tcwv_mean = pl.mean(tcwv_years_mean,axis=0)

# calculate the mean for entire re-analysis period:
tcdq_ann_mean = pl.mean(tcdq_years_sum,axis=0) # tcdq annual mean from tcdq annual TOTALS
tcwv_ann_mean = pl.mean(tcwv_years_sum,axis=0)

tcdq_filt = scifilt.gaussian_filter(tcdq_ann_mean,sigma=(2.5,1),order=0)
tcdq_filt2 = scifilt.gaussian_filter(tcdq_mean,sigma=(2.5,1),order=0)

#----------------------------E-P MONTHLY MEANS-----------------------------

#monthly means, calculate along tcdq_years 0 axis to get means for all 12 months:
# need to take into account different length of each month
months_lenths = [31,28,31,30,31,30,31,31,30,31,30,31]
# apply Gaussian Filter to each year:
tcdq_yrs_filt = scifilt.gaussian_filter(tcdq_years,sigma=(0,0,2.5,1),order=0)
tcdq_mm = pl.mean(tcdq_yrs_filt,axis=0)
tcwv_mm = pl.mean(tcwv_years,axis=0)

# need centred difference to calculate d/dt of monthly mean TCWV

# need empty array for time derivatives
#tcwv_months_dt = pl.zeros_like(tcwv_mm)

#specify 2*dt, no. of seconds in a month x2?
time_step = (365/12)*86400

# loop over months, from 0 to 10
#for m in range(tcwv_monthly_means.shape[0]-1):
    #d/dt[month] = (mean[month+1]-mean[month-1])/2*dt
#    tcwv_months_dt[m] = (tcwv_monthly_means[m+1] - tcwv_monthly_means[m-1])/(2*time_step)
#calculate time derivative for December:
#d/dt[-1] = (mean[0]-mean[-2])/2*dt
#tcwv_months_dt[-1] = (tcwv_monthly_means[0] - tcwv_monthly_means[-2])/(2*time_step)

tcwv_dt = dtcwv_dt(tcwv_mm)

# Calculate E-P monthly means: d/dt(tcwv) + tcdq

EmP_mm = (tcwv_dt + tcdq_mm)*86400
# rescaling by a factor of 86400/30 gives E-P in mm/day
#not completely accurate as not every month has 30 days, will need to split it up eventually
#probably have an array with 30,31 and 28 then loop over each months
# would now like to mask this


#Write this array to a netCDF file?
# need to apply land-sea mask at some point, don't care about land!

# Land-sea mask:
maskfile = Dataset('/panfs/jasmin/era/era-in/netc/ggis/1989/jan1989/ggis198901010000.nc','r')
mask = maskfile.variables['LSM'][:]
maskfile.close()

# best way to mask land values
# make empty array for moisture flux divergence over ocean only:
tcdq_masked = pl.zeros_like(tcdq_filt[0])
for i in range(0,tcdq_mean[0].shape[1]):
    for j in range(0,tcdq_mean[0].shape[2]):
        tcdq_masked[0,i,j] = pl.ma.masked_where(mask[0,0,i,j]>0.5, tcdq_filt2[0,0,i,j])
        # any tcdq value with mask value >0 is NaN, non-ocean values!
tcdq_ocean = pl.ma.masked_invalid(tcdq_masked) #changes NaNs to '-'

#need lat & lon from ONE FILE ONLY:
latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                                year_input[0] + '/' + filenames[0,0],'r')
lon = latlon.variables['longitude'][:]
lat = latlon.variables['latitude'][:]
latlon.close()


scale = (86400*365)/1000
tcdq_zonal = pl.mean(tcdq_filt2*scale,axis=-1)
tcdq_ocz = pl.nanmean(tcdq_ocean*scale,axis=-1)
pl.figure(1)
pl.plot(lat,pl.squeeze(pl.fliplr(tcdq_zonal)),color='k',label=r'Global $\nabla\cdot\overline{\rho q\mathbf{u}}$')
pl.plot(lat,pl.squeeze(tcdq_ocz),color='blue',label=r'Ocean $\nabla\cdot\overline{\rho q\mathbf{u}}$')
pl.legend(loc=0)
pl.xlim(-90,90)
pl.xticks([-90,-60,-30,0,30,60,90])
pl.axhline(y=0,color='k',ls='--')
pl.ylim(-1.5,3)
pl.ylabel('m/yr')
pl.xlabel('Latitude ($^\circ$)')
pl.title('Zonally averaged moisture flux divergence'); pl.close()

# make empty array for moisture flux divergence over ocean only:
"""tcdq_ocean = pl.zeros_like(overall_mean) # same shape as tcdq
for i in range(0,256):
    for j in range(0,512):
        if mask[0,0,i,j] == 0.:
            tcdq_ocean[0,0,i,j] = overall_mean[0,0,i,j]
        else:
            tcdq_ocean[0,0,i,j] = 0. # prefer to have this as float('nan')"""



"""tcdq_months_masked = pl.zeros_like(tcdq_monthly_means[:,0])
tcwv_months_masked = pl.zeros_like(tcwv_monthly_means[:,0])
for m in range(0,tcdq_monthly_means.shape[0]):
    for i in range(0,tcdq_monthly_means.shape[3]):
        for j in range(0,tcdq_monthly_means.shape[4]):
            tcdq_months_masked[m,0,i,j] = pl.ma.masked_where(mask[0,0,i,j]>0., tcdq_monthly_means[m,0,0,i,j])
            tcwv_months_masked[m,0,i,j] = pl.ma.masked_where(mask[0,0,i,j]>0., tcwv_monthly_means[m,0,0,i,j])
tcdq_ocean_months = pl.ma.masked_invalid(tcdq_months_masked)
tcwv_ocean_months = pl.ma.masked_invalid(tcwv_months_masked)"""

#mask land values in E-P monthly means:
#EmP_mm_masked = pl.zeros_like(EmP_monthly_means) #mm stands for 'monthly_means'
#for m in range(0,EmP_monthly_means.shape[0]):
#    for i in range(0,EmP_monthly_means.shape[3]):
#        for j in range(0,EmP_monthly_means.shape[4]):
#            EmP_mm_masked[m,0,0,i,j] = pl.ma.masked_where(mask[0,0,i,j]>0.,EmP_monthly_means[m,0,0,i,j])
#EmP_mm_ocean = pl.ma.masked_invalid(EmP_mm_masked)



#----------------------------CONVERT TO E-P--------------------------------
# divide by density of liquid water, rho = 0.805 kg/m**3
# multiply by number of seconds in a day 60*60*24 = 86400 s
# no need to multiply by 365 days as mean has been calculated across 12 months
EmP_all = (tcdq_ann_mean/(10**3))*86400
EmP_ocean = (tcdq_ocean/(10**3))*86400

#print 'Maximum evaporation = ', EmP.max(), 'm/yr'
#max_index = pl.where(EmP==EmP.max()) #index of max evaporation in EmP array
#max_loc = [max_index[1][0],max_index[2][0]]
#max_latlon = [lat[max_loc[0]],lon[max_loc[1]]] # not sure about these locations
#print 'Maximum precipitation = ', EmP.min(), 'm/yr'
#min_index = pl.where(EmP==EmP.min())
#min_loc = [min_index[1][0],min_index[2][0]]
#min_latlon = [lat[min_loc[0]],lon[min_loc[1]]]

#http://nbviewer.ipython.org/github/Unidata/netcdf4-python/blob/master/examples/writing_netCDF.ipynb
"""newnc = Dataset('test.nc',mode='w',format='NETCDF4')
lat_dim = newnc.createDimension('lat', 256)
lon_dim = newnc.createDimension('lon', 512)
time_dim = newnc.createDimension('time', None)
lat_in = newnc.createVariable('lat', pl.float32, ('lat',))
lat_in.units = 'degrees_north'
lat_in.long_name = 'latitude'
lon_in = newnc.createVariable('lon', pl.float32, ('lon',))
lon_in.units = 'degrees_east'
lon_in.long_name = 'longitude'
time = newnc.createVariable('time', pl.float64, ('time',))
#time.units = '' #what should time units be?
time.long_name = 'time'

TCDQ = newnc.createVariable('tcdq',pl.float64,('time','lat','lon'))
TCDQ.units = 'kg m**-2 s**-1'
TCDQ.standard_name = 'vertical integral of divergence of moisture flux'
#nlats = len(lat_dim); nlons = len(lon_dim); ntimes = 3
#lat[:] = 90. - (180./nlats)*pl.arange(nlats)
#lon[:] = (180./nlats)*pl.arange(nlons)
lat_in[:] = lat # straight from ERA-Interim nc file
lon_in[:] = lon # straight from ERA-Interim nc file
TCDQ[:,:,:] = tcdq_ann_mean

lsm = newnc.createVariable('LSM',pl.float32,('time','lat','lon'))
lsm.units = '0 to 1'
lsm.standard_name = 'Land-Sea Mask'
lsm[:,:,:] = mask
#ncview really doesn't like the mask

TCDQ_ocean = newnc.createVariable('tcdq ocean',pl.float64,('time','lat','lon'))
TCDQ_ocean.units = 'kg m**-2 s**-1'
TCDQ_ocean.standard_name = 'vertical integral of divergence of moisture flux over oceans'
TCDQ_ocean[:,:,:] = tcdq_ocean
#ncview likes a masked array even less

newnc.close()"""

#need to write monthly means to a different netcdf file with 12 time steps
newmask = pl.zeros([12,1,256,512]) # land-sea mask must have 12 time steps
newmask[:] = mask # set every time step of newmask equal to the land-sea mask

"""nc_months = Dataset('monthly_means.nc',mode='w',format='NETCDF4')

lat_dim = nc_months.createDimension('lat', 256)
lon_dim = nc_months.createDimension('lon', 512)
time_dim = nc_months.createDimension('time', 12)

lat_in = nc_months.createVariable('lat', pl.float32, ('lat',))
lat_in.units = 'degrees_north'
lat_in.long_name = 'latitude'
lon_in = nc_months.createVariable('lon', pl.float32, ('lon',))
lon_in.units = 'degrees_east'
lon_in.long_name = 'longitude'

time = nc_months.createVariable('time', pl.float64, ('time',))
time.units = 'months'
time.long_name = 'time'

TCDQ = nc_months.createVariable('tcdq',pl.float64,('time','lat','lon'))
TCDQ.units = 'kg m**-2 s**-1'
TCDQ.standard_name = 'vertical integral of divergence of moisture flux'
#nlats = len(lat_dim); nlons = len(lon_dim); ntimes = 3
lat_in[:] = lat # straight from ERA-Interim nc file
lon_in[:] = lon # straight from ERA-Interim nc file
TCDQ[:,:,:] = tcdq_monthly_means[:,0,0,:,:]

TCWV = nc_months.createVariable('tcwv',pl.float64,('time','lat','lon'))
TCWV.units = 'kg m**-2'
TCWV.standard_name = 'total column water vapour'
TCWV[:,:,:] = tcwv_monthly_means[:,0,0,:,:]

EmP = nc_months.createVariable('EmP',pl.float64,('time','lat','lon'))
EmP.units = 'mm/day'
EmP.standard_name = 'Evaporation minus precipitation'
EmP[:,:,:] = EmP_monthly_means[:,0,0,:,:]

lsm = nc_months.createVariable('LSM',pl.float32,('time','lat','lon'))
lsm.units = '0 to 1'
lsm.standard_name = 'Land-Sea Mask'
lsm[:] = newmask[:,0,:,:]

nc_months.close()"""

#-----------------------PLOT OVERALL MEAN----------------------------------



# can't remember why these are here, need to check
lon_0 = lon.mean()
lat_0 = lat.mean()

pl.figure(2)
m = Basemap(projection='cyl',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
        llcrnrlon=-180,urcrnrlon=180,lat_ts=20)
tcdq_filt, newlon = shiftgrid(180.0, tcdq_filt*(86400/1000), lon, start=False)
m.drawcoastlines()
lons, lats = pl.meshgrid(newlon,lat)
X, Y = m(lons,lats)

norm = MidpointNormalize(midpoint=0)

cmap = pl.get_cmap('seismic')
#cs = m.pcolormesh(X, Y, pl.squeeze(tcdq_filt), cmap=cmap, norm=pl.Normalize(-4,4,clip=False))#,vmin=-4,vmax=4)
cf = m.contourf(X,Y, pl.squeeze(tcdq_filt),cmap=cmap,norm=pl.Normalize(-4,4,clip=False),#vmin=-4,vmax=4,
                levels=[-4,-2,-1,-0.5,-0.25,-0.125,0.125,0.25,0.5,1,2,4],extend="min")
m.fillcontinents(color='gainsboro')
#cmap.set_over(color='maroon')
#cmap.set_under(color='#000047')
#cmap.set_bad('lightgrey')
cbar = m.colorbar(cf,location='bottom', pad = "10%",extend="min",ticks=[-4,-1,-0.25,0,0.25,1,4])
cbar.set_label('m yr$^{-1}$',fontsize=32)
cbar.ax.tick_params(labelsize=32)
#tick_locator = ticker.MaxNLocator(nbins=7)
#cbar.locator = tick_locator
cbar.update_ticks()
pl.title('$e-p$',fontsize=45,position=(0.5,1.01))
m.drawparallels([-60,-30,0,30,60],labels=[1,0,0,1],linewidth=0,ax=None,fontsize=32)
m.drawmeridians([-90,0,90],labels=[1,0,0,1],linewidth=0,ax=None,fontsize=32)
#pl.savefig('tcdq_mean_1979-2014.png')

#-------------------------PLOT MASKED MEAN---------------------------------

# to plot masked tcdq:
#use same lat and lon, same map
# tcdq_ocean, lon = shiftgrid(180.0, tcdq_ocean, lon, start=False)
# same coastlines, meshgrid, norm, cmap
# cs = m.pcolormesh(X, Y, pl.squeeze(tcdq_ocean), cmap=cmap, norm=norm)
#m.fillcontinents(color='white')
# same tick stuff
# title will be same and filename will be different

#-------------------------------PLOT E-P-----------------------------------

#use same lat and lon, same map
# EmP, lon = shiftgrid(180.0, EmP, lon, start=False)
# same coastlines, meshgrid, norm, cmap
# cs = m.pcolormesh(X, Y, pl.squeeze(EmP), cmap=cmap, norm=norm)
#m.fillcontinents(color='white')
# same tick stuff
# title and filename

#latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
#                                year_input[0] + '/' + filenames[0,0],'r')
#lon = latlon.variables['longitude'][:]
#lat = latlon.variables['latitude'][:]
#latlon.close()

#---------------------------PLOT E-P MONTHLY MEANS-------------------------
   
months_names = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
   
fig, ax = pl.subplots(3,4,figsize=(24,11))
for i in range(12):
    latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                                year_input[0] + '/' + filenames[0,0],'r')
    lon = latlon.variables['longitude'][:]
    lat = latlon.variables['latitude'][:]
    latlon.close()
    lon_0 = lon.mean()
    lat_0 = lat.mean()
    pl.subplot(3,4,i+1)#,sharex=fig,sharey=fig) # subplot starts at 1 hence i+1 needed to index subplots
    m = Basemap(projection='cyl',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
                    llcrnrlon=-180,urcrnrlon=180,lat_ts=20)
    EmP_mm[i], newlon = shiftgrid(180.0, EmP_mm[i], lon, start=False)
    m.drawcoastlines(linewidth=0.5)
    lons, lats = pl.meshgrid(newlon,lat)
    X, Y = m(lons, lats)
    #norm = MidpointNormalize(midpoint=0) # calls a class that centres the colourmap at zero
    cmap = pl.get_cmap('seismic')
    #cmap.set_bad('lightgrey')
    # constrain all data within the same bounds to use use same colour map:
    pcm=pl.pcolormesh(X,Y,pl.squeeze(EmP_mm[i]),cmap=cmap,norm=pl.Normalize(-10,10,clip=False))
    m.fillcontinents(color='gainsboro')
    if i in (0,4,8):
        m.drawparallels([-60,-30,0,30,60],labels=[1,0,0,1],linewidth=0,ax=None,fontsize=14)
    if i > 7:
        m.drawmeridians([-90,0,90],labels=[1,0,0,1],linewidth=0,ax=None,fontsize=14)
    pl.title('$E-P$  '+months_names[i],fontsize=16)
    
f = pl.gcf()
colax = f.add_axes([0.9,0.25,0.015,0.5]) # set position of colour bar
#pcmbounds = pl.linspace(-20,20,11)
#cmap.set_over(color='maroon')
#cmap.set_under(color='#000047')
clb = pl.colorbar(pcm, cax=colax,extend='min')#,boundaries=pcmbounds)
clb.set_label('$E-P$  (mm day$^{-1}$)',fontsize=16)
clb.ax.tick_params(labelsize=14) 
pl.subplots_adjust(wspace=0.04,hspace=-0.3,bottom=0.1,top=0.9,right=0.89)
#pl.close()
    
#if __name__=='__main__':
#    main()