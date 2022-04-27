from netCDF4 import MFDataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap, shiftgrid
import os

#nc = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/2013/hgaw*.nc','r')

MainFolder=r"/panfs/jasmin/era/era-in/netc/monthly_means/"
#files = []
for (path, dirs, ncfile) in os.walk(MainFolder):
    for dirs in range(1979,2014):
        print dirs
        list_of_files = os.listdir(os.getcwd())
        for ncfile in list_of_files:
            if ncfile.startswith('hgaw'):
                print ncfile
    #for i in os.listdir(path):
     #   if os.path.isfile(os.path.join(path,i)) and 'hgaw' in i:
      #      files.append(i)
    #for ncfile in files:
        #print ncfile
        #files = [i for i in os.listdir(path) if os.path.join(path,i)) and \
        #    'hgaw' in i]
        #if ncfile[4:]=='hgaw':
            ncfile=os.path.join(path,ncfile)
            ncfile=MFDataset('/panfs/jasmin/era/era-in/netc/monthly_means/*/hgaw201301*.nc')
            tcdq = ncfile.variables['TCDQ'][:]
            lon = ncfile.variables['longitude'][:]
            lat = ncfile.variables['latitude'][:]
            ncfile.close()

"""            big_array=[]
            for j in tcdq:
                big_array.append(j)
                big_array=np.array(big_array)
                mean=np.mean(big_array)

            m = Basemap(projection='merc',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
                llcrnrlon=-180,urcrnrlon=180,lat_ts=20)
            tcdq, lon = shiftgrid(180.0, mean, lon, start=False)
            m.drawcoastlines()
            lon, lat = np.meshgrid(lon,lat)
            X, Y = m(lon,lat)

            cmap = plt.get_cmap('seismic')
            cs = m.pcolormesh(X, Y, mean, cmap=cmap)

            cbar = m.colorbar(cs,location='bottom', pad = "10%")
            cbar.set_label('kg m$^{-2}$ s$^{-1}$')
            tick_locator = ticker.MaxNLocator(nbins=5)
            cbar.locator = tick_locator
            cbar.update_ticks()
            plt.title('Annual Mean vertical integral of divergence of moisture flux')
            plt.show()"""
