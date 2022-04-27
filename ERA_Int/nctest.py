from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.colors import Normalize

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        a, b = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, a, b))

#vertical integral of divergence of moisture flux
#monthly mean for 1200 December 2013
nc = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/2013/hgaw2013121200.nc','r')

tcuq = nc.variables['TCUQ'][:]
tcvq = nc.variables['TCVQ'][:]
tcdq = nc.variables['TCDQ'][:]
lon = nc.variables['longitude'][:]
lat = nc.variables['latitude'][:]

lon_0 = lon.mean()
lat_0 = lat.mean()



m = Basemap(projection='merc',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20)
tcdq, lon = shiftgrid(180.0, tcdq, lon, start=False)
m.drawcoastlines()
lon, lat = np.meshgrid(lon,lat)
X, Y = m(lon,lat)

norm = MidpointNormalize(midpoint=0)

cmap = plt.get_cmap('seismic')
cs = m.pcolormesh(X, Y, np.squeeze(tcdq), cmap=cmap, norm=norm)


cbar = m.colorbar(cs,location='bottom', pad = "10%")
cbar.set_label('kg m$^{-2}$ s$^{-1}$')
tick_locator = ticker.MaxNLocator(nbins=3)
cbar.locator = tick_locator
cbar.update_ticks()
plt.title('Mean vertical integral of divergence of moisture flux 12Z December 2013')
#plt.savefig('tcqd_2013121200.png')
plt.show()
