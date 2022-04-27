#from netCDF4 import MFDataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap, shiftgrid
import os

MainFolder="/panfs/jasmin/era/era-in/netc/monthly_means/"

time = range(0,24,6)

for Y in range(1979,2014):
    for M in range(1,12):
        for T in time(0:3)
            yr = '%d', Y
            mn = '%2.2d', M
            tm = '%4.0d', T
        
        #READ file here?

        full = MainFolder = '/' + yr + 'hgaw' + mn + tm + '.nc'
