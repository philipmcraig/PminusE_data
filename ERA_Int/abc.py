from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap, shiftgrid
import os

MainFolder=r"/panfs/jasmin/era/era-in/netc/monthly_means/"

files=[]
for (path, dirs, ncfile) in os.walk(MainFolder):
    for dirs in xrange(1979,2015):
        print dirs
        list_of_files = os.listdir(os.getcwd())
        for ncfile in list_of_files:
            if ncfile.startswith('hgaw'):
                print ncfile
