import os
from netCDF4 import MFDataset

for i in os.listdir("/panfs/jasmin/era/era-in/netc/monthly_means/2013/"):
    if i.startswith("hgaw"):
        print i
        #F = MFDataset('i')
        #F.close()
        continue
