# calculate annual mean vertical integral of divergence of moisture flux

import sys
from netCDF4 import Dataset, MFDataset

# find /panfs/jasmin/era/era-in/netc/monthly_means/2013/hgaw*nc | grep -v atest* | python mean_test.py

X = sys.stdin.readlines()
print X
print len(X)

#Y = MFDataset("/panfs/jasmin/era/era-in/netc/monthly_means/1979/hgaw*.nc")

for i in range(0,len(X)):
    ncfile[i] = Dataset(,'r')
    tcdq[i] = ncfile[i].variables['TCDQ'][:]
    ncfile.close()
s
