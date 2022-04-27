# -*- coding: utf-8 -*-
"""
Created on Thu May 28 13:35:32 2015

@author: np838619

For reading in the Dai & Trenberth runoff data.
NEVER use the csv module. It is useless. Pandas does exactly what I want &
is similar to R.
Will combine this with ERA-Interim calculations from a netCDF file as some point.
"""

import pandas
import pylab as pl


runoff = pandas.read_csv('DT03_runoff.csv',sep=',',header=1,skiprows=0,
    names=['row','lat','Atl','Pac','Ind','Glo'])