# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 13:36:28 2015

@author: np838619
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas

#read in csv file containing E-P-R data
#file has data for the Atlantic, Pacific & Indian oceans
FWflux = pandas.read_csv('FWflux_data.csv',sep=',',header=1,skiprows=1,index_col=0,
    names=['A_bas','A_bas_err','A_4560','A_4560_err','A_3545','A_3545_err','A_nth','A_nth_err',
           'A_mid','A_mid_err','A_sth','A_sth_err',
           'P_bas','P_bas_err','P_47BS','P_47BS_err','P_3047','P_3047_err','P_nth','P_nth_err',
           'P_sth','P_sth_err',
           'I_bas','I_bas_err','I_nth','I_nth_err','I_mid','I_mid_err','I_sth','I_sth_err'])

A = FWflux.A_bas
P = FWflux.P_bas
I = FWflux.I_bas

Ae = FWflux.A_bas_err
Pe = FWflux.P_bas_err
Ie = FWflux.I_bas_err

#data
# Bottom 2 bands from Schanze (2010), 4th columns, removed because of different Southern Ocean boundary
# Actual values from Schanze are 0.1, 0.01 and 0.57
# These values should be between about 35S and 45N
# S10 Atlantic value changed as I was coming from Greenland down
# V14 presents integrals north of a band, had to do some subtracting to get it from 45N down!

#Ar = [0.016,-,0.24,0.3,-] # Arctic P-E data
#So = [-,0.8,0.67,0.64,0.61]

#B = [A[0],P[0],I[0]] # Baumgartner & Reichel (1975)
#W = [A[1],P[1],I[1]] # Wijffels (1992)
D = [A[2],P[2],I[2]] # Dai & Trenberth (2003)
G = [A[3],P[3],I[3]] # Ganachaud & Wunsch (2003)
T = [A[4],P[4],I[4]] # Talley (2008)
S = [A[5],P[5],I[5]] # Schanze et al. (2010)
V = [A[6],P[6],I[6]] # Valdivieso et al. (2014)
E = [A[7],P[7],I[7]] # ERA-Interim moisture flux divergence & Dai & Trenberth (2002) runoff

flux = np.array([D,T,S,V,E]) # stick them all in one array

D_err = [Ae[2],Pe[2],Ie[2]] # nothing actually in here
G_err = [Ae[3],Pe[3],Ie[3]] # calculated with RMS in paper, I rounded to 2 decimal places
T_err = [Ae[4],Pe[4],Ie[4]] # calculated with RMS, hard to tell what Talley actually did though
S_err = [Ae[5],Pe[5],Ie[5]] # nothing actually in here
V_err = [Ae[6],Pe[6],Ie[6]] # explicitly stated in abstract
E_err = [Ae[7],Pe[7],Ie[7]] # standard deviation from inter-annual variability

errs = np.array([D_err,T_err,S_err,V_err,E_err]) # stick them all in one array

ind = np.arange(3)
width = 0.15 # width of each bar

clrs = ['g', 'c', 'y', 'b', 'magenta']
names = ['Dai & Trenberth (2003)','Talley (2008)',
         'Schanze et al. (2010)','Valdivieso et al. (2014)','ERA-I & DT02 runoff']
m = [] # empty list needed for the first part of the bar instance
fig, ax = plt.subplots()
for i in range(flux.shape[0]):
    q = ax.bar(ind+i*width,flux[i],width,color=clrs[i],yerr=errs[i],
               error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
    m.append(q[0])
    
ax.legend(m,names,loc=1,ncol=2)
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-0.2,1.3))

ax.set_ylabel('$E-P-R$ (Sv)',fontsize=16)
ax.set_title('Net surface water flux')
ax.set_xticks(ind+3*width)
ax.set_xticklabels(('Atlantic','Pacific','Indian'))

plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')

plt.savefig('EPR_hist60.png')