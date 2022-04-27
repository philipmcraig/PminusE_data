from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas

plt.close('all')
#read in csv file containing E-P-R data
#file has data for the Atlantic, Pacific & Indian oceans
FWflux = pandas.read_csv('FWflux_data.csv',sep=',',header=1,skiprows=1,index_col=0,
    names=['A_bas','A_bas_err','A_4560','A_4560_err','A_3545','A_3545_err','A_nth','A_nth_err',
           'A_mid','A_mid_err','A_sth','A_sth_err',
           'P_bas','P_bas_err','P_47BS','P_47BS_err','P_3047','P_3047_err','P_nth','P_nth_err',
           'P_mid','P_mid_err','P_sth','P_sth_err',
           'I_bas','I_bas_err','I_nth','I_nth_err','I_mid','I_mid_err','I_sth','I_sth_err'])

A = FWflux.A_3545
P = FWflux.P_3047
I = FWflux.I_bas

Ae = FWflux.A_3545_err
Pe = FWflux.P_3047_err
Ie = FWflux.I_bas_err

#data
# Bottom 2 bands from Schanze (2010), 4th columns, removed because of different Southern Ocean boundary
# Actual values from Schanze are 0.1, 0.01 and 0.57
# These values should be between about 35S and 45N
# S10 Atlantic value changed as I was coming from Greenland down
# V14 presents integrals north of a band, had to do some subtracting to get it from 45N down!

#Ar = [0.016,-,0.24,0.3,-] # Arctic P-E data
#So = [-,0.8,0.67,0.64,0.61]

D = [A[2],P[2],I[2]] # Dai & Trenberth (2003)
G = [A[3],P[3],I[3]] # Ganachaud & Wunsch (2003)
T = [A[4],P[4],I[4]] # Talley (2008)
S = [A[5],P[5],I[5]] # Schanze et al. (2010)
V = [A[6],P[6],I[6]] # Valdivieso et al. (2014)
E = [A[7],P[7],I[7]] # ERA-Interim moisture flux divergence & Dai & Trenberth (2002) runoff
M = [A[8],P[8],I[8]] # MIT ECCO Model

flux45 = np.array([E,D,S,G,T,V,M])

D_err = [Ae[2],Pe[2],Ie[2]]
G_err = [Ae[3],Pe[3],Ie[3]] # calculated with RMS in paper, I rounded to 2 decimal places
T_err = [Ae[4],Pe[4],Ie[4]] #calculated with RMS, hard to tell what Talley actually did though
S_err = [Ae[5],Pe[5],Ie[5]]
V_err = [Ae[6],Pe[6],Ie[6]] #explicitly stated in abstract
E_err = [Ae[7],Pe[7],Ie[7]] # standard deviation from inter-annual variability
M_err = [Ae[8],Pe[8],Ie[8]]

errs45 = np.array([E_err,D_err,S_err,G_err,T_err,V_err,M_err])

ind = np.arange(3)
width = 0.13 # width of each bar

clrs45 = ['magenta','g','y','r','c', 'b', 'orange']
names45 = ['ERA-I & DT02 runoff','Dai & Trenberth (2003)','Schanze et al. (2010)',
           'Ganachaud & Wunsch (2003)','Talley (2008)','Valdivieso et al. (2014)',
         'ECCOv4']
m1 = []
fig, ax = plt.subplots(2,1,figsize=(10,10))
ax1 = plt.subplot(211)

for i in range(flux45.shape[0]):
    q = ax1.bar(ind+i*width,flux45[i],width,color=clrs45[i],yerr=errs45[i],
               error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
    m1.append(q[0])

#ax1.legend(m,names45,loc=0,ncol=2)
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-0.5,1.5))

ax1.set_ylabel('$E-P-R$ (Sv)',fontsize=20)
plt.setp(ax1.get_yticklabels(),fontsize=17)
ax1.set_xticks(ind+(7/2)*width)
ax1.set_xticklabels(('Atlantic','Pacific','Indian'),fontsize=20)
ax1.annotate('(a)',(0,1.02),xycoords='axes fraction',size=25)
#ax1.legend(m1,names45,loc=2,ncol=2,fontsize=13.5,columnspacing=1)
plt.grid(axis='y')
plt.axvline(x=0.96,color='k',ls=':')
plt.axvline(x=1.96,color='k',ls=':')
#plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.title('35$^\circ$S-45$^\circ$N',fontsize=20)

#-------------------------------------------------------------------------------

A = FWflux.A_bas
P = FWflux.P_bas
# Indian ocean same as above

Ae = FWflux.A_bas_err
Pe = FWflux.P_bas_err
# Indian ocean same as above

D = [A[2],P[2],I[2]] # Dai & Trenberth (2003)
T = [A[4],P[4],I[4]] # Talley (2008)
S = [A[5],P[5],I[5]] # Schanze et al. (2010)
V = [A[6],P[6],I[6]] # Valdivieso et al. (2014)
E = [A[7],P[7],I[7]] # ERA-Interim moisture flux divergence & Dai & Trenberth (2002) runoff
M = [A[8],P[8],I[8]] # ECCO

flux60 = np.array([E,D,S,T,V,M])

D_err = [Ae[2],Pe[2],Ie[2]]
T_err = [Ae[4],Pe[4],Ie[4]] #calculated with RMS, hard to tell what Talley actually did though
S_err = [Ae[5],Pe[5],Ie[5]]
V_err = [Ae[6],Pe[6],Ie[6]] #explicitly stated in abstract
E_err = [Ae[7],Pe[7],Ie[7]] # standard deviation from inter-annual variability
M_err = [Ae[8],Pe[8],Ie[8]]

errs60 = np.array([E_err,D_err,S_err,T_err,V_err,M_err])

clrs60 = ['magenta','g','y','c',  'b', 'orange']
names60 = ['ERA-I & DT02 runoff','Dai & Trenberth (2003)','Schanze et al. (2010)',
           'Talley (2008)','Valdivieso et al. (2014)','ECCOv4']
m2 = []

#fig,ax2 = plt.subplots(figsize=(9,5.5))
ax2 = plt.subplot(212)

for i in range(flux60.shape[0]):
    q = ax2.bar(ind+i*width,flux60[i],width,color=clrs60[i],yerr=errs60[i],
               error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
    m2.append(q[0])

plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-0.5,1.5))
ax2.set_ylabel('$E-P-R$ (Sv)',fontsize=20)
ax2.annotate('(b)',(0,1.02),xycoords='axes fraction',size=25)
plt.grid(axis='y')
plt.axvline(x=0.95,color='k',ls=':')
plt.axvline(x=1.95,color='k',ls=':')
#plt.axhline(y=0.47,color='k',ls='--',lw=2)
plt.setp(ax2.get_yticklabels(),fontsize=17)
ax2.set_xticks(ind+3*width)
ax2.set_xticklabels(('Atlantic\n35$^\circ$S-60$^\circ$N',
                    'Pacific\n30$^\circ$S-BS',
                    'Indian\n>35$^\circ$S'),fontsize=20)
#ax.annotate('(b)',(0,1.02),xycoords='axes fraction',size=25)
ax2.legend(m1,names45,loc=1,ncol=2,fontsize=13.5,columnspacing=1)
plt.title('35$^\circ$S-65$^\circ$N',fontsize=20)

plt.subplots_adjust(top=0.95,bottom=0.14)
plt.tight_layout()
#plt.savefig('/home/np838619/PminusE_data/EPR_barcharts.png',dpi=300)
plt.show()
