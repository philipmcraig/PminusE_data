# -*- coding: utf-8 -*-
"""
last updated 5/12/16 2:23PM (6th commit, cluster doesn't have git)
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas

def NanAsX(fluxes,ax,index,width):
    """
    """
    
    nan_array = np.isnan(fluxes)
    
    for j in range(len(nan_array[index])):
        if nan_array[index,j] == True:
            ax.plot(index*width+width/2+j,0,marker='x',color='k',mew=1,ms=7.5)


#read in csv file containing E-P-R data
#file has data for the Atlantic, Pacific & Indian oceans
FWflux = pandas.read_csv('FWflux_data.csv',sep=',',header=1,skiprows=1,index_col=0,
    names=['A_bas','A_bas_err','A_4560','A_4560_err','A_3545','A_3545_err','A_nth','A_nth_err',
           'A_mid','A_mid_err','A_sth','A_sth_err',
           'P_bas','P_bas_err','P_47BS','P_47BS_err','P_3047','P_3047_err','P_nth','P_nth_err',
           'P_mid','P_mid_err','P_sth','P_sth_err',
           'I_bas','I_bas_err','I_nth','I_nth_err','I_mid','I_mid_err','I_sth','I_sth_err'])


A60 = FWflux.A_4560
An = FWflux.A_nth
Am = FWflux.A_mid
As = FWflux.A_sth
Pbs = FWflux.P_47BS
Pn = FWflux.P_nth
Pm = FWflux.P_mid
Ps = FWflux.P_sth
In = FWflux.I_nth
Im = FWflux.I_mid
Is = FWflux.I_sth

A60_err = FWflux.A_4560_err
An_err = FWflux.A_nth_err
Am_err = FWflux.A_mid_err
As_err = FWflux.A_sth_err
Pbs_err = FWflux.P_47BS_err
Pn_err = FWflux.P_nth_err
Pm_err = FWflux.P_mid_err
Ps_err = FWflux.P_sth_err
In_err = FWflux.I_nth_err
Im_err = FWflux.I_mid_err
Is_err = FWflux.I_sth_err

# Atlantic Ocean E-P

D_A = [A60[2],An[2],Am[2],As[2]] # Dai & Trenberth (2003)
G_A = [A60[3],An[3],Am[3],As[3]] # Ganachaud & Wunsch (2003)
T_A = [A60[4],An[4],Am[4],As[4]] # Talley (2008)
S_A = [A60[5],An[5],Am[5],As[5]] # Schanze et al. (2010)
V_A = [A60[6],An[6],Am[6],As[6]] # Valdivieso et al. (2014)
E_A = [A60[7],An[7],Am[7],As[7]] # ERA-Interim moisture flux divergence with Dai & Trenberth (2002) runoff
M_A = [A60[8],An[8],Am[8],As[8]] # ECCO

Aflux = np.array([E_A,D_A,S_A,G_A,T_A,V_A,M_A])

D_A_err = [A60_err[2],An_err[2],Am_err[2],As_err[2]]
G_A_err = [A60_err[3],An_err[3],Am_err[3],As_err[3]]
T_A_err = [A60_err[4],An_err[4],Am_err[4],As_err[4]]
S_A_err = [A60_err[5],An_err[5],Am_err[5],As_err[5]]
V_A_err = [A60_err[6],An_err[6],Am_err[6],As_err[6]]
E_A_err = [A60_err[7],An_err[7],Am_err[7],Am_err[7]]
M_A_err = [A60_err[8],An_err[8],Am_err[8],Am_err[8]]

Aflux_err = np.array([E_A_err,D_A_err,S_A_err,G_A_err,T_A_err,V_A_err,M_A_err])

plt.subplots(3,1,figsize=(10,12))

ind = np.arange(4) # remember that arange is spelt wrong!
width = 0.13 # width of each bar

clrs = ['magenta','g','y','r','c','b','orange']
names = ['ERAI & DT02 runoff','Dai & Trenberth (2003)','Schanze et al. (2010)',
         'Ganachaud & Wunsch (2003)','Talley (2008)','Valdivieso et al. (2014)',
        'ECCOv4']
m = []


ax1 = plt.subplot(311)
for i in range(Aflux.shape[0]):
    q = ax1.bar(ind+i*width,Aflux[i],width,color=clrs[i],yerr=Aflux_err[i],
               error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
    m.append(q[0])
    NanAsX(Aflux,ax1,i,width)

plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-0.7,0.7))

ax1.set_ylabel('$E-P-R$ (Sv)',fontsize=20)
plt.grid(axis='y')
plt.axvline(x=0.95,color='k',ls=':')
plt.axvline(x=1.95,color='k',ls=':')
plt.axvline(x=2.95,color='k',ls=':',)
plt.setp(ax1.get_yticklabels(),fontsize=17)
ax1.set_title('Atlantic',fontsize=20)
ax1.set_xticks(ind+(3)*width)
ax1.set_xticklabels(('45-60$^\circ$N','24-45$^\circ$N',
                            '16$^\circ$S-24$^\circ$N','35-16$^\circ$S'),fontsize=20)
ax1.annotate('(a)',(0,1.02),xycoords='axes fraction',size=25)
plt.yticks(np.linspace(-0.6,0.6,5))
ax1.legend(m,names,loc=4,fontsize=12,ncol=2,columnspacing=1)
#plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
#plt.gcf().subplots_adjust(bottom=0.15)
#plt.tight_layout()
#plt.savefig('atlantic_bands.png')

####################################

# Pacific Ocean E-P

D_P = [Pbs[2],Pn[2],Pm[2],Ps[2]]
G_P = [Pbs[3],Pn[3],Pm[3],Ps[3]]
T_P = [Pbs[4],Pn[4],Pm[4],Ps[4]]
S_P = [Pbs[5],Pn[5],Pm[5],Ps[5]]
V_P = [Pbs[6],Pn[6],Pm[6],Ps[6]]
E_P = [Pbs[7],Pn[7],Pm[7],Ps[7]]
M_P = [Pbs[8],Pn[8],Pm[8],Ps[8]]

Pflux = np.array([E_P,D_P,S_P,G_P,T_P,V_P,M_P])

D_P_err = [Pbs_err[2],Pn_err[2],Pm_err[2],Ps_err[2]]
G_P_err = [Pbs_err[3],Pn_err[3],Pm_err[3],Ps_err[3]]
T_P_err = [Pbs_err[4],Pn_err[4],Pm_err[4],Ps_err[4]]
S_P_err = [Pbs_err[5],Pn_err[5],Pm_err[5],Ps_err[5]]
V_P_err = [Pbs_err[6],Pn_err[6],Pm_err[6],Ps_err[6]]
E_P_err = [Pbs_err[7],Pn_err[7],Pm_err[7],Ps_err[7]]
M_P_err = [Pbs_err[8],Pn_err[8],Pm_err[8],Ps_err[8]]

Pflux_err = np.array([E_P_err,D_P_err,S_P_err,G_P_err,T_P_err,V_P_err,M_P_err])

ind = np.arange(4)
m = []

ax2 = plt.subplot(312)
for i in range(Pflux.shape[0]):
    q = ax2.bar(ind+i*width,Pflux[i],width,color=clrs[i],yerr=Pflux_err[i],
               error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
    m.append(q[0])
    NanAsX(Pflux,ax2,i,width)
    
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-0.7,0.7))

ax2.set_ylabel('$E-P-R$ (Sv)',fontsize=20)
plt.grid(axis='y')
plt.axvline(x=0.95,color='k',ls=':')
plt.axvline(x=1.95,color='k',ls=':')
plt.axvline(x=2.95,color='k',ls=':')
plt.setp(ax2.get_yticklabels(),fontsize=17)
ax2.set_title('Pacific',fontsize=20)
ax2.set_xticks(ind+(3)*width)
ax2.set_xticklabels(('47$^\circ$N-BS','24-47$^\circ$N','17$^\circ$S-24$^\circ$N',
                                     '30$^\circ$S-17$^\circ$S'),fontsize=20)
ax2.annotate('(b)',(0,1.02),xycoords='axes fraction',size=25)
plt.yticks(np.linspace(-0.6,0.6,5))
#ax2.legend(m,names,loc=2,fontsize=12,columnspacing=1,ncol=2)
#plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
#plt.tight_layout()
#plt.savefig('pacific_bands.png')

###########################

# Indian Ocean E-P

D_I = [In[2],Im[2],Is[2]]
G_I = [In[3],Im[3],Is[3]]
T_I = [In[4],Im[4],Is[4]]
S_I = [In[5],Im[5],Is[5]]
V_I = [In[6],Im[6],Is[6]]
E_I = [In[7],Im[7],Is[7]]
M_I = [In[8],Im[8],Is[8]]

Iflux = np.array([E_I,D_I,S_I,G_I,T_I,V_I,M_I])

D_I_err = [In_err[2],Im_err[2],Is_err[2]]
G_I_err = [In_err[3],Im_err[3],Is_err[3]]
T_I_err = [In_err[4],Im_err[4],Is_err[4]]
S_I_err = [In_err[5],Im_err[5],Is_err[5]]
V_I_err = [In_err[6],Im_err[6],Is_err[6]]
E_I_err = [In_err[7],Im_err[7],Im_err[7]]
M_I_err = [In_err[8],Im_err[8],Im_err[8]]

Iflux_err=np.array([E_I_err,D_I_err,S_I_err,G_I_err,T_I_err,V_I_err,M_I_err])

ind = np.arange(3)
m = []

ax3 = plt.subplot(313)
for i in range(Iflux.shape[0]):
    q = ax3.bar(ind+i*width,Iflux[i],width,color=clrs[i],yerr=Iflux_err[i],
               error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
    m.append(q)
    NanAsX(Iflux,ax3,i,width)
    
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-0.7,0.7))

ax3.set_ylabel('$E-P-R$ (Sv)',fontsize=20)
plt.grid(axis='y')
plt.axvline(x=0.95,color='k',ls=':')
plt.axvline(x=1.95,color='k',ls=':')
plt.setp(ax3.get_yticklabels(),fontsize=17)
ax3.set_title('Indian',fontsize=20)
ax3.set_xticks(ind+(3)*width)
ax3.set_xticklabels(('>8$^\circ$S','20-8$^\circ$S','35-20$^\circ$S'),fontsize=20)
ax3.annotate('(c)',(0,1.02),xycoords='axes fraction',size=25)
plt.yticks(np.linspace(-0.6,0.6,5))
#ax3.legend(m,names,loc=4,ncol=2,fontsize=12,columnspacing=1)
#plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
#plt.tight_layout()
#plt.savefig('indian_bands.png')

############################

plt.subplots_adjust(bottom=0.06,top=0.95,hspace=0.24)
#plt.savefig('/home/np838619/PminusE_data/latbands.png')
plt.show()
