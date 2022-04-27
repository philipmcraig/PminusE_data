#!usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas

#read in csv file containing E-P-R data
#file has data for the Atlantic, Pacific & Indian oceans
FWflux = pandas.read_csv('FWflux_data.csv',sep=',',header=1,skiprows=1,index_col=0,
    names=['A_bas','A_bas_err','A_nth','A_nth_err','A_mid','A_mid_err','A_sth','A_sth_err',
           'P_bas','P_bas_err','P_nth','P_nth_err','P_sth','P_sth_err',
           'I_bas','I_bas_err','I_nth','I_nth_err','I_mid','I_mid_err','I_sth','I_sth_err'])

An = FWflux.A_nth
Am = FWflux.A_mid
As = FWflux.A_sth
Pn = FWflux.P_nth
Ps = FWflux.P_sth
In = FWflux.I_nth
Im = FWflux.I_mid
Is = FWflux.I_sth

An_err = FWflux.A_nth_err
Am_err = FWflux.A_mid_err
As_err = FWflux.A_sth_err
Pn_err = FWflux.P_nth_err
Ps_err = FWflux.P_sth_err
In_err = FWflux.I_nth_err
Im_err = FWflux.I_mid_err
Is_err = FWflux.I_sth_err

# Atlantic Ocean E-P

D_A = [An[2],Am[2],As[2]] # Dai & Trenberth (2003)
G_A = [An[3],Am[3],As[3]] # Ganachaud & Wunsch (2003)
T_A = [An[4],Am[4],As[4]] # Talley (2008)
S_A = [An[5],Am[5],As[5]] # Schanze et al. (2010)
V_A = [An[6],Am[6],As[6]] # Valdivieso et al. (2014)

G_A_err = [An_err[3],Am_err[3],As_err[3]]
T_A_err = [An_err[4],Am_err[4],As_err[4]]
#V_A_err = [0.1,0.11,0.16]

ind = np.arange(3) # remember that arange is spelt wrong!
width = 0.15 # width of each bar

fig, ax = plt.subplots()
Dai = ax.bar(ind,D_A,width,color='g')
Gan = ax.bar(ind+width,G_A,width,color='r',yerr=G_A_err,error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
Tal = ax.bar(ind+2*width,T_A,width,color='c',yerr=T_A_err,error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
Sch = ax.bar(ind+3*width,S_A,width,color='y')
Val = ax.bar(ind+4*width,V_A,width,color='b')#,yerr=V_A_err,error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
# i'm a bit unsure about the V14 error bars at the moment
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-0.7,0.7))

ax.set_ylabel('$E-P-R$ (Sv)')
ax.set_title('Surface water flux:  Atlantic Ocean')
ax.set_xticks(ind+(2)*width)
ax.set_xticklabels(('24-45$^\circ$N \n 26.5-47$^\circ$N (V14)',
    '16$^\circ$S-24$^\circ$N \n 19$^\circ$S-24$^\circ$N (GW03) \n 16$^\circ$S-26.5$^\circ$N V(14)',
                                    '32-16$^\circ$S \n 30-19$^\circ$S (GW03)'))

ax.legend((Dai[0],Gan[0],Tal[0],Sch[0],Val[0]),('Dai & Trenberth (2003)','Ganachaud & Wunsch (2003)',
           'Talley (2008)','Schanze et al. (2010)','Valdivieso et al. (2014)'),loc=4)
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('atlantic_bands.png')

####################################

# Pacific Ocean E-P

D_P = [Pn[2],Ps[2]]
G_P = [Pn[3],Ps[3]]
T_P = [Pn[4],Ps[4]]
S_P = [Pn[5],Ps[5]]
V_P = [Pn[6],Ps[6]]

G_P_err = [Pn_err[3],Ps_err[3]]
T_P_err = [Pn_err[4],Ps_err[4]]
#V_P_err = [Pn_err[6],Ps_err[6]]

ind = np.arange(2)

fig, ax = plt.subplots()
Dai = ax.bar(ind,D_P,width,color='g')
Gan = ax.bar(ind+width,G_P,width,color='r',yerr=G_P_err,error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
Tal = ax.bar(ind+2*width,T_P,width,color='c',yerr=T_P_err,error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
Sch = ax.bar(ind+3*width,S_P,width,color='y')
Val = ax.bar(ind+4*width,V_P,width,color='b')
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-0.7,0.7))

ax.set_ylabel('$E-P-R$ (Sv)')
ax.set_title('Surface water flux: Pacific Ocean')
ax.set_xticks(ind+(2)*width)
ax.set_xticklabels(('24-47$^\circ$N','30$^\circ$S-24$^\circ$N \n 32$^\circ$S-24$^\circ$N (V14)'))

ax.legend((Dai[0],Gan[0],Tal[0],Sch[0],Val[0]),('Dai & Trenberth (2003)','Ganachaud & Wunsch (2003)',
           'Talley (2008)','Schanze et al. (2010)','Valdivieso et al. (2014)'),loc=4)
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.savefig('pacific_bands.png')

###########################

# Indian Ocean E-P

D_I = [In[2],Im[2],Is[2]]
G_I = [In[3],Im[3],Is[3]]
T_I = [In[4],Im[4],Is[4]]
S_I = [In[5],Im[5],Is[5]]
V_I = [In[6],Im[6],Is[6]]

G_I_err = [In_err[3],Im_err[3],Is_err[3]]
T_I_err = [In_err[4],Im_err[4],Is_err[4]]
#V_I_err = [In_err[6],Im_err[6],Is_err[6]]

ind = np.arange(3)
fig, ax = plt.subplots()
Dai = ax.bar(ind,D_I,width,color='g')
Gan = ax.bar(ind+width,G_I,width,color='r',yerr=G_I_err,error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
Tal = ax.bar(ind+2*width,T_I,width,color='c',yerr=T_I_err,error_kw=dict(elinewidth=2,ecolor='k',capsize=5))
Sch = ax.bar(ind+3*width,S_I,width,color='y')
Val = ax.bar(ind+4*width,V_I,width,color='b')
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-0.7,0.7))

ax.set_ylabel('$E-P-R$ (Sv)')
ax.set_title('Surface water flux: Indian Ocean')
ax.set_xticks(ind+(2)*width)
ax.set_xticklabels(('>8$^\circ$S \n >5$^\circ$S (S10)','20-8$^\circ$S \n 25-5$^\circ$S (S10)',
            '32-20$^\circ$S \n 35-25$^\circ$S (S10)'))

ax.legend((Dai[0],Gan[0],Tal[0],Sch[0],Val[0]),('Dai & Trenberth (2003)','Ganachaud & Wunsch (2003)',
          'Talley (2008)','Schanze et al. (2010)','Valdivieso et al. (2014)'),loc=4)
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
plt.savefig('indian_bands.png')

############################

plt.show()
