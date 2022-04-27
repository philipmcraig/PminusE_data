import numpy as np
import matplotlib.pyplot as plt

# P-E for each ocean by band scaled by band surface area, (P-E)/area

# Atlantic Ocean
# 24-45N total surface area is 1.82x10^13 m^2 from Talley (2008)
# 16S-24N total surface area is 3.01x10^13 m^2
# 32-16S total surface area is 0.92x10^13 m^2

a1 = 1.82e10 # scaled to km^2
a2 = 3.01e10
a3 = 0.92e10

N=3
D_A = [-0.22/a1,-0.13/a2,-0.22/a3] # Dai & Trenberth (2002)
G_A = [0/a1,-0.1/a2,-0.36/a3] # Ganachaud & Wunsch (2003)
T_A = [-0.03/a1,-0.1/a2,-0.46/a3] # Talley (2008)
S_A = [-0.43/a1,0.17/a2,-0.23/a3] # Schanze et al. (2010)

ind = np.arange(N) # remember that arange is spelt wrong!
width = 0.15 # width of each bar

fig, ax = plt.subplots()
Dai = ax.bar(ind,D_A,width,color='g')
Gan = ax.bar(ind+width,G_A,width,color='r')
Tal = ax.bar(ind+2*width,T_A,width,color='c')
Sch = ax.bar(ind+3*width,S_A,width,color='y')
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-1e-10,1e-10))

ax.set_ylabel('P-E/area (Sv/km$^{-2}$)') # Sverdrups per square km
ax.set_title('Precipitation minus Evaporation over the Atlantic Ocean \n')
ax.set_xticks(ind+(2)*width)
ax.set_xticklabels(('24-45$^\circ$N','16$^\circ$S-24$^\circ$N \n 19$^\circ$S-24$^\circ$N (GW03)','32-16$^\circ$S \n 30-19$^\circ$S (GW03)'))

ax.legend((Dai[0],Gan[0],Tal[0],Sch[0]),('GW03','DT02','T08','S10'))
plt.savefig('atlantic_scaled.png')

############################################

# Pacific Ocean
# 24-47N total surface area is 2.51x10^13 m^2
# 30S-24N total surface area is 9.23x10^13 m^2

p1 = 2.51e10
p2 = 9.23e10

M=2
D_P = [-0.06/p1,-0.1/p2]
G_P = [0.14/p1,-0.2/p2]
T_P = [0.08/p1,-0.14/p2]
S_P = [0.03/p1,-0.49/p2]

ind = np.arange(M)

fig, ax = plt.subplots()
Dai = ax.bar(ind,D_P,width,color='g')
Gan = ax.bar(ind+width,G_P,width,color='r')
Tal = ax.bar(ind+2*width,T_P,width,color='c')
Sch = ax.bar(ind+3*width,S_P,width,color='y')
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-10e-12,10e-12))

ax.set_ylabel('P-E/area (Sv/km$^{-2}$)')
ax.set_title('Precipitation minus Evaporation over the Pacific Ocean \n')
ax.set_xticks(ind+(2)*width)
ax.set_xticklabels(('24-47$^\circ$N','30$^\circ$S-24$^\circ$N'))

ax.legend((Dai[0],Gan[0],Tal[0],Sch[0]),('DT02','GW03','T08','S10'))
plt.savefig('pacific_scaled.png')

###############################################

# Indian Ocean
# > 8S total surface area is 2.55x10^13 m^2
# 20-8S total surface area is 0.92x10^13 m^2
# 32-20S total surface area is 0.98x10^13 m^2

i1 = 2.55e10
i2 = 0.92e10
i3 = 0.98e10

L=3
D_I = [0.13/i1,-0.08/i2,-0.38/i3]
G_I = [0.1/i1,-0.33/i2,-0.35/i3]
T_I = [0.36/i1,-0.43/i2,-0.31/i3]
S_I = [0.02/i1,-0.33/i2,-0.28/i3]

ind = np.arange(L)
fig, ax = plt.subplots()
Gan = ax.bar(ind,G_I,width,color='r')
Dai = ax.bar(ind+width,D_I,width,color='g')
Tal = ax.bar(ind+2*width,T_I,width,color='c')
Sch = ax.bar(ind+3*width,S_I,width,color='y')
plt.axhline(linewidth=1, color='black') # adds a line at y=0 because I wanted to show x-axis
plt.ylim((-1e-10,1e-10))

ax.set_ylabel('P-E/area (Sv/km$^{-2}$)')
ax.set_title('Precipitation minus Evaporation over the Indian Ocean \n')
ax.set_xticks(ind+(2)*width)
ax.set_xticklabels(('>8$^\circ$S','20-8$^\circ$S','32-20$^\circ$S'))

ax.legend((Dai[0],Gan[0],Tal[0],Sch[0]),('GW03','DT02','T08','S10'))
plt.savefig('indian_scaled.png')

###############################################

plt.show()
