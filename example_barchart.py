import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()

ax = fig.add_subplot(1,1,1)

x = [1,2,3,4]

#data
A = -0.543 # Atlantic loses 0.543e6 freshwater (Wijffels et al.,1992)
P = 0.887 # Pacific gains 0.887 freshwater
I = -0.441 # Indian loses 0.441e6 freshwater
Ar = 0.095 # Arctic gains 0.096e6 freshwater
PminusE = [A,P,I,Ar] # make array of P-E data

ax.bar(x,PminusE,width=0.3)

plt.show()
