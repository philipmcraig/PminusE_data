# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 10:09:38 2015

@author: np838619
"""

import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pylab import *

delta = 0.025
x = y = np.arange(-3.0, 3.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
Z = Z2-Z1  # difference of Gaussians
ax = Axes(plt.gcf(),[0,0,1,1],yticks=[],xticks=[],frame_on=False)
plt.gcf().delaxes(plt.gca())
plt.gcf().add_axes(ax)
im = plt.imshow(Z, cmap=cm.gray)

plt.show()