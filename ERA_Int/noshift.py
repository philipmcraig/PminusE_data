# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:20:07 2015

@author: np838619
"""

t = []
for line in open("oceanconfig.py","r"):
    lend = line.split('\n')
    if "shiftgrid" not in line:
        t.append(lend)

h=t[9:]
for i in range(len(h)):
    exec(h[i][0])