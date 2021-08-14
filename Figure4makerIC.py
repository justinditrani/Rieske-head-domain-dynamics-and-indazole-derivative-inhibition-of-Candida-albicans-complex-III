#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 12:16:55 2020

@author: justin
"""

#This program find changes in the prteitn motion described by compponent 0 in J1036 upon binding of ligand csIC_Z (Inz5 bound), csIC_A(inhibitor free).
#please email justin.ditrani@gmail.com if there are any questions. 

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import pandas
import operator
import collections
import os as os
import statistics as stats

def open_cs(filename):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file = dir_path + '/' + filename
    csfile = np.load(file)
    return csfile if type(csfile) is np.ndarray else np.load(csfile)

def H_A(csarray):
    j = []
    k = []
    for i in range(len(csarray['blob/path'])):
        if csarray['blob/path'][i][1:2] == b'1':
            j = np.append(j,i)
        if csarray['blob/path'][i][1:2] == b'2':
            k = np.append(k,i)
    qq = j[0].astype(int)
    print(qq)
    return csarray[0:qq-1],csarray[qq:]

def MK_cs(cs_PCA,csA_2,UU):
    YY = np.where(np.isin(cs_PCA['uid'][:],csA_2['uid'][UU]).astype(int) == 1)
    return cs_PCA[YY]
    

csALL = open_cs('cryosparc_P26_J877_particles.cs')

csIC = open_cs('cryosparc_P26_J1036_particles.cs')

csIC_Z,csIC_A = H_A(csIC)

##both
fig = plt.figure(figsize=(28,11))
ax = fig.add_subplot(111)


plt.hist(csIC_Z['components_mode_0/value'],bins = 100, histtype = 'stepfilled',linewidth=8,ec = 'k',facecolor = 'steelblue',alpha = 0.7,density=True,label='lnz-5 bound')
plt.hist(csIC_A['components_mode_0/value'],bins = 100, histtype = 'stepfilled',linewidth=8,ec = 'k',facecolor = 'powderblue',alpha = 0.7,density=True,label='APO')

plt.locator_params(axis='y', nbins=2)
plt.locator_params(axis='x', nbins=2)
fs =100
ax.set_ylabel('Population density',size =fs)
ax.set_xlabel('Variability component',size=fs)

ax.tick_params(labelsize=fs)
ax.legend(fontsize =70,loc='upper right')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(-200, 200)

###one
fig = plt.figure(figsize=(28,11))
ax = fig.add_subplot(111)

plt.hist(csIC['components_mode_0/value'][:],bins = 100, histtype = 'stepfilled',ec='k',linewidth=8,color = 'lightsteelblue',alpha = 1)

plt.locator_params(axis='y', nbins=2)
plt.locator_params(axis='x', nbins=2)
fs = 100
ax.set_ylabel('Number of particles',size =fs)
ax.set_xlabel('Variability component',size=fs)

ax.tick_params(labelsize=fs)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(-200, 200)



















