
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 08:50:22 2020

@author: justi
"""

#This program finds correclations between vectrocs from a PCA job (P26_J1582) using the symmetry expanded particle stacks (P26_J877)
#please email justin.ditrani@gmail.com if there are any questions. 

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import pandas
import operator
import collections
import os as os
import statistics as stats
import joypy
import matplotlib.cm as cm

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

csCB = open_cs('cryosparc_P26_J1582_particles.cs')

csCB_H, csCB_A = H_A(csCB)

csA_1 = csALL[range(0,len(csALL),2)]
csA_2 = csALL[range(1,len(csALL),2)]

UUA1 = np.where(np.isin(csA_1['uid'][:],csCB_A['uid']).astype(int) == 1)
UUA2 = np.where(np.isin(csA_2['uid'][:],csCB_A['uid']).astype(int) == 1)
UUCBA12 = np.intersect1d(UUA1,UUA2)

UUcb1 = np.where(np.isin(csCB_A['uid'][:],csA_1['uid'][UUCBA12]).astype(int) == 1)
CB_A1 = csCB_A[UUcb1]
UUcb2 = np.where(np.isin(csCB_A['uid'][:],csA_2['uid'][UUCBA12]).astype(int) == 1)
CB_A2 = csCB_A[UUcb2]

PP1 = np.where(CB_A1['components_mode_0/value']>0)
csBC_CA1 = CB_A1[PP1]
PP2 = np.where(CB_A1['components_mode_0/value']<0)
csBC_BA1 = CB_A1[PP2]


UUCA1 = np.where(np.isin(csA_1['uid'][:],csBC_CA1['uid']).astype(int) == 1)
UUCA2 = np.where(np.isin(csCB_A['uid'][:],csA_2['uid'][UUCA1]).astype(int) == 1)
h_CBc_A = csCB_A[UUCA2]

UUBA1 = np.where(np.isin(csA_1['uid'][:],csBC_BA1['uid']).astype(int) == 1)
UUBA2 = np.where(np.isin(csCB_A['uid'][:],csA_2['uid'][UUBA1]).astype(int) == 1)
h_CBb_A = csCB_A[UUBA2]

#transparency
tp = 0.8
fs = 36
lw = 8

fig = plt.figure(figsize=(28,11))
ax = fig.add_subplot(111)


N, bins, patches = plt.hist(CB_A1['components_mode_0/value'],100, histtype = 'stepfilled',ec='k',linewidth=8,color = 'lightsteelblue',alpha = tp)

for bb, p in zip(bins,patches):
    print(bb)
    if bb > 0:
        p.set_facecolor('rebeccapurple')
    else:
        p.set_facecolor('tab:red')

plt.locator_params(axis='y', nbins=2)
plt.locator_params(axis='x', nbins=2)
fs = 100
ax.set_ylabel('Number of particles',size =fs)
ax.set_xlabel('Variability component',size=fs)

ax.tick_params(labelsize=fs)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(-200, 200)


fig = plt.figure(figsize=(28,11))
ax = fig.add_subplot(111)

N, bins, patches = plt.hist(CB_A1['components_mode_0/value'],100, histtype = 'stepfilled',ec='k',linewidth=8,color = 'lightsteelblue',alpha = tp)


for bb, p in zip(bins,patches):
    print(bb)
    if bb < 0:
        p.set_facecolor('rebeccapurple')
    else:
        p.set_facecolor('tab:red')

plt.locator_params(axis='y', nbins=2)
plt.locator_params(axis='x', nbins=2)
fs = 100
ax.set_ylabel('Number of particles',size =fs)
ax.set_xlabel('Variability component',size=fs)

ax.tick_params(labelsize=fs)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(-200, 200)

fs = 60



APO_CB = {'..': h_CBb_A['components_mode_0/value'],'.': h_CBc_A['components_mode_0/value']}


fig, axes = joypy.joyplot(APO_CB,hist="True",bins = 100, histtype = 'stepfilled',grid='both',density='True',fill=True,linewidth=lw,overlap=-0.1,alpha = tp,x_range=[-150,150],figsize=(10,14),color=['tab:red'])

plt.rc("font", size=fs)
plt.xlabel('$\mathit{b/c1}$ variability component',  fontsize=fs)
plt.text(50, 0.94, "$\mathit{c1-state}$",fontsize=fs,color='tab:red',alpha = 0.7)
plt.text(-180, 0.94, "$\mathit{b-state}$",fontsize=fs,color='tab:red',alpha = 0.7)





fig, axes = joypy.joyplot(APO_CB,hist="True",bins = 100, histtype = 'stepfilled',grid='both',density='True',fill=True,linewidth=lw,overlap=-0.1,alpha = tp,x_range=[-150,150],figsize=(10,14),color=['rebeccapurple'])

plt.rc("font", size=fs)
plt.xlabel('$\mathit{b/c1}$ variability component',  fontsize=fs)
plt.text(50, 0.94, "$\mathit{c1-state}$",fontsize=fs,color='tab:red',alpha = 0.7)
plt.text(-180, 0.94, "$\mathit{b-state}$",fontsize=fs,color='tab:red',alpha = 0.7)












