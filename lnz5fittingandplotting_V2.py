# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 14:54:40 2021

@author: justi
"""



import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np

import os as os
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import matplotlib.cm as cm
import math
from matplotlib.ticker import LogLocator, AutoLocator

def inhibcurve(It,Kd,kcat):
    b = -(Et + It + abs(Kd))
    a = 1
    c = Et*It
    
    EI = (-b - np.sqrt(np.square(b) - 4*a*c))/2
    print(EI)
    Eapo = Et - EI
    print(kcat*Eapo)
    return np.dot(Eapo,kcat)  

Et = 25e-9

dir_path = os.path.dirname(os.path.realpath(__file__))
Az = pd.read_csv(dir_path + '/bc1 activity_V2.csv',skiprows=1)


time = np.arange(0,15*61,15)
concInhib = np.dot([0, 0.655, 1.638, 4.096, 10.24, 25.6, 64, 160, 400, 1000, 2500],1e-9)

hdrs = list(Az.columns.values)

Avp = []
for xx in range(np.size(concInhib)):
    print(xx)
    AA = (Az[hdrs[xx]] + Az[hdrs[xx+12]] + Az[hdrs[xx+24]]+Az[hdrs[xx+36]])/4
    Avp = np.append(Avp,AA)

Av = np.reshape(Avp,(11,60))
tt = 8

rates = []
rSTD = []
for xx in range(np.size(concInhib)):
    reg = LinearRegression().fit(np.reshape(time[0:tt],(-1,1)),Av[xx,0:tt])
    rates = np.append(rates,reg.coef_)

    regERR = LinearRegression().fit(np.reshape(time[0:tt],(-1,1)),Az[hdrs[xx]][0:tt])
    rER1 = regERR.coef_
    regERR = LinearRegression().fit(np.reshape(time[0:tt],(-1,1)),Az[hdrs[xx+12]][0:tt])
    rER2 = regERR.coef_
    regERR = LinearRegression().fit(np.reshape(time[0:tt],(-1,1)),Az[hdrs[xx+24]][0:tt])
    rER3 = regERR.coef_
    
    rSTD = np.append(rSTD,np.std([rER1,rER2,rER3]))
    
BRreg = LinearRegression().fit(np.reshape(time[0:tt],(-1,1)),Br['H10'][0:tt])
BR = BRreg.coef_

BRreg = LinearRegression().fit(np.reshape(time[0:tt],(-1,1)),Br['H10'][0:tt])
BR1 = BRreg.coef_
BRreg = LinearRegression().fit(np.reshape(time[0:tt],(-1,1)),Br['H11'][0:tt])
BR2 = BRreg.coef_
BRreg = LinearRegression().fit(np.reshape(time[0:tt],(-1,1)),Br['H12'][0:tt])
BR3 = BRreg.coef_

BRstd = np.std([BR1,BR2,BR3])

CI = concInhib[0:]
R =(rates[0:]-BR)/(rates[0]-BR)
vals, cov = curve_fit(inhibcurve,CI,R,p0 = [10e-9,1e9])

color_list = cm.tab10(np.linspace(0, 1, 11))

ERRs = np.sqrt(((rSTD[0:]**2+BRstd**2)/(rates[0:]-BR)**2) + (rSTD[0]/rates[0])**2)*R*5

cc=1
fig = plt.figure()
ax = fig.add_subplot(1,1, 1)
plt.errorbar(CI[1:]*1e9, R[1:]*5, ERRs[1:], fmt='.k',capsize=8,capthick=3)

for x,y,z in zip(CI[1:]*1e9,R[1:],color_list[0:]):
    plt.scatter(x,y*5,color=z,s=100,zorder=9)
    cc = cc+1

CIsmooth = np.arange(0,5e-6,3.125e-10)
plt.plot(CIsmooth[2:]*1e9,inhibcurve(CIsmooth, *vals)[2:]*5,color = 'k',linewidth = 3,zorder=10)

dom = max(CI)-min(CI)
ra = max(R)-min(R)
#plt.ylim(min(R)-ra/10,max(R)+ra/10)
#plt.xlim
plt.ylabel('Activity (s⁻¹)',fontsize = 40)
plt.xlabel("lnz-5 (nM)",fontsize = 40)

plt.xticks(fontsize=32)
plt.yticks(fontsize=32)

#ax.set_xscale('log', linthreshx=1e-8,linscalex=1e-8)
ax.set_xscale('log')
#plt.xlim(0.6,3000)
plt.minorticks_on()
xaxis = plt.gca().xaxis
xaxis.set_minor_locator(LogLocator(subs=np.arange(0, 10)))
#plt.gca().set_xticks([0, 10, 100 , 1000])
plt.tick_params(which=  'minor',length = 8, width = 1)
plt.tick_params(which=  'major',length = 14, width = 2)
##insettime
ext = 10
fig = plt.figure()
ax = fig.add_subplot(1,1, 1)

for y,c,SL in zip(Av[1:],color_list[0:],rates[0:]):
    plt.plot(time[0:tt+ext],y[0:tt+ext]-y[0],'o',color=c,ms=25)
    plt.plot(time[0:tt+ext],time[0:tt+ext]*SL,'-',color='k',linewidth = 10)
    
plt.ylabel("$A_{550}$",fontsize = 100)
plt.xlabel("time (s)",fontsize = 100)

plt.xticks([0,250],fontsize=80)
plt.yticks([0,0.08],fontsize=80)    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
