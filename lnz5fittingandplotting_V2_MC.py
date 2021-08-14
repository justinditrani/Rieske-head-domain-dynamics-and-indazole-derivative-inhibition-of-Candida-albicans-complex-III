import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np

import os as os
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import matplotlib.cm as cm
import math
from matplotlib.ticker import LogLocator, AutoLocator


#MC simulation to get errors for fitting of inhibition data.
#please email justin.ditrani@gmail.com if there are any questions. 


def inhibcurve(It,Kd,kcat):
    b = -(Et + It + abs(Kd))
    a = 1
    c = Et*It
    
    EI = (-b - np.sqrt(np.square(b) - 4*a*c))/2

    Eapo = Et - EI

    return np.dot(Eapo,kcat)  

Et = 25e-9

dir_path = os.path.dirname(os.path.realpath(__file__))
Az = pd.read_csv(dir_path + '/bc1 activity_redo.csv',skiprows=1)


time = np.arange(0,15*61,15)
concInhib = np.dot([0, 0.655, 1.638, 4.096, 10.24, 25.6, 64, 160, 400, 1000, 2500],1e-9)

hdrs = list(Az.columns.values)

Avp = []
for xx in range(np.size(concInhib)):

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

SimC = inhibcurve(CI, *vals)[0:]
RMSD = np.sum(((SimC-R)**2)/np.size(CI))**(0.5)

NN = 10000
vv = []
for i in range(NN):
    Array = np.random.randn(np.size(SimC))*(RMSD)
    SimE = SimC + Array 
    
    vals, cov = curve_fit(inhibcurve,CI,SimE,p0 = [10e-9,1e9])
    vv = np.append(vv,vals[0])
    
SDEV = np.std(vv)
AVE = np.mean(vv)









