# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 13:58:22 2024

@author: dejonm1
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def Casfunc(dist,a,n):
    return a/(dist**n)

def fitCasfunc(dist,data,initialguess):
    
    popt,pcov = curve_fit(Casfunc,dist,data,initialguess,maxfev = 20000)

    #Calculate the std errors in the parameters
    perr = np.sqrt(np.diag(pcov))
    
    #Calculate how good of a job the fitting did (r_squared)
    residuals = data - Casfunc(dist, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((data-np.mean(data))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return popt,pcov,perr,r_squared

Data = pd.read_hdf('Casimir_pressure_dist_v4.h5')

initialguess = [22537944,3.2]

popt,pcov,perr,r_squared = fitCasfunc(Data.iloc[:,0],np.pi*Data.iloc[:,1],initialguess)

#Optimal values:
popt1 = [-6.95778690e+07, 3.19273178e+00]
perr1 = [4.04349563e+05, 2.37283016e-03]
r_squared1 = 0.9999791131853106

dist_forplot = np.linspace(Data.iloc[0,0],Data.iloc[-1,0],1001)
Casmodel = Casfunc(dist_forplot,*popt1)
#%%
fig1 = plt.figure(figsize = (6,3))
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Plate separation (nm)')
ax1.set_ylabel('Pressure (kPa)')

ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'Lifshitz formula (BCS)')
#ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,3]/1e3, marker = 'o', alpha = 0.6, label = 'Plasma model')
ax1.plot(dist_forplot,Casmodel/1e3, color = 'black', label = 'Fit $\propto d^{-3.1927}$')

ax1.legend(fontsize = 'small',frameon = False)
ax1.set_xlim((9,30))

fig1.savefig('Plot_casimir_pressure_full_v2.png', dpi = 300, bbox_inches = 'tight')
fig1.savefig('Plot_casimir_pressure_full_v2.svg', format = 'svg', bbox_inches = 'tight')