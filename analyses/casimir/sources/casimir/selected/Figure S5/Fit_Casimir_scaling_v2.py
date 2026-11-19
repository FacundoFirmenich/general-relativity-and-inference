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
Data2 = pd.read_hdf('Casimir_pressure_dist_v6_RRR_1.h5')
Data3 = pd.read_hdf('Casimir_pressure_dist_v6_RRR_10.h5')
Data4 = pd.read_hdf('Casimir_pressure_dist_v6_epsilonr_1.h5')
Data5 = pd.read_hdf('Casimir_pressure_dist_v6_epsilonr_1_06.h5')
Data6 = pd.read_hdf('Casimir_pressure_dist_v6_omegap_10eV.h5')
Data7 = pd.read_hdf('Casimir_pressure_dist_v6_omegap_12_5eV.h5')
Data8 = pd.read_hdf('Casimir_pressure_dist_v6_omegap_13_5eV.h5')
Data9 = pd.read_hdf('Casimir_pressure_dist_v6_omegap_15eV.h5')

initialguess = [22537944,3.2]

popt,pcov,perr,r_squared = fitCasfunc(Data.iloc[:,0],np.pi*Data.iloc[:,1],initialguess)
popt2,pcov2,perr2,r_squared2 = fitCasfunc(Data2.iloc[:,0],np.pi*Data2.iloc[:,1],initialguess)
popt3,pcov3,perr3,r_squared3 = fitCasfunc(Data3.iloc[:,0],np.pi*Data3.iloc[:,1],initialguess)
popt4,pcov4,perr4,r_squared4 = fitCasfunc(Data4.iloc[:,0],np.pi*Data4.iloc[:,1],initialguess)
popt5,pcov5,perr5,r_squared5 = fitCasfunc(Data5.iloc[:,0],np.pi*Data5.iloc[:,1],initialguess)
popt6,pcov6,perr6,r_squared6 = fitCasfunc(Data6.iloc[:,0],np.pi*Data6.iloc[:,1],initialguess)
popt7,pcov7,perr7,r_squared7 = fitCasfunc(Data7.iloc[:,0],np.pi*Data7.iloc[:,1],initialguess)
popt8,pcov8,perr8,r_squared8 = fitCasfunc(Data8.iloc[:,0],np.pi*Data8.iloc[:,1],initialguess)
popt9,pcov9,perr9,r_squared9 = fitCasfunc(Data9.iloc[:,0],np.pi*Data9.iloc[:,1],initialguess)

#Optimal values:
popt1 = [-6.95778690e+07, 3.19273178e+00]
perr1 = [4.04349563e+05, 2.37283016e-03]
r_squared1 = 0.9999791131853106

dist_forplot = np.linspace(Data.iloc[0,0],Data.iloc[-1,0],1001)
Casmodel = Casfunc(dist_forplot,*popt1)
#%%
fig1 = plt.figure(figsize = (6,3))
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Distance (nm)')
ax1.set_ylabel('Pressure (kPa)')

ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'Lifshitz formula (BCS)')
ax1.plot(Data2.iloc[:,0],np.pi*Data2.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'RRR = 1')
ax1.plot(Data3.iloc[:,0],np.pi*Data3.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'RRR = 10')
ax1.plot(Data4.iloc[:,0],np.pi*Data4.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\epsilon_r$ = 1')
ax1.plot(Data5.iloc[:,0],np.pi*Data5.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\epsilon_r$ = 1.06')
ax1.plot(Data6.iloc[:,0],np.pi*Data6.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\Omega_p$ = 10 eV')
ax1.plot(Data7.iloc[:,0],np.pi*Data7.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\Omega_p$ = 12.5 eV')
ax1.plot(Data8.iloc[:,0],np.pi*Data8.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\Omega_p$ = 13.5 eV')
ax1.plot(Data9.iloc[:,0],np.pi*Data9.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\Omega_p$ = 15 eV')
#ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,3]/1e3, marker = 'o', alpha = 0.6, label = 'Plasma model')
ax1.plot(dist_forplot,Casmodel/1e3, color = 'black', label = 'Fit $\propto d^{-3.1927}$')

ax1.legend(fontsize = 'small',frameon = False)
ax1.set_xlim((9,30))

#fig1.savefig('Plot_casimir_pressure_full_v2.png', dpi = 300, bbox_inches = 'tight')
#fig1.savefig('Plot_casimir_pressure_full_v2.svg', format = 'svg', bbox_inches = 'tight')

#%%
fig1 = plt.figure(figsize = (6,3))
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Distance (nm)')
ax1.set_ylabel('Pressure (kPa)')

ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'Nominal simulation (BCS)')
ax1.plot(Data2.iloc[:,0],np.pi*Data2.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'RRR = 1')
ax1.plot(Data3.iloc[:,0],np.pi*Data3.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'RRR = 10')
#ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,3]/1e3, marker = 'o', alpha = 0.6, label = 'Plasma model')
ax1.plot(dist_forplot,Casmodel/1e3, color = 'black', label = 'Fit $\propto d^{-3.1927}$')

ax1.legend(fontsize = 'small',frameon = False)
ax1.set_xlim((9,30))

fig1.savefig('Plot_casimir_pressure_RRR.png', dpi = 300, bbox_inches = 'tight')
fig1.savefig('Plot_casimir_pressure_RRR.svg', format = 'svg', bbox_inches = 'tight')

#%%
fig1 = plt.figure(figsize = (6,3))
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Distance (nm)')
ax1.set_ylabel('Pressure (kPa)')

ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'Nominal simulation (BCS)')
ax1.plot(Data4.iloc[:,0],np.pi*Data4.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$ϵ_r$ = 1')
ax1.plot(Data5.iloc[:,0],np.pi*Data5.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$ϵ_r$ = 1.06')
#ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,3]/1e3, marker = 'o', alpha = 0.6, label = 'Plasma model')
ax1.plot(dist_forplot,Casmodel/1e3, color = 'black', label = 'Fit $\propto d^{-3.1927}$')

ax1.legend(fontsize = 'small',frameon = False)
ax1.set_xlim((9,30))

fig1.savefig('Plot_casimir_pressure_epsilon.png', dpi = 300, bbox_inches = 'tight')
fig1.savefig('Plot_casimir_pressure_epsilon.svg', format = 'svg', bbox_inches = 'tight')


#%%
fig1 = plt.figure(figsize = (6,3))
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Distance (nm)')
ax1.set_ylabel('Pressure (kPa)')

ax1.plot(Data6.iloc[:,0],np.pi*Data6.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\Omega_p$ = 10 eV')
ax1.plot(Data7.iloc[:,0],np.pi*Data7.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\Omega_p$ = 12.5 eV')
ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'Nominal simulation (BCS)')
ax1.plot(Data8.iloc[:,0],np.pi*Data8.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\Omega_p$ = 13.5 eV')
ax1.plot(Data9.iloc[:,0],np.pi*Data9.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = r'$\Omega_p$ = 15 eV')
#ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,3]/1e3, marker = 'o', alpha = 0.6, label = 'Plasma model')
ax1.plot(dist_forplot,Casmodel/1e3, color = 'black', label = 'Fit $\propto d^{-3.1927}$')

ax1.legend(fontsize = 'small',frameon = False)
ax1.set_xlim((9,30))

fig1.savefig('Plot_casimir_pressure_omegap.png', dpi = 300, bbox_inches = 'tight')
fig1.savefig('Plot_casimir_pressure_omegap.svg', format = 'svg', bbox_inches = 'tight')
