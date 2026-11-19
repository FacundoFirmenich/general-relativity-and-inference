# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:14:44 2024

@author: dejonm1
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def parabola(dist,k,offs,yoff):
    return k*(dist-offs)**2 + yoff

def fitparabola(dist,data,initialguess):
    pmin = [0,0,-100]
    pmax = [100,35,10]
    popt,pcov = curve_fit(parabola,dist,data,initialguess,maxfev = 8000, bounds = (pmin,pmax))

    #Calculate the std errors in the parameters
    perr = np.sqrt(np.diag(pcov))
    
    #Calculate how good of a job the fitting did (r_squared)
    residuals = data - parabola(dist, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((data-np.mean(data))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return popt,pcov,perr,r_squared



ds = np.linspace(12,30,101) #Distances, in nm
xaxis = np.linspace(5,35,100001)
targetf = 10.01e6*np.ones(len(ds))
omegas = 2*np.pi*10.58e6*np.ones(len(ds))#2*np.pi*10.001e6 #Mechanical mode frequency

#These values are the result of a (rather long) manual optimization and interation process.
#They ensure that the final frequency of an oscillator with these starting frequencies will
#will end up at 2*pi*10.01 MHz (where we have measured it in our case)
#Works in the range of ds, so 12-30 nm.
omegas = omegas - 2*np.pi*np.array([-2.59887511e+07, -2.48705907e+07, -2.38033474e+07, -2.27830275e+07,
       -2.18080835e+07, -2.08756828e+07, -1.99833818e+07, -1.91291453e+07,
       -1.83109755e+07, -1.75267842e+07, -1.67751148e+07, -1.60543606e+07,
       -1.53627974e+07, -1.46991581e+07, -1.40619205e+07, -1.34501827e+07,
       -1.28625837e+07, -1.22982098e+07, -1.17559451e+07, -1.12346931e+07,
       -1.07339212e+07, -1.02527392e+07, -9.79040380e+06, -9.34550842e+06,
       -8.91856338e+06, -8.50789127e+06, -8.11356585e+06, -7.73465906e+06,
       -7.37097908e+06, -7.02189243e+06, -6.68668765e+06, -6.36519917e+06,
       -6.05693196e+06, -5.76159531e+06, -5.47852870e+06, -5.20738720e+06,
       -4.94801294e+06, -4.69996364e+06, -4.46270854e+06, -4.23612235e+06,
       -4.01967297e+06, -3.81342056e+06, -3.61632822e+06, -3.42884866e+06,
       -3.24998481e+06, -3.07965828e+06, -2.91756765e+06, -2.76358999e+06,
       -2.61697893e+06, -2.47753877e+06, -2.34505159e+06, -2.21917157e+06,
       -2.09961213e+06, -1.98602097e+06, -1.87820858e+06, -1.77570769e+06,
       -1.67835474e+06, -1.58585713e+06, -1.49805376e+06, -1.41456166e+06,
       -1.33524449e+06, -1.25986116e+06, -1.18818598e+06, -1.11992463e+06,
       -1.05504229e+06, -9.93362780e+05, -9.34537012e+05, -8.78564396e+05,
       -8.25209183e+05, -7.74398185e+05, -7.25921551e+05, -6.79653752e+05,
       -6.35538010e+05, -5.93397187e+05, -5.53039392e+05, -5.14684894e+05,
       -4.77862076e+05, -4.42698049e+05, -4.09127195e+05, -3.76932068e+05,
       -3.46039403e+05, -3.16504659e+05, -2.88194671e+05, -2.61049282e+05,
       -2.35011725e+05, -2.10028343e+05, -1.86048402e+05, -1.63023946e+05,
       -1.40877281e+05, -1.19631468e+05, -9.91822543e+04, -7.94944490e+04,
       -6.05905284e+04, -4.23519526e+04, -2.48296483e+04, -7.93937916e+03,
        8.34580683e+03,  2.40281483e+04,  3.91570724e+04,  5.37338559e+04,
        6.78038732e+04])

omegas = np.minimum(omegas,2*np.pi*5e7)

n = -3.19273178 #Extracted from fit to Casimir_pressure_dist_v3.h5, 95% confidence interval is +-2.37283016e-03
Caspres = -6.95778690e7/(xaxis**-n) #Fitted function to Casimir_pressure_dist_v3.h5, to cover any distance d (in nm).
#95% confidence interval in these the parameter is +- 4.04349563e+05
Casforce = Caspres*np.pi*(0.5*11.3e-6)**2
Caspotential = np.cumsum(Casforce*(xaxis[-1] - xaxis[0])*1e-9/len(xaxis))
Casplot = -Caspotential +Caspotential[-1]

xstable = np.zeros(len(ds))
popts = np.zeros((3,len(ds)))
perrs = np.zeros((3,len(ds)))

for idx,d in enumerate(ds):
    hbar = 1.0545718e-34
    kb = 1.38064852e-23
    omega_1 = omegas[idx]
    meff = 3.96e-14 #This is SI unit, which is kg. Louise calculated it based on the thickness/size of the drum + density of Al.
    xzpf = np.sqrt(hbar/(2*meff*omega_1))
    
    k = omega_1**2*meff
    
    mechpot = 0.5*k*(xaxis*1e-9-d*1e-9)**2
    
    totpot = (mechpot+Casplot)*1e15
    totpotdiff = np.diff(totpot)
    negdiff = np.where(totpotdiff < 0)[0]
    
    if len(negdiff) == 0:
        xstable[idx] = 0
    else:
        xind = np.argmin(totpot[negdiff]) + negdiff[0]
        xstable[idx] = xaxis[xind]
        if len(negdiff) > 300:
            stbpart = xaxis[xind - 100:xind + 100]
            stbvals = totpot[xind - 100:xind + 100]
            
            initialguess = [0.5*k*1e-3,xstable[idx],np.min(stbvals)]
            
            popts[:,idx],pcov,perrs[:,idx],r_squared = fitparabola(stbpart,stbvals,initialguess)

#%%
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Initial separation (nm)')
ax1.set_ylabel('Stable position (nm)')

ax1.plot(ds,xstable, linestyle = 'None', marker = '.', label = 'Simulated')
ax1.plot(ds,ds-xstable, linestyle = 'dashed', color = 'grey', label = 'Pos. shift')
ax1.plot(ds,ds, linestyle = 'dotted', color = 'black')

ax1.legend()

fig1.savefig('Stable_position_calculation_for_final_frequency.png', dpi = 300, bbox_inches = 'tight')
fig1.savefig('Stable_position_calculation_for_final_frequency.svg', format = 'svg', bbox_inches = 'tight')
#ax1.set_ylim(-10,10)

#%%
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

ax2.set_xlabel('Separation (nm)')
ax2.set_ylabel('Energy (fJ)')

ax2.plot(stbpart,stbvals)
pb1 = parabola(stbpart,0.5*k*1e-3,xstable[idx],np.min(stbpart))
pb2 = parabola(stbpart,*popts[:,idx])
ax2.plot(stbpart,pb1)
ax2.plot(stbpart,pb2)

#%%
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

ax3.set_xlabel('Initial separation (nm)')
ax3.set_ylabel('Frequency (Hz)')

finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

ax3.plot(ds,finalf, linestyle = 'None', marker = '.', label = 'Final frequency')
ax3.plot(ds,omegas/2/np.pi, linestyle = 'dotted', color = 'black', label = 'Start frequency')
ax3.plot(ds,targetf, linestyle = 'dashed', color = 'black', label = 'Target')
ax3.legend()
ax3.set_ylim(0,4e7)
#ax3.set_yscale('log')
#ax3.set_xscale('log')

fig3.savefig('Perturbed_versus_unperturbed_frequency_distance_2.png', dpi = 300, bbox_inches = 'tight')
fig3.savefig('Perturbed_versus_unperturbed_frequency_distance_2.svg', format = 'svg', bbox_inches = 'tight')

#%%
sumarr = np.vstack([ds,omegas,xstable])

np.savetxt("Optimized_mech_freqs.txt",sumarr)

#%%
test = np.loadtxt("Optimized_mech_freqs.txt")