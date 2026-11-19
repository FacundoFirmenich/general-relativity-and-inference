# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 11:44:25 2024

@author: dejonm1
"""

import numpy as np
import matlab.engine
import matplotlib.pyplot as plt

#%%

#Initialize Matlab engine. (takes some time)
eng = matlab.engine.start_matlab()
eng.cd(r'C:\LocalUserData\User-data\dejonm1\MatCont\MatCont7p4', nargout = 0)

#Pre-calculated 'initial frequencies' that result in a final frequency of 10.01 MHz (as observed). 
#This file also contains the spacing offset for the local minimum, necessary for MatCont
#See 'Precalculate_mechfreqs.py' and v2 of that script as well.
init_freq = np.loadtxt("Optimized_mech_freqs.txt")

#Some parameters which are fixed in this system, because they can be measured or calculated 
#(or cannot feasibly be measured more accurately, like the effective mass)
hbar = 1.055e-34 #Planck's reduced constant in SI units
omega_l = 2*np.pi*5.461831e9 #Cavity frequency (in Hz), typically 'laser frequency' in optomechanics.
meff = 3.96e-14 #In kg, from COMSOL simulation
gamma_m = 2*np.pi*170 #in Hz, Mechanical damping/loss rate (without Casimir force)
area = np.pi*(0.5*11.3e-6)**2 #Area of the drum, from Louise's design file and verified by AFM
omega_1 = 2*np.pi*10.001e6 #Measured mechanical mode frequency (thus includes Casimir effects)
kappa_e = 2*np.pi*251e3 #External cavity decay rate (in Hz), from cavity calibration
kappa_i = 2*np.pi*297e3 #Internal cavity decay rate (in Hz), from cavity calibration
kappa = kappa_e + kappa_i #Cavity decay rate
delta = 0 #Detuning cavity-drive (in Hz)
g0_1 = 2*np.pi*150 #Optomechanical coupling rate (Hz), from sideband coupling calibration
xzpf = np.sqrt(hbar/(2*meff*omega_1)) #Zero-point fluctuations of the mechanical resonator.
l_power = 1e-12 #Drive power at cavity frequency, in Watt
l_power2 = 2e-13 #Drive power at red sideband, in Watt
drive = np.sqrt(l_power*kappa_e / (hbar*omega_l))
drive_2 = np.sqrt(l_power2*kappa_e / (hbar*omega_l))

#I want only one parameter in the whole system. This is the initial separation between the two plates,
#Without considering the Casimir force. It is the offset between the mechanical zero and the other plate.
d = 18 #Distance between two drums in nm

omega_m = np.interp(d,init_freq[0,:],init_freq[1,:]) #in Hz, Mechanical frequency (without Casimir force)
#omega_m = 2*np.pi*16.3e6 #in Hz, Mechanical frequency (without Casimir force)

xinit = np.interp(d,init_freq[0,:],init_freq[2,:])-d #Initial guess for offset value of top plate position due to Casimir force, in nm
#xinit = -2.89 #Initial guess for offset value of top plate position due to Casimir force

#Caspres = 2175.904 #Casimir pressure from BCS model in Pascal
n = -3.19273178 #Extracted from fit to Casimir_pressure_dist_v3.h5, 95% confidence interval is +-2.37283016e-03
Caspres = 6.95778690e7/(d**-n) #Fitted function to Casimir_pressure_dist_v3.h5, to cover any distance d (in nm).
#95% confidence interval in this parameter is +- 4.04349563e+05

#%%

#Initialize array 
opt_dist = np.zeros(len(init_freq[0,:]))

#For every distance d in that we want to consider, 
#Take the frequency and initial guess for the steady position
#And then let it converge to the 'final' steady state in Matcont.
for idx,d in enumerate(init_freq[0,:]):
    omega_m = init_freq[1,idx]
    xinit = init_freq[2,idx] - d
    
    #Caspres = 2175.904 #Casimir pressure from BCS model in Pascal
    n = -3.19273178 #Extracted from fit to Casimir_pressure_dist_v3.h5, 95% confidence interval is +-2.37283016e-03
    Caspres = 6.95778690e7/(d**-n) #Fitted function to Casimir_pressure_dist_v3.h5, to cover any distance d (in nm).
    
    #Test the convergence of the solution in MatCont for the parameters (mainly Caspres and xinit)
    retval = eng.Test_convergence(omega_m,meff,gamma_m,n,d,Caspres,area,xinit)
    Conv_vals = np.array(retval._data).reshape((4,11))
    opt_dist[idx] = Conv_vals[0,-1] #Revised offset value for the top plate position due to the Casimir force, in nm
    
#%%
recalc_freq = np.vstack([init_freq,opt_dist + init_freq[0,:]])
np.savetxt("Optimized_mech_freqs_v2.txt",recalc_freq)

fig1= plt.figure()
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Initial separation (nm)')
ax1.set_ylabel('Stable position (nm)')

ax1.plot(recalc_freq[0,:],recalc_freq[0,:], linestyle = 'dotted', color = 'black')
ax1.plot(recalc_freq[0,:],recalc_freq[2,:], marker = '.', linestyle = 'None', label = 'Precalculated')
ax1.plot(recalc_freq[0,:],recalc_freq[3,:], marker = '.', linestyle = 'None', label = 'Recalculated')

ax1.legend()

fig1.savefig('Comparison_precalculated_recalculated_stable_positions.png', dpi = 300, bbox_inches = 'tight')
