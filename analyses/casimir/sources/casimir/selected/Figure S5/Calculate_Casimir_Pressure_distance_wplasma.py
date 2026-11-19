# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 16:30:07 2023

@author: dejonm1
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad

#I'll attempt to numerically calculate the Casimir pressure, based on Jurgen's thesis and whatever other references are necessary.
#Introduce some constants/material parameters
#Also based on Bimonte2019

#We are not using these constants, because it's based on a theory paper and fucking theorists use natural units without writing it down :<
epsilon_0 = 8.854e-12 #in F/m, dielectric permittivity of vacuum
kB = 1.381e-23 #Boltzmann constant, m^2 kg / (s^2 K) 
hbar = 1.055e-34 #Reduced planck constant, J*s
eV = 1.602e-19 #Electronvolt, in J
c = 2.9979e8 #Light speed, in m/s

#Material parameters for aluminium
Omega_p = 13 #in eV/hbar, plasma frequency of Aluminium, from Handbook of Optical Constants of Solids, E.D. Palik, 1985, via Bimonte 2019
gamma_0 = 0.1 #in eV/hbar, Phenomenological relaxation rate (dissipation of current, through collisions of electrons), source idem.
RRR = 1 #Residual resistance ratio, resistivity at 300K / resistivity just above superconducting transition
gamma = gamma_0/RRR
Tc = 1.2 #critical temperature of Aluminium, in K
epsilon_r = 1.03 #relative permittivity of Aluminium, taking into account the core interband transitions, based on Handbook of Optical Constants of Solids, E.D. Palik, 1985, via Bimonte 2019

#Calculate the zero-temperature gap.
c1 = 1.764
c2 = 0.9963
c3 = 0.7735
#Delta_0 = c1*kB*Tc*np.sqrt(1 - 0/Tc)*(c2 + c3*0/Tc) /eV
Delta_0 = 1.76*kB*Tc /eV

#To translate the superconducting gap energy into photon frequency:
#Delta_0*eV / (hbar 2 np.pi)
#For Al, this corresponds to about 44 GHz

#According to the literature, we only have to evaluate the system at the Matsubara frequencies
#These frequencies depend on the temperature
a = 40e-9 #Nominal gap spacing, 40 nm
T = 0.25 #Temperature 250 mK, this is where the highest difference of pressure exists
Delta = c1*kB*Tc*np.sqrt(1 - T/Tc)*(c2 + c3*T/Tc) /eV #Actual gap energy (temperature-dependent)
lmax = 500 #Maximum Matsubara frequency to evaluate.
l = np.linspace(0,lmax,lmax+1)
ksil = 2*np.pi*kB*T*l/eV 
#At zero frequency, the calculation for g converges but it's difficult to numerically integrate
#So pick a value slightly higher than zero.
ksil[0] = 2*np.pi*kB*T*0.001/eV

def integrandofg(varepsilon,ksi,T):
    Tc = 1.2 #critical temperature of Aluminium, in K
    kB = 1.381e-23 #Boltzmann constant, m^2 kg / (s^2 K) 
    eV = 1.602e-19 #Electronvolt, in J
    
    #Integral only valid for T<Tc, else it should return zero bc heavyside function
    if T<Tc:
        #Temperature dependent gap, known from BCS theory
        #Source: Tinkham, Introduction to superconductivity, 1996
        #constants
        c1 = 1.764
        c2 = 0.9963
        c3 = 0.7735
        Delta = c1*kB*Tc*np.sqrt(1 - T/Tc)*(c2 + c3*T/Tc)/eV #Gap is expressed in eV for some theory reason.
    
        #Definitions
        z = 1j*ksi
        
        #Eqs B2-B4 from Jurgen's thesis, or A2-A6 from Bimonte's 2019 paper.
        E = np.sqrt(varepsilon**2 + Delta**2)
        Qplus = np.sqrt((E+z)**2 - Delta**2)
        Aplus = E*(E+z) + Delta**2
        
        Gplus = (varepsilon**2 * Qplus + (Qplus + 1j*gamma)*Aplus) / (Qplus * (varepsilon**2 - (Qplus + 1j*gamma)**2))
    
        #output the integrand of the integral (A1) in Bimonte's 2019 paper, B1 in Jurgen's thesis.
        integrand = np.tanh(E/(2*kB*T))*np.real(Gplus)/E
    else:
        integrand = 0    
    return integrand

def integrandofp(kort,ksil,epsilon_l,a):
    c = 2.9979e8 #Light speed, in m/s
    #a = 40e-9 #Spacing, in nm
    
    s_l = np.sqrt(epsilon_l * ksil**2 / c**2 + kort**2)
    q_l = np.sqrt(ksil**2 / c**2 + kort**2)
    rte = (q_l - s_l)/(q_l + s_l)
    rtm = (epsilon_l * q_l - s_l)/(epsilon_l * q_l + s_l)
    
    #Throws issues if rte goes to 0, take elementwise maximum with small value instead
    #rte = np.maximum(rte,1e-30)
    
    #integrand = kort*q_l*(1/(np.exp(2*a*q_l)/(rte**2) - 1))
    #integrand = kort*q_l*(1/(np.exp(2*a*q_l)/(rtm**2) - 1))
    integrand = -kB*T * kort*q_l*(1/(np.exp(2*a*q_l)/(rte**2) - 1) + 1/(np.exp(2*a*q_l)/(rtm**2) - 1)) / np.pi
    return integrand

#%%
As = np.linspace(10e-9,100e-9,101)

p_bcs = np.zeros(len(As))
p_drude = np.zeros(len(As))
p_plasma = np.zeros(len(As))

lmax = 200000 #Maximum Matsubara frequency to evaluate.
l = np.linspace(0,lmax,lmax+1)

#Calculate Matsubara frequencies
ksil = 2*np.pi*kB*T*l/eV 
#At zero frequency, the calculation for g converges but it's difficult to numerically integrate
#So pick a value slightly higher than zero.
ksil[0] = 2*np.pi*kB*T*0.005/eV

#Initialize empty arrays to contain the integrated value of g, and the error
gval = np.zeros(len(ksil))
gerr = np.zeros(len(ksil))
for idx,ksi in enumerate(ksil):
    
    if T < Tc:
        #Temperature-dependent gap
        Delta = c1*kB*Tc*np.sqrt(1 - T/Tc)*(c2 + c3*T/Tc) /eV
        #Heuristic, sharp feature in g scales with sqrt(ksi*Delta).
        #See also https://iopscience.iop.org/article/10.1088/1751-8113/43/14/145304/pdf
        #I tested the convergence, for this bounds the error seems small enough
        #for all 0.00001*Delta < ksi < 1000 Delta_0     
        bound = 300*np.sqrt(ksi*Delta)
        gval[idx],gerr[idx] = quad(integrandofg,-bound,bound, args=(ksi,T),limit = 50)
    else:
        gval[idx] = 0

#Calculate the permittivities from BCS theory and from the Drude model
epsilon_bcs = epsilon_r + (Omega_p**2/ksil) * (1 / (ksil + gamma) + gval/(ksil)) #Eq 8 of Bimonte 2019
epsilon_drude = epsilon_r + Omega_p**2 / (ksil * (ksil + gamma)) #Eq 5 of Bimonte 2019
epsilon_plasma = epsilon_r + Omega_p**2  / ksil**2

#Initialize empty arrays to contain the integraded value of p, and the error
pbcs = np.zeros(len(ksil))
perr1 = np.zeros(len(ksil))
pdrude = np.zeros(len(ksil))
perr2 = np.zeros(len(ksil))
pplasma = np.zeros(len(ksil))
perr3 = np.zeros(len(ksil))

for idx2,a in enumerate(As):
    #iterate over all the Matsubara frequencies
    for idx,l in enumerate(ksil):
        #Rescale the frequencies into actual Hz
        ksi = ksil[idx]*eV/hbar
        #Evaluate the pressures for the BCS and Drude models separately
        #That is, integrate over all different momenta k_ort
        #Bounds are somewhat heuristic, integral should cut off well before plasema frequencies
        pbcs[idx],perr1[idx] = quad(integrandofp,1e-9,1e9,args = (ksi,epsilon_bcs[idx],a), limit = 50)
        pdrude[idx],perr2[idx] = quad(integrandofp,1e-9,1e9,args = (ksi,epsilon_drude[idx],a), limit = 50)
        pplasma[idx],perr3[idx] = quad(integrandofp,1e-9,1e9,args = (ksi,epsilon_plasma[idx],a), limit = 50)
        #Give half weight to the zero Matsubara frequency for theory reasons, see
        #Bimonte 2019
        if idx == 0:
            pbcs[idx] = 0.5*pbcs[idx]
            pdrude[idx] = 0.5*pdrude[idx]
            pplasma[idx] = 0.5*pplasma[idx]
    
    #Sum over all the Matsurbara frequencies    
    p_bcs[idx2] = np.sum(pbcs)
    p_drude[idx2] = np.sum(pdrude)
    p_plasma[idx2] = np.sum(pplasma)

#%%

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Gap (nm)')
ax1.set_ylabel('$\Delta$P (mPa)')

ax1.plot(As,(p_bcs - p_drude)*1e3)
ax1.plot(As,(p_bcs - p_plasma)*1e3)

#%%
Casmodel = -155.22*230e-26/((As)**3.2)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

ax2.set_xlabel('Gap (nm)')
ax2.set_ylabel('Pressure (Pa)')

ax2.plot(As*1e9,p_drude, label = 'Drude')
ax2.plot(As*1e9,p_plasma, label = 'Plasma')
ax2.plot(As*1e9,p_bcs, label = 'BCS', linestyle = 'dashed')
#ax2.plot(As*1e9,Casmodel)

#ax2.set_xlim(0,10)
#ax2.set_xscale('log')
#ax2.set_yscale('log')

ax2.plot()
ax2.legend()
#ax2.set_ylim((-600,0))

#%%
dataCas = pd.DataFrame(data = None, index = None, columns = None)

dataCas['Gap (nm)'] = As*1e9
dataCas['P Drude (Pa)'] = p_drude
dataCas['P BCS (Pa)'] = p_bcs
dataCas['P Plasma (Pa)'] = p_plasma

dataCas.to_hdf('Casimir_pressure_dist_v6_RRR_1.h5', key = 'data')