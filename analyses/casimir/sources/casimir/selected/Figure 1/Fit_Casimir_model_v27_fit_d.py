# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:35:01 2024

@author: dejonm1
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import json
from scipy.optimize import curve_fit
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import bisplrep, bisplev
from odeintw import odeintw
from scipy.fft import fftfreq,fftshift,fft
import matlab.engine
from matplotlib import gridspec

def Lorentzian(x, I, x0, gamma):
    #Cauchy-distribution (Lorentzian power spectrum)
    return I / ( 1 + ( ( x-x0 ) / gamma ) ** 2 )

def Powerfactor(omega_1,kappa,delta,drive,g0_1,drive_2,xzpf):
    numamps = 101

    xmotions = np.logspace(-12,-9,numamps)
    avcavpow = np.zeros(numamps)#Average cavity power
    #avcaplot = np.zeros((2000,numamps)) #uncomment to get the cavity field time dependence

    #Iterate over all different motion amplitudes.
    for idx,xmotion in enumerate(xmotions):
        #Initial conditions for the simulation
        x0 = -xmotion/xzpf + 0j
        p0 = 0 + 0j
        a0 = 0+100j
        
        #Load parameters and initial conditions
        params = [omega_1, 0, kappa, delta, drive, g0_1, drive_2, 0]
        coords0 = [x0,p0,a0]
        
        #Set up parameters of the simulation
        npts = 160000
        deltat = 1/(100*omega_1)
        t = np.linspace(0.,npts*deltat,npts)

        #Run the simulation
        psol1 = odeintw(Cavandonedrum2,coords0,t,args = (params,))

        #Take the average cavity power (over the last mechanical periods)
        lbpow = np.abs((psol1[-2000:,2] + np.conj(psol1[-2000:,2]))/2)**2
        avcavpow[idx] = np.mean(lbpow)
        #avcaplot[:,idx] = lbpow #Uncomment to get the cavity field time dependence
        
    return np.array([xmotions,avcavpow/avcavpow[0]])

def Cavandonedrum2(coords,t,params):
    #Coordinates of two membranes and one optical field
    x1,p1,a = coords
    #All the parameters for this system except for noise.
    omega_1, gamma_1, kappa, delta, drive, g0_1, drive_2, VNAdiff = params
    
    #Time derivatives,
    derivs = [omega_1*p1,
              -omega_1*x1-gamma_1*p1-g0_1*np.conj(a)*a,
              -(1j*delta + kappa/2)*a + 1j*g0_1*x1*a - drive/np.sqrt(2) - drive_2/np.sqrt(2) * np.exp(-1j*(omega_1+VNAdiff)*t)]
    return derivs

def SimCas_w_prescribedamplitude(d,ampl_opto):
    Fext = 0.1 #Initial value for Fext (driving force)
    
    meff = 3.96e-14 #In kg, from COMSOL simulation
    gamma_m = 2*np.pi*140#2*np.pi*170 #in Hz, Mechanical damping/loss rate (without Casimir force)
    area = np.pi*(0.5*11.3e-6)**2 #Area of the drum, from Louise's design file and verified by AFM
    omega_m = np.interp(d,init_freq[0,:],init_freq[1,:]) #in Hz, Mechanical frequency (without Casimir force)
    xinit = np.interp(d,init_freq[0,:],init_freq[2,:])-d #Initial guess for offset value of top plate position due to Casimir force, in nm
    n = -3.19273178 #Extracted from fit to Casimir_pressure_dist_v3.h5, 95% confidence interval is +-2.37283016e-03
    Caspres = 6.95778690e7/(d**-n) #Fitted function to Casimir_pressure_dist_v3.h5, to cover any distance d (in nm).
    
    #Initial simulation of the Casimir force for the initial Fext guess.
    output = eng.Casimir_single(omega_m,meff,gamma_m,n,d,Caspres,area,xinit,Fext*1e-12)
    cascurve = np.array(output._data).reshape((2,int(output.size[1]/2)))
    
    #Determine the maximum amplitude
    ampl_cas = np.max(cascurve[0,:])
    
    #If the maximum Casimir amplitude is different from the maximum amplitude from the optomechanics simulation
    while np.abs(ampl_opto/ampl_cas - 1) > 1e-5:
        #Re-estimate Fext
        if np.abs(ampl_opto/ampl_cas - 1) < 1e-2:
            Fext = ((ampl_opto/ampl_cas - 1)*0.8 + 1) *Fext
        else:
            Fext = (ampl_opto/ampl_cas)*Fext  
        print('Fext: {:.9f} pN'.format(Fext))
        
        #Re-do the simulation
        output = eng.Casimir_single(omega_m,meff,gamma_m,n,d,Caspres,area,xinit,Fext*1e-12)
        cascurve = np.array(output._data).reshape((2,int(output.size[1]/2)))
        
        #Re-determine the maximum amplitude
        ampl_cas = np.max(cascurve[0,:])
    
    return cascurve


#%%

#import all the data 

mat = scipy.io.loadmat('2024-03-19-1759.mat')
filtfreqs = np.load("2024-03-19-1759_filtfreqs.npy")
filtdata_up = np.load("2024-03-19-1759_filtdata_up.npy")
filtdata_down = np.load("2024-03-19-1759_filtdata_down.npy")

freqs = mat["freq"]

#Dim 0 = all frequencies
#Dim 1 = which sideband (red,blue,2red,2blue,3red,3blue)
#Dim 2 = VNA power
#Dim 3 = MW power
#Dim 4 = Navg
Data_down = mat["data_down"] #Frequency swept up
Data_up = mat["data_up"] #Frequency swept down

MWpowers = mat["MWpowers"]
VNApowers = mat["VNApowers"]

VNApoints = mat["numfreqs"][0,0]
VNAspan = mat["VNAspan"][0,0]

VNAfreqsred = mat["allfreqs_red"]
VNAfreqsblue = mat["allfreqs_blue"]
fcenters = mat["frequency_centers"]

# --- #

mat2 = scipy.io.loadmat('2024-03-21-1025.mat')
filtfreqs2 = np.load("2024-03-21-1025_filtfreqs.npy")
filtdata_up2 = np.load("2024-03-21-1025_filtdata_up.npy")
filtdata_down2 = np.load("2024-03-21-1025_filtdata_down.npy")

freqs2 = mat2["freq"]

#Dim 0 = all frequencies
#Dim 1 = which sideband (red,blue,2red,2blue,3red,3blue)
#Dim 2 = VNA power
#Dim 3 = MW power
#Dim 4 = Navg
Data_down2 = mat2["data_down"] #Frequency swept up
Data_up2 = mat2["data_up"] #Frequency swept down

MWpowers2 = mat2["MWpowers"]
VNApowers2 = mat2["VNApowers"]

VNApoints2 = mat2["numfreqs"][0,0]
VNAspan2 = mat2["VNAspan"][0,0]

VNAfreqsred2 = mat2["allfreqs_red"]
VNAfreqsblue2 = mat2["allfreqs_blue"]
fcenters2 = mat2["frequency_centers"]

# --- #

mat3 = scipy.io.loadmat('2024-03-21-1858.mat')
filtfreqs3 = np.load("2024-03-21-1858_filtfreqs.npy")
filtdata_up3 = np.load("2024-03-21-1858_filtdata_up.npy")
filtdata_down3 = np.load("2024-03-21-1858_filtdata_down.npy")
freqs3 = mat3["freq"]

#Dim 0 = all frequencies
#Dim 1 = which sideband (red,blue,2red,2blue,3red,3blue)
#Dim 2 = VNA power
#Dim 3 = MW power
#Dim 4 = Navg
Data_down3 = mat3["data_down"] #Frequency swept up
Data_up3 = mat3["data_up"] #Frequency swept down

MWpowers3 = mat3["MWpowers"]
VNApowers3 = mat3["VNApowers"]

VNApoints3 = mat3["numfreqs"][0,0]
VNAspan3 = mat3["VNAspan"][0,0]

VNAfreqsred3 = mat3["allfreqs_red"]
VNAfreqsblue3 = mat3["allfreqs_blue"]
fcenters3 = mat3["frequency_centers"]

# --- #

mat4 = scipy.io.loadmat('2024-03-22-0918.mat')
filtfreqs4 = np.load("2024-03-22-0918_filtfreqs.npy")
filtdata_up4 = np.load("2024-03-22-0918_filtdata_up.npy")
filtdata_down4 = np.load("2024-03-22-0918_filtdata_down.npy")

freqs4 = mat4["freq"]

#Dim 0 = all frequencies
#Dim 1 = which sideband (red,blue,2red,2blue,3red,3blue)
#Dim 2 = VNA power
#Dim 3 = MW power
#Dim 4 = Navg
Data_down4 = mat4["data_down"] #Frequency swept up
Data_up4 = mat4["data_up"] #Frequency swept down

MWpowers4 = mat4["MWpowers"]
VNApowers4 = mat4["VNApowers"]

VNApoints4 = mat4["numfreqs"][0,0]
VNAspan4 = mat4["VNAspan"][0,0]

VNAfreqsred4 = mat4["allfreqs_red"]
VNAfreqsblue4 = mat4["allfreqs_blue"]
fcenters4 = mat4["frequency_centers"]

#%%

#Initialize Matlab engine. (takes some time)
eng = matlab.engine.start_matlab()
eng.cd(r'C:\LocalUserData\User-data\dejonm1\MatCont\MatCont7p4', nargout = 0)

#Pre-calculated 'initial frequencies' that result in a final frequency of 10.01 MHz (as observed). 
#This file also contains the spacing offset for the local minimum, necessary for MatCont
#See 'Precalculate_mechfreqs.py' and v2 of that script as well.
#init_freq = np.loadtxt("Optimized_mech_freqs.txt")
#
#I put all the 'init_freqs' stable positions through MatCont to save some time. 
#With the revised numbers, I can just interpolate the initial position instead of having to recalculate for every function call.
init_freq = np.loadtxt("Optimized_mech_freqs_v2.txt")

#%%
#Some parameters which are fixed in this system, because they can be measured or calculated 
#(or cannot feasibly be measured more accurately, like the effective mass)
hbar = 1.055e-34 #Planck's reduced constant in SI units
omega_l = 2*np.pi*5.461831e9 #Cavity frequency (in Hz), typically 'laser frequency' in optomechanics.
meff = 3.96e-14 #In kg, from COMSOL simulation
gamma_m = 2*np.pi*140#2*np.pi*170 #in Hz, Mechanical damping/loss rate (without Casimir force)
area = np.pi*(0.5*11.3e-6)**2 #Area of the drum, from Louise's design file and verified by AFM
omega_1 = 2*np.pi*10.011e6 #Measured mechanical mode frequency (thus includes Casimir effects)
kappa_e = 2*np.pi*251e3 #External cavity decay rate (in Hz), from cavity calibration
kappa_i = 2*np.pi*297e3 #Internal cavity decay rate (in Hz), from cavity calibration
kappa = kappa_e + kappa_i #Cavity decay rate
delta = 0#2*np.pi*5e3#0 #Detuning cavity-drive (in Hz)
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

n = -3.19273178 #Extracted from fit to Casimir_pressure_dist_v3.h5, 95% confidence interval is +-2.37283016e-03
Caspres = 6.95778690e7/(d**-n) #Fitted function to Casimir_pressure_dist_v3.h5, to cover any distance d (in nm).
#95% confidence interval in this parameter is +- 4.04349563e+05

#%%

MWpower_sims = np.linspace(-40,-12.5,12)
VNApower_sims = np.linspace(-40,-10,13)

Peakpowers = np.loadtxt('Fit_Casimir_model_v27_Peakpowers.txt')
simpowers_VNAsweep_40 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-40dBm_MW.txt')
simpowers_VNAsweep_35 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-35dBm_MW.txt')
simpowers_VNAsweep_30 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-30dBm_MW.txt')
simpowers_VNAsweep_27_5 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-27_5dBm_MW.txt')
simpowers_VNAsweep_25 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-25dBm_MW.txt')
simpowers_VNAsweep_22_5 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-22_5dBm_MW.txt')
#simpowers_VNAsweep_20 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-20dBm_MW.txt')
simpowers_VNAsweep_20 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-20dBm_MW_altfreq.txt')
simpowers_VNAsweep_17_5 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-17_5dBm_MW.txt')
simpowers_VNAsweep_15 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-15dBm_MW.txt')

#%%
VNApowers_inW = 10**((-30+VNApower_sims-22.2-47.4)/10)
MWpowers_inW = 10**((-30+MWpower_sims-3-47.4)/10)
#
PF_file = np.loadtxt("Power_factors_9_12e-11.txt")
xmotions = PF_file[0,:]
lpower2s = PF_file[1,:]
All_PF = PF_file[2:,:]

#Somewhat complicated rescaling is necessary; All_PF was calculated for MW power of 1e-12.
#So take the MW power of the current run, in W, divide it by 1e-12, and then divide the relevant lpower2s entry by that value.
#This gives you the correct 'adjusted VNA power' to use (the power factor only accounts for the relative powers between the two)
#e.g. for -20 dBm MW drive, MWpowers_inW[9]/1e-12 = 161.18
#Then for the -20 dBm VNA drive ( = 1.096e-12 W), we should actually take the 6.7608e-15 W lpower2s curve. 

#%%

Fext = 0.130915
d = 18 #Distance between two drums in nm
omega_m = np.interp(d,init_freq[0,:],init_freq[1,:]) #in Hz, Mechanical frequency (without Casimir force)
#omega_m = 2*np.pi*16.3e6 #in Hz, Mechanical frequency (without Casimir force)
xinit = np.interp(d,init_freq[0,:],init_freq[2,:])-d #Initial guess for offset value of top plate position due to Casimir force, in nm
#xinit = -2.89 #Initial guess for offset value of top plate position due to Casimir force
#Caspres = 2175.904 #Casimir pressure from BCS model in Pascal
n = -3.19273178 #Extracted from fit to Casimir_pressure_dist_v3.h5, 95% confidence interval is +-2.37283016e-03
Caspres = 6.95778690e7/(d**-n) #Fitted function to Casimir_pressure_dist_v3.h5, to cover any distance d (in nm).

#Fext is rescaled by 1e-12 to make the actual value closer to 1.
output = eng.Casimir_single(omega_m,meff,gamma_m,n,d,Caspres,area,xinit,Fext*1e-12)
cascurve = np.array(output._data).reshape((2,int(output.size[1]/2)))

ampl_opto = simpowers_VNAsweep_20[4,12] #amplitude of the optomechanical motion, in nm
ampl_cas = np.max(cascurve[0,:])

#%%
#d = 19.0

#cascurve_6a = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[8,12])
#omega_m6a = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 18.5

#cascurve_6b = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[8,12])
#omega_m6b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_6c = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[8,12])
omega_m6c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 0
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -11.0#16.5 #This is in log-scale
noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -10.9#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -10.7#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs3[:,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift
data1 = filtdata_up3[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down3[:,sbnum,4,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_6a[0,:]*np.sqrt(np.interp((cascurve_6a[0,:])*1e-9,xmotions/2,All_PF[1,:]))

#cascurve_f = (cascurve_6a[1,:]/2/np.pi*omega_m6a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_6b[0,:]*np.sqrt(np.interp((cascurve_6b[0,:])*1e-9,xmotions/2,All_PF[1,:]))

#cascurve_f2 = (cascurve_6b[1,:]/2/np.pi*omega_m6b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_6c[0,:]*np.sqrt(np.interp((cascurve_6c[0,:])*1e-9,xmotions/2,All_PF[1,:]))

cascurve_f3 = (cascurve_6c[1,:]/2/np.pi*omega_m6c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%
#d = 19.0

#cascurve_1a = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[3,12])
#omega_m1a = np.interp(d,init_freq[0,:],init_freq[1,:])

#d = 18.5

#cascurve_1b = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[3,12])
#omega_m1b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_1c = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[3,12])
omega_m1c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 4
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -11.9#16.5 #This is in log-scale
noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.6#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,3,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,3,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,3,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = 'o', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = 'o', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_1a[0,:]*np.interp((cascurve_1a[0,:])*1e-9,xmotions,All_PF[20,:])

#cascurve_f = (cascurve_1a[1,:]/2/np.pi*omega_m1a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_1b[0,:]*np.interp((cascurve_1b[0,:])*1e-9,xmotions,All_PF[20,:])

#cascurve_f2 = (cascurve_1b[1,:]/2/np.pi*omega_m1b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_1c[0,:]*np.interp((cascurve_1c[0,:])*1e-9,xmotions,All_PF[20,:])

cascurve_f3 = (cascurve_1c[1,:]/2/np.pi*omega_m1c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-1000,1000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-1000,0,1000])
ax2.set_xticklabels(['-1000','$+\omega_m$','+1000'])

#%%
#d = 19.0

#cascurve_2a = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[4,12])
#omega_m2a = np.interp(d,init_freq[0,:],init_freq[1,:])

d = 18.0

cascurve_2b = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[4,12])
omega_m2b = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 17.0

#cascurve_2c = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[4,12])
#omega_m2c = np.interp(d,init_freq[0,:],init_freq[1,:])
#%%
#d = 17.5

#cascurve_2d = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[4,12])
#omega_m2d = np.interp(d,init_freq[0,:],init_freq[1,:])

#d = 18.5

#cascurve_2e = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[4,12])
#omega_m2e = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 4
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -11.9#16.5 #This is in log-scale
noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*4.24 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.99 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.6#16.5 #This is in log-scale

hoffset4 = 2*np.pi*5.85 #in kHz, due to rescaling. (17 nm)
voffset4 = -11.6#16.5 #This is in log-scale

hoffset5 = 2*np.pi*0.76 #in kHz, due to rescaling. (17 nm)
voffset5 = -11.7#16.5 #This is in log-scale

fig2 = plt.figure(figsize = (6.4,3))
ax2 = fig2.add_subplot(111)

sbnum = 1 
freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Measured power (dBm)')

ax2.plot(freqs1,data1, marker = 'o', color = 'C2', label = 'Measured', alpha = 0.6)
ax2.plot(freqs1,data2, marker = 'o', color = 'C2', alpha = 0.6)

#comp_cascurve1 = cascurve_2a[0,:]*np.interp((cascurve_2a[0,:])*1e-9,xmotions,All_PF[76,:])

#cascurve_f = (cascurve_2a[1,:]/2/np.pi*omega_m2a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

comp_cascurve2 = cascurve_2b[0,:]*np.interp((cascurve_2b[0,:])*1e-9,xmotions,All_PF[76,:])

cascurve_f2 = (cascurve_2b[1,:]/2/np.pi*omega_m2b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

#comp_cascurve3 = cascurve_2c[0,:]*np.interp((cascurve_2c[0,:])*1e-9,xmotions,All_PF[76,:])

#cascurve_f3 = (cascurve_2c[1,:]/2/np.pi*omega_m2c-omega_1-hoffset3*1e3)/2/np.pi
#cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#comp_cascurve4 = cascurve_2d[0,:]*np.interp((cascurve_2d[0,:])*1e-9,xmotions,All_PF[76,:])

#cascurve_f4 = (cascurve_2d[1,:]/2/np.pi*omega_m2d-omega_1-hoffset4*1e3)/2/np.pi
#cascurve_a4 = 20*np.log10(np.sqrt(comp_cascurve4**2 + noisefl**2))+voffset4

#comp_cascurve5 = cascurve_2e[0,:]*np.interp((cascurve_2e[0,:])*1e-9,xmotions,All_PF[76,:])

#cascurve_f5 = (cascurve_2e[1,:]/2/np.pi*omega_m2e-omega_1-hoffset5*1e3)/2/np.pi
#cascurve_a5 = 20*np.log10(np.sqrt(comp_cascurve5**2 + noisefl**2))+voffset5

#ax2.plot(cascurve_f,cascurve_a, color = 'black', alpha = 0.2, label = '17 nm')

#ax2.plot(cascurve_f4,cascurve_a4, color = 'black', linestyle = 'dotted', alpha = 0.4, label = 'd = 17.5 nm')

ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid', label = 'd = 18 nm')

#ax2.plot(cascurve_f5,cascurve_a5, color = 'black', linestyle = 'dashed', alpha = 0.4, label = 'd = 18.5 nm')

#ax2.plot(cascurve_f3,cascurve_a3, color = 'black', alpha = 0.2, label = '19 nm')

ax2.set_ylim(-62,-35)
ax2.set_xlim(-1400,400)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-1000,-500,0])
ax2.set_xticklabels(['-1000','-500','$+\omega_m$'])

ax2.legend()

fig2.savefig('Casimir_fit_d_plot.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Casimir_fit_d_plot.svg', format = 'svg', bbox_inches = 'tight')

#%%
#d = 19.0

#cascurve_3a = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[0,12])
#omega_m3a = np.interp(d,init_freq[0,:],init_freq[1,:])

d = 18.0

cascurve_3b = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[0,12])
omega_m3b = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 17.0

#cascurve_3c = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[0,12])
#omega_m3c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 4
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.38 #in kHz, due to rescaling. (18.1 nm)
voffset = -11.9#16.5 #This is in log-scale
noisefl = 2.5e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*4.21 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.96 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.6#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_3a[0,:]*np.interp((cascurve_3a[0,:])*1e-9,xmotions,All_PF[1,:])

#cascurve_f = (cascurve_3a[1,:]/2/np.pi*omega_m3a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

comp_cascurve2 = cascurve_3b[0,:]*np.interp((cascurve_3b[0,:])*1e-9,xmotions,All_PF[1,:])

cascurve_f2 = (cascurve_3b[1,:]/2/np.pi*omega_m3b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

#comp_cascurve3 = cascurve_3c[0,:]*np.interp((cascurve_3c[0,:])*1e-9,xmotions,All_PF[1,:])

#cascurve_f3 = (cascurve_3c[1,:]/2/np.pi*omega_m3c-omega_1-hoffset3*1e3)/2/np.pi
#cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

#ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%
#d = 19.0

#cascurve_4a = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[1,12])
#omega_m4a = np.interp(d,init_freq[0,:],init_freq[1,:])

d = 18.0

cascurve_4b = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[1,12])
omega_m4b = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 17.0

#cascurve_4c = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[1,12])
#omega_m4c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 4
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.38 #in kHz, due to rescaling. (18.1 nm)
voffset = -11.9#16.5 #This is in log-scale
noisefl = 2.5e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*4.21 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.96 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.6#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,1,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,1,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,1,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_4a[0,:]*np.interp((cascurve_4a[0,:])*1e-9,xmotions,All_PF[7,:])

#cascurve_f = (cascurve_4a[1,:]/2/np.pi*omega_m4a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

comp_cascurve2 = cascurve_4b[0,:]*np.interp((cascurve_4b[0,:])*1e-9,xmotions,All_PF[7,:])

cascurve_f2 = (cascurve_4b[1,:]/2/np.pi*omega_m4b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

#comp_cascurve3 = cascurve_4c[0,:]*np.interp((cascurve_4c[0,:])*1e-9,xmotions,All_PF[7,:])

#cascurve_f3 = (cascurve_4c[1,:]/2/np.pi*omega_m4c-omega_1-hoffset3*1e3)/2/np.pi
#cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

#ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%
#d = 19.0

#cascurve_5a = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[2,12])
#omega_m5a = np.interp(d,init_freq[0,:],init_freq[1,:])

d = 18.0

cascurve_5b = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[2,12])
omega_m5b = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 17.0

#cascurve_5c = SimCas_w_prescribedamplitude(d,simpowers_VNAsweep_20[2,12])
#omega_m5c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 4
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.38 #in kHz, due to rescaling. (18.1 nm)
voffset = -11.9#16.5 #This is in log-scale
noisefl = 2.5e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*4.21 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.96 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.6#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,2,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,2,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_5a[0,:]*np.interp((cascurve_5a[0,:])*1e-9,xmotions,All_PF[13,:])

#cascurve_f = (cascurve_5a[1,:]/2/np.pi*omega_m5a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

comp_cascurve2 = cascurve_5b[0,:]*np.interp((cascurve_5b[0,:])*1e-9,xmotions,All_PF[13,:])

cascurve_f2 = (cascurve_5b[1,:]/2/np.pi*omega_m5b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

#comp_cascurve3 = cascurve_5c[0,:]*np.interp((cascurve_5c[0,:])*1e-9,xmotions,All_PF[13,:])

#cascurve_f3 = (cascurve_5c[1,:]/2/np.pi*omega_m5c-omega_1-hoffset3*1e3)/2/np.pi
#cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

#ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%

testamp = 10**(VNApower_sims/20)*1.94
#d = 19.0

#cascurve_7a = SimCas_w_prescribedamplitude(d,testamp[12])
#omega_m7a = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 18.5

#cascurve_7b = SimCas_w_prescribedamplitude(d,testamp[12])
#omega_m7b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_7c = SimCas_w_prescribedamplitude(d,testamp[12])
omega_m7c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 0
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -12.2#16.5 #This is in log-scale
noisefl = 2.2e-3 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -12.1#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.9#16.5 #This is in log-scale

fig2 = plt.figure(figsize = (3.2,3))
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs3[:,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift
data1 = filtdata_up3[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down3[:,sbnum,8,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Measured power (dBm)')

ax2.plot(freqs1,data1, marker = 'o', linestyle = 'None', color = 'C9', alpha = 0.6, label = 'Measured')
ax2.plot(freqs1,data2, marker = 'o', linestyle = 'None', color = 'C9', alpha = 0.6)
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_7a[0,:]*np.sqrt(np.interp((cascurve_7a[0,:])*1e-9,xmotions/2,All_PF[76,:]))

#cascurve_f = (cascurve_7a[1,:]/2/np.pi*omega_m7a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_7b[0,:]*np.sqrt(np.interp((cascurve_7b[0,:])*1e-9,xmotions/2,All_PF[76,:]))

#cascurve_f2 = (cascurve_7b[1,:]/2/np.pi*omega_m7b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_7c[0,:]*np.sqrt(np.interp((cascurve_7c[0,:])*1e-9,xmotions/2,All_PF[77,:]))

cascurve_f3 = (cascurve_7c[1,:]/2/np.pi*omega_m7c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

comp_cascurve4 = cascurve_7c[0,:]#*np.sqrt(np.interp((cascurve_7c[0,:])*1e-9,xmotions/2,All_PF[76,:]))

cascurve_f4 = (cascurve_7c[1,:]/2/np.pi*omega_m7c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a4 = 20*np.log10(np.sqrt(comp_cascurve4**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f4,cascurve_a4, color = 'black', linestyle = 'dashed', alpha = 0.7, label = 'No powerfactor')
ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'solid', label = 'Corrected')


ax2.set_ylim(-47,-26)
ax2.set_xlim(-10000,2000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-8000,-4000,0])
ax2.set_xticklabels(['-8000','-4000','$+\omega_m$'])

ax2.legend(loc = 'lower left', frameon = False)
fig2.savefig('Casimir_powerfactor_plot.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Casimir_powerfactor_plot.svg', format = 'svg', bbox_inches = 'tight')

#%%

testamp = 10**(VNApower_sims/20)*1.94
#d = 19.0

#cascurve_8a = SimCas_w_prescribedamplitude(d,testamp[11])
#omega_m8a = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 18.5

#cascurve_8b = SimCas_w_prescribedamplitude(d,testamp[11])
#omega_m8b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_8c = SimCas_w_prescribedamplitude(d,testamp[11])
omega_m8c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 0
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -12.2#16.5 #This is in log-scale
noisefl = 1.5e-3 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -12.1#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.9#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs3[:,sbnum,7,mwpower][1:-1]-fcenters3[0,sbnum]-fshift
data1 = filtdata_up3[:,sbnum,7,mwpower][1:-1]
data2 = filtdata_down3[:,sbnum,7,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_8a[0,:]*np.sqrt(np.interp((cascurve_8a[0,:])*1e-9,xmotions/2,All_PF[70,:]))

#cascurve_f = (cascurve_8a[1,:]/2/np.pi*omega_m8a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_8b[0,:]*np.sqrt(np.interp((cascurve_8b[0,:])*1e-9,xmotions/2,All_PF[70,:]))

#cascurve_f2 = (cascurve_8b[1,:]/2/np.pi*omega_m8b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_8c[0,:]*np.sqrt(np.interp((cascurve_8c[0,:])*1e-9,xmotions/2,All_PF[70,:]))

cascurve_f3 = (cascurve_8c[1,:]/2/np.pi*omega_m8c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-10000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%

testamp = 10**(VNApower_sims/20)*1.94
#d = 19.0

#cascurve_9a = SimCas_w_prescribedamplitude(d,testamp[10])
#omega_m9a = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 18.5

#cascurve_9b = SimCas_w_prescribedamplitude(d,testamp[10])
#omega_m9b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_9c = SimCas_w_prescribedamplitude(d,testamp[10])
omega_m9c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 0
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -12.2#16.5 #This is in log-scale
noisefl = 1.0e-3 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -12.1#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.9#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs3[:,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift
data1 = filtdata_up3[:,sbnum,6,mwpower][1:-1]
data2 = filtdata_down3[:,sbnum,6,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_9a[0,:]*np.sqrt(np.interp((cascurve_9a[0,:])*1e-9,xmotions/2,All_PF[63,:]))

#cascurve_f = (cascurve_9a[1,:]/2/np.pi*omega_m9a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_9b[0,:]*np.sqrt(np.interp((cascurve_9b[0,:])*1e-9,xmotions/2,All_PF[63,:]))

#cascurve_f2 = (cascurve_9b[1,:]/2/np.pi*omega_m9b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_9c[0,:]*np.sqrt(np.interp((cascurve_9c[0,:])*1e-9,xmotions/2,All_PF[63,:]))

cascurve_f3 = (cascurve_9c[1,:]/2/np.pi*omega_m9c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-10000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%
testamp = 10**(VNApower_sims/20)*1.94
#d = 19.0

#cascurve_10a = SimCas_w_prescribedamplitude(d,testamp[9])
#omega_m10a = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 18.5

#cascurve_10b = SimCas_w_prescribedamplitude(d,testamp[9])
#omega_m10b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_10c = SimCas_w_prescribedamplitude(d,testamp[9])
omega_m10c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 0
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -12.2#16.5 #This is in log-scale
noisefl = 1.0e-3 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -12.1#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.9#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs3[:,sbnum,5,mwpower][1:-1]-fcenters3[0,sbnum]-fshift
data1 = filtdata_up3[:,sbnum,5,mwpower][1:-1]
data2 = filtdata_down3[:,sbnum,5,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_10a[0,:]*np.sqrt(np.interp((cascurve_10a[0,:])*1e-9,xmotions/2,All_PF[57,:]))

#cascurve_f = (cascurve_10a[1,:]/2/np.pi*omega_m10a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_10b[0,:]*np.sqrt(np.interp((cascurve_10b[0,:])*1e-9,xmotions/2,All_PF[57,:]))

#cascurve_f2 = (cascurve_10b[1,:]/2/np.pi*omega_m10b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_10c[0,:]*np.sqrt(np.interp((cascurve_10c[0,:])*1e-9,xmotions/2,All_PF[57,:]))

cascurve_f3 = (cascurve_10c[1,:]/2/np.pi*omega_m10c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-10000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%
testamp = 10**(VNApower_sims/20)*1.94
#d = 19.0

#cascurve_11a = SimCas_w_prescribedamplitude(d,testamp[8])
#omega_m11a = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 18.5

#cascurve_11b = SimCas_w_prescribedamplitude(d,testamp[8])
#omega_m11b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_11c = SimCas_w_prescribedamplitude(d,testamp[8])
omega_m11c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 0
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -12.2#16.5 #This is in log-scale
noisefl = 7e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -12.1#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.9#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs3[:,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift
data1 = filtdata_up3[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down3[:,sbnum,4,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_11a[0,:]*np.sqrt(np.interp((cascurve_11a[0,:])*1e-9,xmotions/2,All_PF[51,:]))

#cascurve_f = (cascurve_11a[1,:]/2/np.pi*omega_m11a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_11b[0,:]*np.sqrt(np.interp((cascurve_11b[0,:])*1e-9,xmotions/2,All_PF[51,:]))

#cascurve_f2 = (cascurve_11b[1,:]/2/np.pi*omega_m11b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_11c[0,:]*np.sqrt(np.interp((cascurve_11c[0,:])*1e-9,xmotions/2,All_PF[51,:]))

cascurve_f3 = (cascurve_11c[1,:]/2/np.pi*omega_m11c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-10000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%
testamp = 10**(VNApower_sims/20)*1.94
#d = 19.0

#cascurve_12a = SimCas_w_prescribedamplitude(d,testamp[7])
#omega_m12a = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 18.5

#cascurve_12b = SimCas_w_prescribedamplitude(d,testamp[7])
#omega_m12b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_12c = SimCas_w_prescribedamplitude(d,testamp[7])
omega_m12c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 0
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -12.2#16.5 #This is in log-scale
noisefl = 5e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -12.1#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.9#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs3[:,sbnum,3,mwpower][1:-1]-fcenters3[0,sbnum]-fshift
data1 = filtdata_up3[:,sbnum,3,mwpower][1:-1]
data2 = filtdata_down3[:,sbnum,3,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_12a[0,:]*np.sqrt(np.interp((cascurve_12a[0,:])*1e-9,xmotions/2,All_PF[44,:]))

#cascurve_f = (cascurve_12a[1,:]/2/np.pi*omega_m12a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_12b[0,:]*np.sqrt(np.interp((cascurve_12b[0,:])*1e-9,xmotions/2,All_PF[44,:]))

#cascurve_f2 = (cascurve_12b[1,:]/2/np.pi*omega_m12b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_12c[0,:]*np.sqrt(np.interp((cascurve_12c[0,:])*1e-9,xmotions/2,All_PF[44,:]))

cascurve_f3 = (cascurve_12c[1,:]/2/np.pi*omega_m12c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-10000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%
testamp = 10**(VNApower_sims/20)*1.94
#d = 19.0

#cascurve_13a = SimCas_w_prescribedamplitude(d,testamp[6])
#omega_m13a = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 18.5

#cascurve_13b = SimCas_w_prescribedamplitude(d,testamp[6])
#omega_m13b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_13c = SimCas_w_prescribedamplitude(d,testamp[6])
omega_m13c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 0
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -12.2#16.5 #This is in log-scale
noisefl = 4e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -12.1#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.9#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs3[:,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift
data1 = filtdata_up3[:,sbnum,2,mwpower][1:-1]
data2 = filtdata_down3[:,sbnum,2,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_13a[0,:]*np.sqrt(np.interp((cascurve_13a[0,:])*1e-9,xmotions/2,All_PF[38,:]))

#cascurve_f = (cascurve_13a[1,:]/2/np.pi*omega_m13a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_13b[0,:]*np.sqrt(np.interp((cascurve_13b[0,:])*1e-9,xmotions/2,All_PF[38,:]))

#cascurve_f2 = (cascurve_13b[1,:]/2/np.pi*omega_m13b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_13c[0,:]*np.sqrt(np.interp((cascurve_13c[0,:])*1e-9,xmotions/2,All_PF[38,:]))

cascurve_f3 = (cascurve_13c[1,:]/2/np.pi*omega_m13c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-10000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])

#%%
testamp = 10**(VNApower_sims/20)*1.94
#d = 19.0

#cascurve_14a = SimCas_w_prescribedamplitude(d,testamp[5])
#omega_m14a = np.interp(d,init_freq[0,:],init_freq[1,:])


#d = 18.5

#cascurve_14b = SimCas_w_prescribedamplitude(d,testamp[5])
#omega_m14b = np.interp(d,init_freq[0,:],init_freq[1,:])


d = 18.0

cascurve_14c = SimCas_w_prescribedamplitude(d,testamp[5])
omega_m14c = np.interp(d,init_freq[0,:],init_freq[1,:])


#%%
#Plot the data to verify
mwpower = 0
sbnum = 1

#hoffset = 2*np.pi*4.22 #in kHz, due to rescaling. (18.1 nm)
hoffset = 2*np.pi*0.41 #in kHz, due to rescaling. (18.1 nm)
voffset = -12.2#16.5 #This is in log-scale
noisefl = 3e-4 #This is in linear scale
fshift = 500

hoffset2 = 2*np.pi*0.77 #in kHz, due to rescaling. (18 nm)
voffset2 = -12.1#16.5 #This is in log-scale

hoffset3 = 2*np.pi*4.24 #in kHz, due to rescaling. (17 nm)
voffset3 = -11.9#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs3[:,sbnum,1,mwpower][1:-1]-fcenters3[0,sbnum]-fshift
data1 = filtdata_up3[:,sbnum,1,mwpower][1:-1]
data2 = filtdata_down3[:,sbnum,1,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')
#ax2.plot(blue1_f-omega_1/2/np.pi,blue1_o + voffset, color = 'C2')
#ax2.plot(blue1_f2-omega_1/2/np.pi,blue1_o2 + voffset, color = 'C3')

#comp_cascurve1 = cascurve_14a[0,:]*np.sqrt(np.interp((cascurve_14a[0,:])*1e-9,xmotions/2,All_PF[32,:]))

#cascurve_f = (cascurve_14a[1,:]/2/np.pi*omega_m14a-omega_1-hoffset*1e3)/2/np.pi
#cascurve_a = 20*np.log10(np.sqrt(comp_cascurve1**2 + noisefl**2))+voffset

#comp_cascurve2 = cascurve_14b[0,:]*np.sqrt(np.interp((cascurve_14b[0,:])*1e-9,xmotions/2,All_PF[32,:]))

#cascurve_f2 = (cascurve_14b[1,:]/2/np.pi*omega_m14b-omega_1-hoffset2*1e3)/2/np.pi
#cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

comp_cascurve3 = cascurve_14c[0,:]*np.sqrt(np.interp((cascurve_14c[0,:])*1e-9,xmotions/2,All_PF[32,:]))

cascurve_f3 = (cascurve_14c[1,:]/2/np.pi*omega_m14c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

#ax2.plot(cascurve_f,cascurve_a, color = 'black', linestyle = 'dotted')

#ax2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

ax2.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'dashed')

ax2.set_ylim(-84,-20)
ax2.set_xlim(-10000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
#%%

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

ax3.set_xlabel('VNA power (dBm)')
ax3.set_ylabel('Amplitude (nm)')

testamp = 10**(VNApower_sims/20)*1.94
ax3.plot(VNApower_sims,simpowers_VNAsweep_20[:,12])
ax3.plot(VNApower_sims,testamp)

ax3.set_yscale('log')

#%%
mwpower = 4
sbnum = 1
fshift = 500

fig4 = plt.figure(figsize = (7,4))
ax4 = fig4.add_subplot(111)

ax4b = plt.twinx(ax4)

ax4.set_xlabel('Frequency (Hz)')
ax4.set_ylabel('Measured power (dBm)')

ax4b.set_ylabel('Displacement power (nm$^2$)')
ax4b.set_yscale('log')

hoffset2 = 2*np.pi*4.21 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale
noisefl = 2.5e-4

alphaval = 0.6
alphaval2 = 1.0

ROI1 = 80
ROI2 = 420

'''
#-37.5 dBm
mwpower = 4
ax4.plot(filtfreqs2[:,sbnum,1,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,1,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C1')
ax4.plot(filtfreqs2[:,sbnum,1,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,1,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C1')

comp_cascurve2 = cascurve_4b[0,:]*np.sqrt(np.interp((cascurve_4b[0,:])*1e-9,xmotions/2,All_PF[7,:]))
cascurve_f2 = (cascurve_4b[1,:]/2/np.pi*omega_m4b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = 20*np.log10(np.sqrt(comp_cascurve2**2 + noisefl**2))+voffset2

ax4.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid')

#-32.5 dBm
mwpower = 4
hoffset3 = 2*np.pi*4.21 
voffset3 = -11.8

ax4.plot(filtfreqs2[:,sbnum,3,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,3,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C3')
ax4.plot(filtfreqs2[:,sbnum,3,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,3,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C3')

comp_cascurve3 = cascurve_1c[0,:]*np.sqrt(np.interp((cascurve_1c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_1c[1,:]/2/np.pi*omega_m1c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

ax4.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'solid')

#-27.5 dBm
mwpower = 0
noisefl = 3e-4 #This is in linear scale

ax4.plot(filtfreqs3[:,sbnum,1,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,1,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C5')
ax4.plot(filtfreqs3[:,sbnum,1,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,1,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C5')

comp_cascurve3 = cascurve_14c[0,:]*np.sqrt(np.interp((cascurve_14c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_14c[1,:]/2/np.pi*omega_m14c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

ax4.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'solid')

#-22.5 dBm
mwpower = 0
noisefl = 5e-4 #This is in linear scale

ax4.plot(filtfreqs3[:,sbnum,3,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,3,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C7')
ax4.plot(filtfreqs3[:,sbnum,3,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,3,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C7')

comp_cascurve3 = cascurve_12c[0,:]*np.sqrt(np.interp((cascurve_12c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_12c[1,:]/2/np.pi*omega_m12c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

ax4.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'solid')

#-17.5 dBm
mwpower = 0
noisefl = 8e-4 #This is in linear scale

ax4.plot(filtfreqs3[:,sbnum,5,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,5,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C9')
ax4.plot(filtfreqs3[:,sbnum,5,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,5,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C9')

comp_cascurve3 = cascurve_10c[0,:]*np.sqrt(np.interp((cascurve_10c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_10c[1,:]/2/np.pi*omega_m10c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

ax4.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'solid')

#-12.5 dBm
mwpower = 0
noisefl = 1.5e-3 #This is in linear scale

ax4.plot(filtfreqs3[:,sbnum,7,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,7,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C1')
ax4.plot(filtfreqs3[:,sbnum,7,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,7,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C1')

comp_cascurve3 = cascurve_8c[0,:]*np.sqrt(np.interp((cascurve_8c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_8c[1,:]/2/np.pi*omega_m8c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = 20*np.log10(np.sqrt(comp_cascurve3**2 + noisefl**2))+voffset3

ax4.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'solid')
'''
#-10 dBm
mwpower = 0
noisefl = 2.2e-3

ax4.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', label = '-10 dBm', alpha = alphaval)
ax4.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', alpha = alphaval)

comp_cascurve3 = cascurve_7c[0,:]*np.sqrt(np.interp((cascurve_7c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_7c[1,:]/2/np.pi*omega_m7c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-15 dBm
mwpower = 0
noisefl = 1.0e-3 #This is in linear scale

ax4.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', label = '-15 dBm', alpha = alphaval)
ax4.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', alpha = alphaval)

comp_cascurve3 = cascurve_9c[0,:]*np.sqrt(np.interp((cascurve_9c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_9c[1,:]/2/np.pi*omega_m9c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-20 dBm
mwpower = 0
noisefl = 7e-4 #This is in linear scale

ax4.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', label = '-20 dBm', alpha = alphaval)
ax4.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', alpha = alphaval)

comp_cascurve3 = cascurve_11c[0,:]*np.sqrt(np.interp((cascurve_11c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_11c[1,:]/2/np.pi*omega_m11c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-25 dBm
mwpower = 0
noisefl = 4e-4 #This is in linear scale

ax4.plot(filtfreqs3[ROI1:ROI2,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,2,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C3', label = '-25 dBm', alpha = alphaval)
ax4.plot(filtfreqs3[ROI1:ROI2,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,2,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C3', alpha = alphaval)

comp_cascurve3 = cascurve_13c[0,:]*np.sqrt(np.interp((cascurve_13c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_13c[1,:]/2/np.pi*omega_m13c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-30 dBm
mwpower = 4
noisefl = 2.5e-4
ax4.plot(filtfreqs2[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C2', label = '-30 dBm', alpha = alphaval)
ax4.plot(filtfreqs2[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C2', alpha = alphaval)

comp_cascurve2 = cascurve_2b[0,:]*np.sqrt(np.interp((cascurve_2b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_2b[1,:]/2/np.pi*omega_m2b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f2[ROI3:flippedbit[0]],cascurve_a2[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f2[flippedbit[-1]:ROI4],cascurve_a2[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-35 dBm
mwpower = 4
noisefl = 2.5e-4
ax4.plot(filtfreqs2[ROI1:ROI2,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[ROI1:ROI2,sbnum,2,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C1', label = '-35 dBm', alpha = alphaval)
ax4.plot(filtfreqs2[ROI1:ROI2,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[ROI1:ROI2,sbnum,2,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C1', alpha = alphaval)

comp_cascurve2 = cascurve_5b[0,:]*np.sqrt(np.interp((cascurve_5b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_5b[1,:]/2/np.pi*omega_m5b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f2[ROI3:flippedbit[0]],cascurve_a2[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f2[flippedbit[-1]:ROI4],cascurve_a2[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-40 dBm
mwpower = 4
noisefl = 2.5e-4
ax4.plot(filtfreqs2[ROI1:ROI2,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[ROI1:ROI2,sbnum,0,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm', alpha = alphaval)
ax4.plot(filtfreqs2[ROI1:ROI2,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[ROI1:ROI2,sbnum,0,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C0', alpha = alphaval)

comp_cascurve2 = cascurve_3b[0,:]*np.sqrt(np.interp((cascurve_3b[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f2 = (cascurve_3b[1,:]/2/np.pi*omega_m3b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

ax4b.plot(cascurve_f2[ROI3:ROI4],cascurve_a2[ROI3:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2, label = 'Theory')

#Ylim matchin
ax4ylow = -79
ax4yhigh = -28

ax4scalel = (10**((ax4yhigh-ax4ylow)/10))
ax4disp_min = 1.9e-7
ax4b.set_ylim(ax4disp_min,ax4disp_min*ax4scalel)

fig4.legend(loc="lower right", bbox_to_anchor=(0.41,0.60), bbox_transform=ax4.transAxes, fontsize = 'small', frameon = False, ncol = 2, columnspacing = 1)

ax4.set_ylim(ax4ylow,ax4yhigh)
ax4.set_xlim(-9000,9500)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax4.set_xticks([-7500,-5000,-2500,0,2500,5000,7500])
ax4.set_xticklabels(['-7500','-5000','-2500','$+\omega_m$','+2500','+5000','+7500'])

fig4.savefig('Casimir_-20dBm_main.png', dpi = 300, bbox_inches = 'tight')
fig4.savefig('Casimir_-20dBm_main.svg', format = 'svg', bbox_inches = 'tight')

#%%
mwpower = 4
sbnum = 1
fshift = 500

fig4 = plt.figure(figsize = (2.5,3))
ax4 = fig4.add_subplot(111)

ax4b = plt.twinx(ax4)

#ax4.set_xlabel('Frequency (Hz)')
#ax4.set_ylabel('Measured power (dBm)')

#ax4b.set_ylabel('Displacement power (nm$^2$)')
ax4b.set_yscale('log')

hoffset2 = 2*np.pi*4.21 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale
noisefl = 2.5e-4
alphaval = 0.6
alphaval2 = 1.0

'''
#-10 dBm
mwpower = 0
noisefl = 2.2e-3

ax4.plot(filtfreqs3[:,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,8,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C6')
ax4.plot(filtfreqs3[:,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,8,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C6')

comp_cascurve3 = cascurve_7c[0,:]*np.sqrt(np.interp((cascurve_7c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_7c[1,:]/2/np.pi*omega_m7c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2*voffset2a

ax4b.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'solid')

#-15 dBm
mwpower = 0
noisefl = 1.0e-3 #This is in linear scale

ax4.plot(filtfreqs3[:,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,6,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C5')
ax4.plot(filtfreqs3[:,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,6,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C5')

comp_cascurve3 = cascurve_9c[0,:]*np.sqrt(np.interp((cascurve_9c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_9c[1,:]/2/np.pi*omega_m9c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2*voffset2a

ax4b.plot(cascurve_f3,cascurve_a3, color = 'black', linestyle = 'solid')

#-20 dBm
mwpower = 0
noisefl = 7e-4 #This is in linear scale


ax4.plot(filtfreqs3[:,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C4', alpha = alphaval)
ax4.plot(filtfreqs3[:,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C4', alpha = alphaval)

comp_cascurve3 = cascurve_11c[0,:]*np.sqrt(np.interp((cascurve_11c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_11c[1,:]/2/np.pi*omega_m11c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2*voffset2a

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f3[0:flippedbit[0]],cascurve_a3[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[-1]::],cascurve_a3[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)
'''
#-25 dBm
mwpower = 0
noisefl = 4e-4 #This is in linear scale

ax4.plot(filtfreqs3[:,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', alpha = alphaval)
ax4.plot(filtfreqs3[:,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', alpha = alphaval)

comp_cascurve3 = cascurve_13c[0,:]*np.sqrt(np.interp((cascurve_13c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_13c[1,:]/2/np.pi*omega_m13c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f3[0:flippedbit[0]],cascurve_a3[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f3[flippedbit[-1]::],cascurve_a3[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-30 dBm
mwpower = 4
noisefl = 2.5e-4
ax4.plot(filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', alpha = alphaval)
ax4.plot(filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', alpha = alphaval)

comp_cascurve2 = cascurve_2b[0,:]*np.sqrt(np.interp((cascurve_2b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_2b[1,:]/2/np.pi*omega_m2b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f2[0:flippedbit[0]],cascurve_a2[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f2[flippedbit[-1]::],cascurve_a2[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-35 dBm
mwpower = 4
noisefl = 2.5e-4
ax4.plot(filtfreqs2[:,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval)
ax4.plot(filtfreqs2[:,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval)

comp_cascurve2 = cascurve_5b[0,:]*np.sqrt(np.interp((cascurve_5b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_5b[1,:]/2/np.pi*omega_m5b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ax4b.plot(cascurve_f2[0:flippedbit[0]],cascurve_a2[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax4b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax4b.plot(cascurve_f2[flippedbit[-1]::],cascurve_a2[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-40 dBm
mwpower =4
noisefl = 2.5e-4
ax4.plot(filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval)
ax4.plot(filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval)

comp_cascurve2 = cascurve_3b[0,:]*np.sqrt(np.interp((cascurve_3b[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f2 = (cascurve_3b[1,:]/2/np.pi*omega_m3b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

ax4b.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid', alpha = alphaval2)

#Ylim matchin
ax4ylow = -60
ax4yhigh = -32

ax4scalel = (10**((ax4yhigh-ax4ylow)/10))
ax4disp_min = 1.5e-5
ax4b.set_ylim(ax4disp_min,ax4disp_min*ax4scalel)

ax4.set_ylim(ax4ylow,ax4yhigh)
ax4.set_xlim(-1800,500)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax4.set_xticks([-1500,-1000,-500,0,500])
ax4.set_xticklabels(['-1500','-1000','-500','$+\omega_m$','+500'])

fig4.savefig('Casimir_-20dBm_inset.png', dpi = 300, bbox_inches = 'tight')
fig4.savefig('Casimir_-20dBm_inset.svg', format = 'svg', bbox_inches = 'tight')

#%%
mwpower = 4
sbnum = 1
fshift = 500

fig5 = plt.figure(figsize = (7,4))

gs = gridspec.GridSpec(1, 2, width_ratios=[2,3], wspace = 0.0) 
ax5a = fig5.add_subplot(gs[0])
ax5b = fig5.add_subplot(gs[1])

ax5a2 = plt.twinx(ax5a)
ax5b2 = plt.twinx(ax5b)

ax5a.set_xlabel('Frequency (Hz)')
ax5a.set_ylabel('Measured power (dBm)')

ax5a2.set_yscale('log')

ax5b2.set_ylabel('Displacement power (nm$^2$)')
ax5b2.set_yscale('log')

hoffset2 = 2*np.pi*4.21 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale
noisefl = 2.5e-4
alphaval = 0.6
alphaval2 = 1.0

alphaval3 = 0.1
alphaval4 = 0.2

ROI1 = 80
ROI2 = 420

#-10 dBm
mwpower = 0
noisefl = 1e-3

ax5a.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', label = '-10 dBm', alpha = alphaval)
ax5a.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', alpha = alphaval)

ax5b.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', alpha = alphaval3)
ax5b.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', alpha = alphaval3)


comp_cascurve3 = cascurve_7c[0,:]*np.sqrt(np.interp((cascurve_7c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_7c[1,:]/2/np.pi*omega_m7c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax5a2.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax5a2.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax5a2.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax5b2.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval4)
ax5b2.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval4)
ax5b2.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval4)



#-15 dBm
mwpower = 0
noisefl = 5e-4 #This is in linear scale

ax5a.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', label = '-15 dBm', alpha = alphaval)
ax5a.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', alpha = alphaval)

ax5b.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', alpha = alphaval3)
ax5b.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', alpha = alphaval3)

comp_cascurve3 = cascurve_9c[0,:]*np.sqrt(np.interp((cascurve_9c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_9c[1,:]/2/np.pi*omega_m9c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax5a2.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax5a2.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax5a2.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax5b2.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval4)
ax5b2.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval4)
ax5b2.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval4)


#-20 dBm
mwpower = 0
noisefl = 4e-4 #This is in linear scale

ax5a.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', label = '-20 dBm', alpha = alphaval)
ax5a.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', alpha = alphaval)

ax5b.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', alpha = alphaval3)
ax5b.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', alpha = alphaval3)

comp_cascurve3 = cascurve_11c[0,:]*np.sqrt(np.interp((cascurve_11c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_11c[1,:]/2/np.pi*omega_m11c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax5a2.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax5a2.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax5a2.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax5b2.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval4)
ax5b2.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval4)
ax5b2.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval4)


#-25 dBm
mwpower = 0
noisefl = 2e-4 #This is in linear scale

ax5a.plot(filtfreqs3[:,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', label = '-25 dBm', alpha = alphaval)
ax5a.plot(filtfreqs3[:,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', alpha = alphaval)

ax5b.plot(filtfreqs3[:,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', alpha = alphaval)
ax5b.plot(filtfreqs3[:,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', alpha = alphaval)


comp_cascurve3 = cascurve_13c[0,:]*np.sqrt(np.interp((cascurve_13c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_13c[1,:]/2/np.pi*omega_m13c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax5a2.plot(cascurve_f3[0:flippedbit[0]],cascurve_a3[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax5a2.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax5a2.plot(cascurve_f3[flippedbit[-1]::],cascurve_a3[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax5b2.plot(cascurve_f3[0:flippedbit[0]],cascurve_a3[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax5b2.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax5b2.plot(cascurve_f3[flippedbit[-1]::],cascurve_a3[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)


#-30 dBm
mwpower = 4
noisefl = 2e-4
ax5a.plot(filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', label = '-30 dBm', alpha = alphaval)
ax5a.plot(filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', alpha = alphaval)

ax5b.plot(filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', alpha = alphaval)
ax5b.plot(filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', alpha = alphaval)


comp_cascurve2 = cascurve_2b[0,:]*np.sqrt(np.interp((cascurve_2b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_2b[1,:]/2/np.pi*omega_m2b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ax5a2.plot(cascurve_f2[0:flippedbit[0]],cascurve_a2[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax5a2.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax5a2.plot(cascurve_f2[flippedbit[-1]::],cascurve_a2[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax5b2.plot(cascurve_f2[0:flippedbit[0]],cascurve_a2[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax5b2.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax5b2.plot(cascurve_f2[flippedbit[-1]::],cascurve_a2[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)


#-35 dBm
mwpower = 4
noisefl = 2e-4
ax5a.plot(filtfreqs2[:,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval, label = '-35 dBm')
ax5a.plot(filtfreqs2[:,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval)

ax5b.plot(filtfreqs2[:,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval)
ax5b.plot(filtfreqs2[:,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval)


comp_cascurve2 = cascurve_5b[0,:]*np.sqrt(np.interp((cascurve_5b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_5b[1,:]/2/np.pi*omega_m5b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ax5a2.plot(cascurve_f2[0:flippedbit[0]],cascurve_a2[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax5a2.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax5a2.plot(cascurve_f2[flippedbit[-1]::],cascurve_a2[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax5b2.plot(cascurve_f2[0:flippedbit[0]],cascurve_a2[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax5b2.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax5b2.plot(cascurve_f2[flippedbit[-1]::],cascurve_a2[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-40 dBm
mwpower =4
noisefl = 2e-4
ax5a.plot(filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval, label = '-40 dBm')
ax5a.plot(filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval)

ax5b.plot(filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval)
ax5b.plot(filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval)


comp_cascurve2 = cascurve_3b[0,:]*np.sqrt(np.interp((cascurve_3b[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f2 = (cascurve_3b[1,:]/2/np.pi*omega_m3b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

ax5a2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid', alpha = alphaval2, label = 'Theory')
ax5b2.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid', alpha = alphaval2)

ax5a.set_xlim(-9000,-2000)
ax5b.set_xlim(-2000,800)

ax5a.spines['right'].set_visible(False)
ax5b.spines['left'].set_visible(False)
ax5a2.spines['right'].set_visible(False)
ax5b2.spines['left'].set_visible(False)

ax5a2.set_yticks([])
ax5a2.minorticks_off()
ax5b.set_yticks([])

ax5a.set_xticks([-8000,-6000,-4000,-2000])
ax5a.set_xticklabels(['-8000','-6000','-4000','-2000'])

ax5b.set_xticks([-1500,-1000,-500,0,500])
ax5b.set_xticklabels(['-1500','-1000','-500','$+\omega_m$','+500'])

d = .010  # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax5a.transAxes, color='k', clip_on=False)
ax5a.plot((1-d, 1+d), (-d, +d), **kwargs)
ax5a.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

kwargs.update(transform=ax5b.transAxes)  # switch to the bottom axes
ax5b.plot((-d, +d), (1-d, 1+d), **kwargs)
ax5b.plot((-d, +d), (-d, +d), **kwargs)

#Ylim matchin
ax5ylow = -76
ax5yhigh = -29

ax5scalel = (10**((ax5yhigh-ax5ylow)/10))
ax5disp_min = 3.8e-7
ax5a2.set_ylim(ax5disp_min,ax5disp_min*ax5scalel)
ax5b2.set_ylim(ax5disp_min,ax5disp_min*ax5scalel)

ax5a.set_ylim(ax5ylow,ax5yhigh)
ax5b.set_ylim(ax5ylow,ax5yhigh)

fig5.legend(loc="lower right", bbox_to_anchor=(1.02,-0.03), bbox_transform=ax5b.transAxes, frameon = False, ncol = 2, columnspacing = 1, fontsize = 'small')

fig5.savefig('Casimir_-20dBm_main_2.png', dpi = 300, bbox_inches = 'tight')
fig5.savefig('Casimir_-20dBm_main_2.svg', format = 'svg', bbox_inches = 'tight')

#%%

mwpower = 4
sbnum = 1
fshift = 500

fig6 = plt.figure(figsize = (6.4,3))
ax6a = fig6.add_subplot(111)

ax6b = plt.twinx(ax6a)

ax6a.set_xlabel('Frequency (Hz)')
ax6a.set_ylabel('Measured power (dBm)')

ax6b.set_ylabel('Displacement power (nm$^2$)')
ax6b.set_yscale('log')

hoffset2 = 2*np.pi*4.21 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale
noisefl = 2.5e-4
alphaval = 1.0
alphaval2 = 1.0

ROI1 = 80
ROI2 = 420

#-10 dBm
mwpower = 0
noisefl = 1e-3

ax6a.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', label = '-10 dBm', alpha = alphaval)
ax6a.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', alpha = alphaval)

comp_cascurve3 = cascurve_7c[0,:]*np.sqrt(np.interp((cascurve_7c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_7c[1,:]/2/np.pi*omega_m7c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax6b.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax6b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax6b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-15 dBm
mwpower = 0
noisefl = 5e-4 #This is in linear scale

ax6a.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', label = '-15 dBm', alpha = alphaval)
ax6a.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', alpha = alphaval)

comp_cascurve3 = cascurve_9c[0,:]*np.sqrt(np.interp((cascurve_9c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_9c[1,:]/2/np.pi*omega_m9c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax6b.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax6b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax6b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-20 dBm
mwpower = 0
noisefl = 4e-4 #This is in linear scale

ax6a.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', label = '-20 dBm', alpha = alphaval)
ax6a.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', alpha = alphaval)

comp_cascurve3 = cascurve_11c[0,:]*np.sqrt(np.interp((cascurve_11c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_11c[1,:]/2/np.pi*omega_m11c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 20000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax6b.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax6b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax6b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-25 dBm
mwpower = 0
noisefl = 2.5e-4 #This is in linear scale

ax6a.plot(filtfreqs3[:,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[:,sbnum,2,mwpower][1:-1], marker = '.', color = 'C3', label = '-25 dBm', alpha = alphaval)
ax6a.plot(filtfreqs3[:,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[:,sbnum,2,mwpower][1:-1], marker = '.', color = 'C3', alpha = alphaval)

comp_cascurve3 = cascurve_13c[0,:]*np.sqrt(np.interp((cascurve_13c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_13c[1,:]/2/np.pi*omega_m13c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ax6b.plot(cascurve_f3[0:flippedbit[0]],cascurve_a3[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax6b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax6b.plot(cascurve_f3[flippedbit[-1]::],cascurve_a3[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-30 dBm
mwpower = 4
noisefl = 2.5e-4

ax6a.plot(filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,4,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C2', label = '-30 dBm', alpha = alphaval)
ax6a.plot(filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,4,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C2', alpha = alphaval)

comp_cascurve2 = cascurve_2b[0,:]*np.sqrt(np.interp((cascurve_2b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_2b[1,:]/2/np.pi*omega_m2b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ax6b.plot(cascurve_f2[0:flippedbit[0]],cascurve_a2[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax6b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax6b.plot(cascurve_f2[flippedbit[-1]::],cascurve_a2[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-35 dBm
mwpower = 4
noisefl = 2.5e-4

ax6a.plot(filtfreqs2[:,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,2,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C1', alpha = alphaval, label = '-35 dBm')
ax6a.plot(filtfreqs2[:,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,2,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C1', alpha = alphaval)

comp_cascurve2 = cascurve_5b[0,:]*np.sqrt(np.interp((cascurve_5b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_5b[1,:]/2/np.pi*omega_m5b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ax6b.plot(cascurve_f2[0:flippedbit[0]],cascurve_a2[0:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax6b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax6b.plot(cascurve_f2[flippedbit[-1]::],cascurve_a2[flippedbit[-1]::], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-40 dBm
mwpower =4
noisefl = 2.5e-4

ax6a.plot(filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[:,sbnum,0,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C0', alpha = alphaval, label = '-40 dBm')
ax6a.plot(filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[:,sbnum,0,mwpower][1:-1], marker = '.', linestyle = 'None', color = 'C0', alpha = alphaval)

comp_cascurve2 = cascurve_3b[0,:]*np.sqrt(np.interp((cascurve_3b[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f2 = (cascurve_3b[1,:]/2/np.pi*omega_m3b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

ax6b.plot(cascurve_f2,cascurve_a2, color = 'black', linestyle = 'solid', alpha = alphaval2, label = 'Theory')

#Ylim matchin
ax6ylow = -84
ax6yhigh = -29

ax6scalel = (10**((ax6yhigh-ax6ylow)/10))
ax6disp_min = 6.0e-8
ax6b.set_ylim(ax6disp_min,ax6disp_min*ax6scalel)

ax6a.set_ylim(ax6ylow,ax6yhigh)

fig6.legend(loc="lower right", bbox_to_anchor=(0.48,0.60), bbox_transform=ax6a.transAxes, frameon = False, ncol = 2, columnspacing = 1, fontsize = 'small')
ax6a.set_xlim(-9000,5000)
ax6a.set_xticks([-8000,-6000,-4000,-2000,0,2000,4000])
ax6a.set_xticklabels(['-8000','-6000','-4000','-2000','$+\omega_m$','+2000','+4000'])

fig6.savefig('Casimir_-20dBm_main_3.png', dpi = 300, bbox_inches = 'tight')
fig6.savefig('Casimir_-20dBm_main_3.svg', format = 'svg', bbox_inches = 'tight')


#%%
mwpower = 4
sbnum = 1
fshift = 500

fig7 = plt.figure() #Main
fig8 = plt.figure() #Inset

ax7a = fig7.add_subplot(111)
ax8a = fig8.add_subplot(111)

ax7b = plt.twinx(ax7a)
ax8b = plt.twinx(ax8a)

ax7a.set_xlabel('Frequency (Hz)')
ax7a.set_ylabel('Measured power (dBm)')

ax7b.set_ylabel('Displacement power (nm$^2$)')
ax7b.set_yscale('log')

ax8a.set_xlabel('Frequency (Hz)')
ax8a.set_ylabel('Measured power (dBm)')

ax8b.set_ylabel('Displacement power (nm$^2$)')
ax8b.set_yscale('log')

hoffset2 = 2*np.pi*4.21 #in kHz, due to rescaling. (18 nm)
voffset2 = -11.8#16.5 #This is in log-scale
noisefl = 2.5e-4
alphaval = 0.6
alphaval2 = 1.0

alphaval3 = 0.1
alphaval4 = 0.2

ROI1 = 80
ROI2 = 420

ROI5 = 200
ROI6 = 280
#-10 dBm
mwpower = 0
noisefl = 1e-3

ax7a.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', label = '-10 dBm', alpha = alphaval)
ax7a.plot(filtfreqs3[ROI1:ROI2,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', alpha = alphaval)

ax8a.plot(filtfreqs3[ROI5:ROI6,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI5:ROI6,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', alpha = alphaval3, label = '-10 dBm')
ax8a.plot(filtfreqs3[ROI5:ROI6,sbnum,8,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI5:ROI6,sbnum,8,mwpower][1:-1], marker = '.', color = 'C6', alpha = alphaval3)


comp_cascurve3 = cascurve_7c[0,:]*np.sqrt(np.interp((cascurve_7c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_7c[1,:]/2/np.pi*omega_m7c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 12000)[0][0]
ROI4 = np.where(cascurve_f3 > -12000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

tst1 = flippedbit[np.where(cascurve_f3[flippedbit] < -12000)[0][-1]]
tst2 = np.where(cascurve_f3 < -12000)[0][0]

ROI7 = np.where(cascurve_f3 < 3000)[0][0]
ROI8 = np.where(cascurve_f3 > -3000)[0][-1]

tst7 = flippedbit[np.where(cascurve_f3[flippedbit] < -3000)[0][-1]]
tst8 = np.where(cascurve_f3 < -3000)[0][0]


ax7b.plot(cascurve_f3[ROI3:tst2],cascurve_a3[ROI3:tst2], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax7b.plot(cascurve_f3[tst1:flippedbit[-1]],cascurve_a3[tst1:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax7b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax8b.plot(cascurve_f3[ROI7:tst8],cascurve_a3[ROI7:tst8], color = 'black', linestyle = 'solid', alpha = alphaval4)
ax8b.plot(cascurve_f3[tst7:flippedbit[-1]],cascurve_a3[tst7:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval4)
ax8b.plot(cascurve_f3[flippedbit[-1]:ROI8],cascurve_a3[flippedbit[-1]:ROI8], color = 'black', linestyle = 'solid', alpha = alphaval4)

#-15 dBm
mwpower = 0
noisefl = 5e-4 #This is in linear scale

ax7a.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', label = '-15 dBm', alpha = alphaval)
ax7a.plot(filtfreqs3[ROI1:ROI2,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', alpha = alphaval)

ax8a.plot(filtfreqs3[ROI5:ROI6,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI5:ROI6,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', alpha = alphaval3, label = '-15 dBm')
ax8a.plot(filtfreqs3[ROI5:ROI6,sbnum,6,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI5:ROI6,sbnum,6,mwpower][1:-1], marker = '.', color = 'C5', alpha = alphaval3)

comp_cascurve3 = cascurve_9c[0,:]*np.sqrt(np.interp((cascurve_9c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_9c[1,:]/2/np.pi*omega_m9c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 12000)[0][0]
ROI4 = np.where(cascurve_f3 > -12000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

tst1 = flippedbit[np.where(cascurve_f3[flippedbit] < -12000)[0][-1]]
tst2 = np.where(cascurve_f3 < -12000)[0][0]

ROI7 = np.where(cascurve_f3 < 3000)[0][0]
ROI8 = np.where(cascurve_f3 > -3000)[0][-1]

tst7 = flippedbit[np.where(cascurve_f3[flippedbit] < -3000)[0][-1]]
tst8 = np.where(cascurve_f3 < -3000)[0][0]

ax7b.plot(cascurve_f3[ROI3:tst2],cascurve_a3[ROI3:tst2], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax7b.plot(cascurve_f3[tst1:flippedbit[-1]],cascurve_a3[tst1:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax7b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax8b.plot(cascurve_f3[ROI7:tst8],cascurve_a3[ROI7:tst8], color = 'black', linestyle = 'solid', alpha = alphaval4)
ax8b.plot(cascurve_f3[tst7:flippedbit[-1]],cascurve_a3[tst7:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval4)
ax8b.plot(cascurve_f3[flippedbit[-1]:ROI8],cascurve_a3[flippedbit[-1]:ROI8], color = 'black', linestyle = 'solid', alpha = alphaval4)

#-20 dBm
mwpower = 0
noisefl = 4e-4 #This is in linear scale

ax7a.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', label = '-20 dBm', alpha = alphaval)
ax7a.plot(filtfreqs3[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', alpha = alphaval)

ax8a.plot(filtfreqs3[ROI5:ROI6,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI5:ROI6,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', alpha = alphaval3, label = '-20 dBm')
ax8a.plot(filtfreqs3[ROI5:ROI6,sbnum,4,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI5:ROI6,sbnum,4,mwpower][1:-1], marker = '.', color = 'C4', alpha = alphaval3)

comp_cascurve3 = cascurve_11c[0,:]*np.sqrt(np.interp((cascurve_11c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_11c[1,:]/2/np.pi*omega_m11c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 12000)[0][0]
ROI4 = np.where(cascurve_f3 > -12000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

tst1 = flippedbit[np.where(cascurve_f3[flippedbit] < -11000)[0][-1]]
tst2 = np.where(cascurve_f3 < -12000)[0][0]

ROI7 = np.where(cascurve_f3 < 3000)[0][0]
ROI8 = np.where(cascurve_f3 > -3000)[0][-1]

tst7 = flippedbit[np.where(cascurve_f3[flippedbit] < -3000)[0][-1]]
tst8 = np.where(cascurve_f3 < -3000)[0][0]

ax7b.plot(cascurve_f3[ROI3:tst2],cascurve_a3[ROI3:tst2], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax7b.plot(cascurve_f3[tst1:flippedbit[-1]],cascurve_a3[tst1:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax7b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax8b.plot(cascurve_f3[ROI7:tst8],cascurve_a3[ROI7:tst8], color = 'black', linestyle = 'solid', alpha = alphaval4)
ax8b.plot(cascurve_f3[tst7:flippedbit[-1]],cascurve_a3[tst7:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval4)
ax8b.plot(cascurve_f3[flippedbit[-1]:ROI8],cascurve_a3[flippedbit[-1]:ROI8], color = 'black', linestyle = 'solid', alpha = alphaval4)

#-25 dBm
mwpower = 0
noisefl = 2e-4 #This is in linear scale

ax7a.plot(filtfreqs3[ROI1:ROI2,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI1:ROI2,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', label = '-25 dBm', alpha = alphaval)
ax7a.plot(filtfreqs3[ROI1:ROI2,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI1:ROI2,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', alpha = alphaval)

ax8a.plot(filtfreqs3[ROI5:ROI6,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_up3[ROI5:ROI6,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', alpha = alphaval, label = '-25 dBm')
ax8a.plot(filtfreqs3[ROI5:ROI6,sbnum,2,mwpower][1:-1]-fcenters3[0,sbnum]-fshift,filtdata_down3[ROI5:ROI6,sbnum,2,mwpower][1:-1], marker = 'o', color = 'C3', alpha = alphaval)


comp_cascurve3 = cascurve_13c[0,:]*np.sqrt(np.interp((cascurve_13c[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f3 = (cascurve_13c[1,:]/2/np.pi*omega_m13c-omega_1-hoffset3*1e3)/2/np.pi
cascurve_a3 = (np.sqrt(comp_cascurve3**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f3 < 12000)[0][0]
ROI4 = np.where(cascurve_f3 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f3))
flippedbit = np.where(freqsign>0)[0]

ROI7 = np.where(cascurve_f3 < 3000)[0][0]
ROI8 = np.where(cascurve_f3 > -3000)[0][-1]

tst7 = flippedbit[np.where(cascurve_f3[flippedbit] < -3000)[0][-1]]
tst8 = np.where(cascurve_f3 < -3000)[0][0]

ax7b.plot(cascurve_f3[ROI3:flippedbit[0]],cascurve_a3[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax7b.plot(cascurve_f3[flippedbit[0]:flippedbit[-1]],cascurve_a3[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax7b.plot(cascurve_f3[flippedbit[-1]:ROI4],cascurve_a3[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax8b.plot(cascurve_f3[ROI7:tst8],cascurve_a3[ROI7:tst8], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax8b.plot(cascurve_f3[tst7:flippedbit[-1]],cascurve_a3[tst7:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax8b.plot(cascurve_f3[flippedbit[-1]:ROI8],cascurve_a3[flippedbit[-1]:ROI8], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-30 dBm
mwpower = 4
noisefl = 2e-4
ax7a.plot(filtfreqs2[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', label = '-30 dBm', alpha = alphaval)
ax7a.plot(filtfreqs2[ROI1:ROI2,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[ROI1:ROI2,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', alpha = alphaval)

ax8a.plot(filtfreqs2[ROI5:ROI6,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[ROI5:ROI6,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', alpha = alphaval, label = '-30 dBm')
ax8a.plot(filtfreqs2[ROI5:ROI6,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[ROI5:ROI6,sbnum,4,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C2', alpha = alphaval)


comp_cascurve2 = cascurve_2b[0,:]*np.sqrt(np.interp((cascurve_2b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_2b[1,:]/2/np.pi*omega_m2b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f2 < 12000)[0][0]
ROI4 = np.where(cascurve_f2 > -20000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ROI7 = np.where(cascurve_f2 < 3000)[0][0]
ROI8 = np.where(cascurve_f2 > -3000)[0][-1]

ax7b.plot(cascurve_f2[ROI3:flippedbit[0]],cascurve_a2[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax7b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax7b.plot(cascurve_f2[flippedbit[-1]:ROI4],cascurve_a2[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax8b.plot(cascurve_f2[ROI7:flippedbit[0]],cascurve_a2[ROI7:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax8b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax8b.plot(cascurve_f2[flippedbit[-1]:ROI8],cascurve_a2[flippedbit[-1]:ROI8], color = 'black', linestyle = 'solid', alpha = alphaval2)


#-35 dBm
mwpower = 4
noisefl = 2.5e-4
ax7a.plot(filtfreqs2[ROI1:ROI2,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[ROI1:ROI2,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval, label = '-35 dBm')
ax7a.plot(filtfreqs2[ROI1:ROI2,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[ROI1:ROI2,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval)

ax8a.plot(filtfreqs2[ROI5:ROI6,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[ROI5:ROI6,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval, label = '-35 dBm')
ax8a.plot(filtfreqs2[ROI5:ROI6,sbnum,2,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[ROI5:ROI6,sbnum,2,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C1', alpha = alphaval)


comp_cascurve2 = cascurve_5b[0,:]*np.sqrt(np.interp((cascurve_5b[0,:])*1e-9,xmotions/2,All_PF[76,:]))
cascurve_f2 = (cascurve_5b[1,:]/2/np.pi*omega_m5b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f2 < 12000)[0][0]
ROI4 = np.where(cascurve_f2 > -25000)[0][-1]

freqsign = np.sign(np.diff(cascurve_f2))
flippedbit = np.where(freqsign>0)[0]

ROI7 = np.where(cascurve_f2 < 3000)[0][0]
ROI8 = np.where(cascurve_f2 > -4000)[0][-1]

ax7b.plot(cascurve_f2[ROI3:flippedbit[0]],cascurve_a2[ROI3:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax7b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax7b.plot(cascurve_f2[flippedbit[-1]:ROI4],cascurve_a2[flippedbit[-1]:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2)

ax8b.plot(cascurve_f2[ROI7:flippedbit[0]],cascurve_a2[ROI7:flippedbit[0]], color = 'black', linestyle = 'solid', alpha = alphaval2)
ax8b.plot(cascurve_f2[flippedbit[0]:flippedbit[-1]],cascurve_a2[flippedbit[0]:flippedbit[-1]], color = 'black', linestyle = 'dotted', alpha = alphaval2)
ax8b.plot(cascurve_f2[flippedbit[-1]:ROI8],cascurve_a2[flippedbit[-1]:ROI8], color = 'black', linestyle = 'solid', alpha = alphaval2)

#-40 dBm
mwpower =4
noisefl = 2.5e-4
ax7a.plot(filtfreqs2[ROI1:ROI2,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[ROI1:ROI2,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval, label = '-40 dBm')
ax7a.plot(filtfreqs2[ROI1:ROI2,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[ROI1:ROI2,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval)

ax8a.plot(filtfreqs2[ROI5:ROI6,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_up2[ROI5:ROI6,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval, label = '-40 dBm')
ax8a.plot(filtfreqs2[ROI5:ROI6,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift,filtdata_down2[ROI5:ROI6,sbnum,0,mwpower][1:-1], marker = 'o', linestyle = 'None', color = 'C0', alpha = alphaval)


comp_cascurve2 = cascurve_3b[0,:]*np.sqrt(np.interp((cascurve_3b[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f2 = (cascurve_3b[1,:]/2/np.pi*omega_m3b-omega_1-hoffset2*1e3)/2/np.pi
cascurve_a2 = (np.sqrt(comp_cascurve2**2 + noisefl**2))**2

ROI3 = np.where(cascurve_f2 < 15000)[0][0]
ROI4 = np.where(cascurve_f2 > -20000)[0][-1]

ROI7 = np.where(cascurve_f2 < 4000)[0][0]
ROI8 = np.where(cascurve_f2 > -5000)[0][-1]

ax7b.plot(cascurve_f2[ROI3:ROI4],cascurve_a2[ROI3:ROI4], color = 'black', linestyle = 'solid', alpha = alphaval2, label = 'Theory')
ax8b.plot(cascurve_f2[ROI7:ROI8],cascurve_a2[ROI7:ROI8], color = 'black', linestyle = 'solid', alpha = alphaval2, label = 'Theory')

#Ylim matchin
ax7ylow = -84#-76
ax7yhigh = -29

#Ylim matchin
ax8ylow = -76
ax8yhigh = -29

ax8scalel = (10**((ax5yhigh-ax5ylow)/10))
ax8disp_min = 3.8e-7

ax7scalel = (10**((ax7yhigh-ax7ylow)/10))
ax7disp_min = 5.7e-8
ax7b.set_ylim(ax7disp_min,ax7disp_min*ax7scalel)
ax8b.set_ylim(ax8disp_min,ax8disp_min*ax8scalel)

ax7a.set_ylim(ax7ylow,ax7yhigh)
ax8a.set_ylim(ax8ylow,ax8yhigh)

ax7a.set_xlim(-10000,5000)
ax7a.set_xticks([-10000,-8000,-6000,-4000,-2000,0,2000,4000])
ax7a.set_xticklabels(['-10000','-8000','-6000','-4000','-2000','$+\omega_m$','+2000','+4000'])

ax8a.set_xlim(-2000,1000)
ax8a.set_xticks([-2000,-1500,-1000,-500,0,500,1000])
ax8a.set_xticklabels(['-2000','-1500','-1000','-500','$+\omega_m$','+500','+1000'])

fig7.legend(frameon = False, ncol = 2, columnspacing = 1, fontsize = 'small')
fig8.legend(frameon = False, ncol = 2, columnspacing = 1, fontsize = 'small')

fig7.savefig('Casimir_-20dBm_main_v4.png', dpi = 300, bbox_inches = 'tight')
fig7.savefig('Casimir_-20dBm_main_v4.svg', format = 'svg', bbox_inches = 'tight')

fig8.savefig('Casimir_-20dBm_inset_v4.png', dpi = 300, bbox_inches = 'tight')
fig8.savefig('Casimir_-20dBm_inset_v4.svg', format = 'svg', bbox_inches = 'tight')

