# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 15:32:16 2025

@author: dejonm1
"""

import numpy as np
#import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.io
import matlab.engine

#%%
'''
First we calculate what the amplitude of the pressure should be such that the pressure takes the same value at a specific distance
I picked this distance to be 15.1 nm, as it is the same value as in the paper.
'''

ds = np.linspace(12,30,1001) #Distances, in nm

#The Casimir force exponent from theory
n1 = -3.19273178  
#Some different possible Casimir exponents
n2 = -2 
n3 = -3
n4 = -4

#Amplitude of the Casimir pressure from theory
Pcas_1 = -1275e-24
#Amplitudes of the Casimir pressure for other Casimir exponents such that they all yield the same force at 15.1 nm exactly.
Pcas_2 = -2.717e-12
Pcas_3 = -4.101e-20
Pcas_4 = -6.191e-28

Pressure_Casimir_1 = Pcas_1 * (ds*1e-9)**n1
Pressure_Casimir_2 = Pcas_2 * (ds*1e-9)**n2
Pressure_Casimir_3 = Pcas_3 * (ds*1e-9)**n3
Pressure_Casimir_4 = Pcas_4 * (ds*1e-9)**n4

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Plate separation (nm)')
ax1.set_ylabel('Pressure (Pa)')

ax1.plot(ds,-Pressure_Casimir_1, label = 'n = 3.193', color = 'black')
ax1.plot(ds,-Pressure_Casimir_2, label = 'n = 2') 
ax1.plot(ds,-Pressure_Casimir_3, label = 'n = 3')
ax1.plot(ds,-Pressure_Casimir_4, label = 'n = 4')

#ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.legend()

fig1.savefig('Plot_of_Casimir_pressure_various_n_equalized_at_15_1_nm.png', dpi = 300, bbox_inches = 'tight')
fig1.savefig('Plot_of_Casimir_pressure_various_n_equalized_at_15_1_nm.svg', format = 'svg', bbox_inches = 'tight')

#%%
#Next, there are some pre-simulations that we need to take care of.
#The nonlinear contribution from the Casimir force decreases the mechanical frequency (significantly). Since we cannot turn the Casimir force off, we only see the
#'dressed' frequency in our measurements, which is 10.001 MHz. The problem is that for each different value of the plate separation of the mechanical zero (d, in the paper,
#is the distance that the mechanical resonator would be if there were no Casimir force), the 'initial frequency' it takes to end up at the correct 'dressed frequency' is different.

def parabola(dist,k,offs,yoff):
    #This function calculates a parabola
    return k*(dist-offs)**2 + yoff

def fitparabola(dist,data,initialguess):
    #This function fits a parabola to some data, using scipy's curve_fit function.
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

#%%

ds = np.linspace(12,30,101) #Distances, in nm
xaxis = np.linspace(5,35,100001) #Distances at which the potential is evaluated.

Caspres_mechfreq = -1.275e-21*((xaxis*1e-9)**n1) #Casimir pressure for the mechanical frequency evaluation
Casforce = Caspres_mechfreq*np.pi*(0.5*11.3e-6)**2 #Casimir force based on the above and the drum parameters
Caspotential = np.cumsum(Casforce*(xaxis[-1] - xaxis[0])*1e-9/len(xaxis)) #Casimir potential based on the above.
Casplot = -Caspotential +Caspotential[-1] #Casimir potential when plotting

targetf = 10.01e6*np.ones(len(ds)) #Target frequency, which is the 'dressed frequency' that we measure.

omegas = 2*np.pi*40e6*np.ones(len(ds)) #Sets an initial guess for all the mechanical frequencies
#%%
#This cell needs to be run many times by hand, as I couldn't be bothered to write an optimization script for the frequencies.

#This line adjusts the mechanical frequencies (comment before running this cell for the first time, and uncomment after
#the script and paste the output from f_adjust instead of the array)


omegas = omegas - 2*np.pi*np.array([-1.24639618e+02,  5.05649248e+01, -2.28412623e+01, -5.69193390e+01,
        2.67652876e+01, -1.22413746e+02, -4.44480295e+01,  6.30214814e+01,
       -9.34760442e+01, -2.76683963e+02, -4.24547400e+01,  1.03765271e+02,
        1.31903597e+02,  1.48890708e+00,  2.39947696e+01,  4.89438760e+01,
        9.85061049e+01, -5.14460843e+01, -4.04939955e+01, -2.15781880e+01,
        1.94465064e-01, -1.45438126e-01, -2.48553370e+01, -6.26404012e-02,
        5.92124522e+01,  5.41620948e+01, -1.75536049e+01,  5.84824503e+00,
       -2.21861862e+01,  2.74893109e+01,  2.05738473e+01, -4.75334212e+00,
        3.15701486e+01,  3.02516525e-01,  1.07587630e+00,  8.34645471e-01,
       -1.61452201e-03,  2.79843509e-01,  3.01212547e+01,  1.15926327e-01,
       -9.47519773e+01, -8.85181515e+01, -8.93618337e+01,  7.76891647e-01,
        1.02466015e+00, -3.53875563e-01,  2.74754759e+00, -3.85652077e-01,
       -2.27441261e+01,  6.98838507e+00, -5.35179925e+00,  2.21501531e+00,
        9.85164165e-01,  1.67080303e+00,  4.20289623e+00,  6.83242305e-01,
        3.70596220e+00,  5.22442253e+00, -8.22521035e+01, -3.81527831e+00,
        3.55404953e+00,  4.20743894e+00,  1.88029656e+00, -4.19216136e+01,
       -1.98231574e+00, -5.28930353e+01,  5.12631321e+00,  3.18919950e+00,
        5.11744913e+00,  5.93218445e+00,  4.03906582e+01,  6.13496478e+00,
        1.68246376e+01,  4.42294746e+00,  2.73274935e+00,  4.89546802e+00,
        7.85485877e+00,  4.93902129e+00,  1.32956872e+01,  6.80192292e+00,
        6.52617428e+00,  6.85005099e+00,  6.12124719e+00,  6.45503977e+00,
        5.89121042e+00,  6.21424490e+00,  6.23116155e+00, -1.70277563e+00,
        6.11717118e+00,  5.89454978e+00,  5.91897427e+00,  5.90169354e+00,
        5.64463765e+00,  5.73141618e+00,  5.53589160e+00,  5.60563623e+00,
        5.33359486e+00,  5.38457387e+00,  5.54008546e+00,  5.45471866e+00,
        5.42864723e+00])

omegas = np.minimum(omegas,2*np.pi*5e7)

xstable = np.zeros(len(ds))
popts = np.zeros((3,len(ds)))
perrs = np.zeros((3,len(ds)))

#This fits a parabola to the calculate Casimir potential.
for idx,d in enumerate(ds):
    hbar = 1.0545718e-34
    kb = 1.38064852e-23
    omega_1 = omegas[idx]
    meff = 3.96e-14 #This is SI unit, which is kg. 
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
            
finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

#The 0.1*finalf/omegas is a scale factor, and you could replace it by a number like 0.1 or so.
f_adjust = (10*(finalf/omegas)**2)*(finalf - targetf)

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

ax3.set_xlabel('Initial separation (nm)')
ax3.set_ylabel('Frequency (Hz)')

finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

ax3.plot(ds,finalf, linestyle = 'None', marker = '.', label = 'Final frequency')
ax3.plot(ds,omegas/2/np.pi, linestyle = 'dotted', color = 'black', label = 'Start frequency')
ax3.plot(ds,targetf, linestyle = 'dashed', color = 'black', label = 'Target')
ax3.legend()
#ax3.set_ylim(0,4e7)

#Prints the output from f_adjust, for easy copying.
f_adjust

#%%

#After running this part of the script enough times, you should get the correct set of initial frequencies that result in the right 'dressed frequency'
omegas = np.array([2.29144409e+08, 2.22137534e+08, 2.15448043e+08, 2.09058071e+08,
       2.02949111e+08, 1.97107182e+08, 1.91516598e+08, 1.86163830e+08,
       1.81037941e+08, 1.76124706e+08, 1.71416249e+08, 1.66899523e+08,
       1.62566341e+08, 1.58409147e+08, 1.54418220e+08, 1.50585458e+08,
       1.46904585e+08, 1.43369252e+08, 1.39972489e+08, 1.36708132e+08,
       1.33571459e+08, 1.30556992e+08, 1.27660480e+08, 1.24876592e+08,
       1.22201217e+08, 1.19630992e+08, 1.17162203e+08, 1.14790336e+08,
       1.12513988e+08, 1.10327982e+08, 1.08230109e+08, 1.06218380e+08,
       1.04289599e+08, 1.02440963e+08, 1.00669734e+08, 9.89745745e+07,
       9.73525773e+07, 9.58009654e+07, 9.43178713e+07, 9.29017686e+07,
       9.15475789e+07, 9.02575182e+07, 8.90273328e+07, 8.78554507e+07,
       8.67387523e+07, 8.56753504e+07, 8.46627771e+07, 8.37007835e+07,
       8.27850575e+07, 8.19146834e+07, 8.10881817e+07, 8.03023605e+07,
       7.95560342e+07, 7.88470240e+07, 7.81732561e+07, 7.75342780e+07,
       7.69267128e+07, 7.63494193e+07, 7.58000680e+07, 7.52807852e+07,
       7.47858868e+07, 7.43153907e+07, 7.38679763e+07, 7.34420242e+07,
       7.30379692e+07, 7.26518290e+07, 7.22859335e+07, 7.19365100e+07,
       7.16039674e+07, 7.12868088e+07, 7.09837051e+07, 7.06957972e+07,
       7.04198190e+07, 7.01576694e+07, 6.99065873e+07, 6.96666307e+07,
       6.94368144e+07, 6.92178113e+07, 6.90078719e+07, 6.88068147e+07,
       6.86143554e+07, 6.84302199e+07, 6.82534350e+07, 6.80839450e+07,
       6.79215904e+07, 6.77655924e+07, 6.76158546e+07, 6.74721322e+07,
       6.73339943e+07, 6.72013188e+07, 6.70736231e+07, 6.69508626e+07,
       6.68328054e+07, 6.67190786e+07, 6.66096452e+07, 6.65041647e+07,
       6.64026079e+07, 6.63046659e+07, 6.62101826e+07, 6.61191431e+07,
       6.60312699e+07])

xstable = np.array([ 9.3053,  9.4559,  9.6074,  9.7601,  9.9134, 10.0679, 10.2236,
       10.3799, 10.538 , 10.697 , 10.8575, 11.0189, 11.1818, 11.3465,
       11.5124, 11.6798, 11.849 , 12.02  , 12.1925, 12.3665, 12.5426,
       12.7205, 12.9002, 13.082 , 13.2656, 13.4513, 13.6391, 13.829 ,
       14.0207, 14.2148, 14.4104, 14.6087, 14.8085, 15.0107, 15.2147,
       15.4205, 15.6281, 15.8375, 16.0484, 16.2611, 16.475 , 16.6904,
       16.9073, 17.1251, 17.3438, 17.5634, 17.7839, 18.0047, 18.2264,
       18.4481, 18.6701, 18.8924, 19.1147, 19.337 , 19.5593, 19.781 ,
       20.0027, 20.2241, 20.4446, 20.6651, 20.885 , 21.1043, 21.323 ,
       21.5411, 21.7583, 21.9752, 22.1912, 22.4066, 22.6211, 22.835 ,
       23.048 , 23.2607, 23.4725, 23.6834, 23.8937, 24.1034, 24.3125,
       24.5207, 24.7283, 24.9353, 25.1417, 25.3472, 25.5524, 25.757 ,
       25.9607, 26.1641, 26.3669, 26.5688, 26.7707, 26.9717, 27.1724,
       27.3725, 27.572 , 27.7712, 27.9698, 28.1681, 28.3658, 28.5632,
       28.7603, 28.9568, 29.153])

#Store the data
sumarr = np.vstack([ds,omegas,xstable])

np.savetxt("Optimized_mech_freqs_for_n=3_193_equalizedat15_1nm.txt",sumarr)

#%%

ds = np.linspace(12,30,101) #Distances, in nm
xaxis = np.linspace(5,35,100001) #Distances at which the potential is evaluated.

Caspres_mechfreq = -2.717e-12*((xaxis*1e-9)**n2) #Casimir pressure for the mechanical frequency evaluation
Casforce = Caspres_mechfreq*np.pi*(0.5*11.3e-6)**2 #Casimir force based on the above and the drum parameters
Caspotential = np.cumsum(Casforce*(xaxis[-1] - xaxis[0])*1e-9/len(xaxis)) #Casimir potential based on the above.
Casplot = -Caspotential +Caspotential[-1] #Casimir potential when plotting

targetf = 10.01e6*np.ones(len(ds)) #Target frequency, which is the 'dressed frequency' that we measure.

omegas = 2*np.pi*30e6*np.ones(len(ds)) #Sets an initial guess for all the mechanical frequencies
#%%
#This cell needs to be run many times by hand, as I couldn't be bothered to write an optimization script for the frequencies.

#This line adjusts the mechanical frequencies (comment before running this cell for the first time, and uncomment after
#the script and paste the output from f_adjust instead of the array)


omegas = omegas - 2*np.pi*np.array([-3.34177844e+01, -1.18490467e+01,  1.72314748e+00,  2.69068690e+01,
        8.09310962e+01,  6.97893116e+01,  3.56655015e+01,  1.44130548e+01,
        3.66282553e+01, -4.47884158e+01, -2.90986595e+01, -2.73967586e+01,
       -4.85900552e+01, -1.26301573e+01, -1.41946539e+01,  1.01325332e+01,
       -8.35308618e+00, -1.63683620e+01, -9.13719305e+01, -4.90711293e+01,
       -1.19299991e+01,  1.10032918e+01, -1.97615281e+01,  3.68598601e+01,
        9.82634039e+00,  2.24453809e+01, -5.93534639e+00,  1.52663051e+00,
       -1.14996166e+01,  3.60672647e+00,  2.43731671e+01,  2.34599116e+01,
        5.44532413e+00,  1.53751982e+01,  2.86346400e+01,  1.93219013e+00,
       -4.42990993e+01, -2.91395370e+01, -5.79323695e+01, -3.90408343e+01,
       -5.25856543e+01,  7.32336048e+01,  9.57291582e-01,  4.15818751e+01,
       -1.19835816e+00, -8.37699302e+00, -7.48870653e-01,  1.83210272e+01,
        8.33741824e+00,  9.04200616e+00,  5.11747622e+00, -1.86494578e+01,
        4.00299742e+00,  2.94784441e+00,  9.96523040e-01,  3.67348836e+00,
       -8.00183278e+00, -4.62847566e-01,  4.13785410e-01,  1.44290429e+01,
        3.25375437e+00,  3.16831424e+00,  3.95987887e+02, -5.57711068e+01,
       -9.08471551e+01, -5.94732313e+01,  5.48010787e-01,  5.71209243e-01,
        1.47265236e+02, -2.81455651e+01,  6.39219590e-01, -1.33576239e-01,
        2.51437035e+00,  8.26521252e+01,  3.60602316e+00,  4.82675113e+01,
        7.63955220e-01,  3.10981597e+00,  7.79374019e-01, -1.53472325e+01,
        8.48917457e-01,  7.14592504e+00, -9.62116434e+00, -1.13095830e+01,
        7.74784509e+00,  5.94056535e+00,  3.32404408e+00,  1.03300222e+00,
        1.01141655e+00,  2.93835925e+00,  3.09644708e+00,  3.05078684e+00,
        1.95639046e+00,  1.79306442e+00,  8.58040106e-01,  1.02671694e+00,
        8.93435306e-01,  9.05201623e-01,  1.06272624e+00,  1.72352539e+00,
        9.60759474e-01])

omegas = np.minimum(omegas,2*np.pi*5e7)

xstable = np.zeros(len(ds))
popts = np.zeros((3,len(ds)))
perrs = np.zeros((3,len(ds)))

#This fits a parabola to the calculate Casimir potential.
for idx,d in enumerate(ds):
    hbar = 1.0545718e-34
    kb = 1.38064852e-23
    omega_1 = omegas[idx]
    meff = 3.96e-14 #This is SI unit, which is kg. 
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
            
finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

#The 0.1*finalf/omegas is a scale factor, and you could replace it by a number like 0.1 or so.
f_adjust = (10*(finalf/omegas)**2)*(finalf - targetf)

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

ax3.set_xlabel('Initial separation (nm)')
ax3.set_ylabel('Frequency (Hz)')

finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

ax3.plot(ds,finalf, linestyle = 'None', marker = '.', label = 'Final frequency')
ax3.plot(ds,omegas/2/np.pi, linestyle = 'dotted', color = 'black', label = 'Start frequency')
ax3.plot(ds,targetf, linestyle = 'dashed', color = 'black', label = 'Target')
ax3.legend()
#ax3.set_ylim(0,4e7)

#Prints the output from f_adjust, for easy copying.
f_adjust

#%%

#After running this part of the script enough times, you should get the correct set of initial frequencies that result in the right 'dressed frequency'
omegas = np.array([1.64613034e+08, 1.61042537e+08, 1.57606442e+08, 1.54297649e+08,
       1.51111264e+08, 1.48041533e+08, 1.45082547e+08, 1.42230041e+08,
       1.39479138e+08, 1.36825568e+08, 1.34264936e+08, 1.31793842e+08,
       1.29408519e+08, 1.27105541e+08, 1.24882174e+08, 1.22734519e+08,
       1.20661029e+08, 1.18657810e+08, 1.16722931e+08, 1.14854880e+08,
       1.13050142e+08, 1.11306419e+08, 1.09623118e+08, 1.07997037e+08,
       1.06427593e+08, 1.04911462e+08, 1.03448075e+08, 1.02036334e+08,
       1.00673249e+08, 9.93582136e+07, 9.80893080e+07, 9.68658893e+07,
       9.56864343e+07, 9.45490449e+07, 9.34531974e+07, 9.23976045e+07,
       9.13797725e+07, 9.04006947e+07, 8.94572411e+07, 8.85493096e+07,
       8.76758091e+07, 8.68352494e+07, 8.60277138e+07, 8.52494609e+07,
       8.45026759e+07, 8.37839454e+07, 8.30929068e+07, 8.24287615e+07,
       8.17902830e+07, 8.11776100e+07, 8.05883332e+07, 8.00219693e+07,
       7.94779398e+07, 7.89554098e+07, 7.84531277e+07, 7.79701392e+07,
       7.75065554e+07, 7.70609472e+07, 7.66324933e+07, 7.62208654e+07,
       7.58243911e+07, 7.54437333e+07, 7.50746991e+07, 7.47243080e+07,
       7.43847788e+07, 7.40591932e+07, 7.37465490e+07, 7.34442619e+07,
       7.31517165e+07, 7.28725180e+07, 7.26029686e+07, 7.23425158e+07,
       7.20918163e+07, 7.18494871e+07, 7.16166304e+07, 7.13910635e+07,
       7.11745638e+07, 7.09649477e+07, 7.07626474e+07, 7.05668680e+07,
       7.03783207e+07, 7.01951472e+07, 7.00191002e+07, 6.98486915e+07,
       6.96835593e+07, 6.95240095e+07, 6.93695108e+07, 6.92201082e+07,
       6.90754232e+07, 6.89352684e+07, 6.87995937e+07, 6.86680573e+07,
       6.85404860e+07, 6.84167726e+07, 6.82969138e+07, 6.81806272e+07,
       6.80677536e+07, 6.79581782e+07, 6.78519377e+07, 6.77486087e+07,
       6.76483513e+07])

xstable = np.array([ 8.4092,  8.555 ,  8.702 ,  8.8502,  8.9996,  9.1505,  9.3026,
        9.4565,  9.6116,  9.7682,  9.926 , 10.0856, 10.247 , 10.4096,
       10.574 , 10.7402, 10.9079, 11.0774, 11.2487, 11.4218, 11.5964,
       11.7728, 11.9513, 12.1313, 12.3134, 12.497 , 12.6827, 12.8699,
       13.0589, 13.2497, 13.442 , 13.6361, 13.832 , 14.0294, 14.228 ,
       14.4284, 14.63  , 14.8334, 15.0377, 15.2432, 15.4499, 15.6578,
       15.8669, 16.0766, 16.2875, 16.4993, 16.7114, 16.9244, 17.1377,
       17.3519, 17.5664, 17.7815, 17.9966, 18.212 , 18.4277, 18.6437,
       18.8594, 19.0754, 19.2914, 19.5074, 19.7234, 19.9391, 20.1545,
       20.3699, 20.585 , 20.8001, 21.0149, 21.2294, 21.4433, 21.6572,
       21.8708, 22.0841, 22.2968, 22.5092, 22.721 , 22.9328, 23.144 ,
       23.3549, 23.5652, 23.7752, 23.9849, 24.1943, 24.4031, 24.6113,
       24.8195, 25.0271, 25.2344, 25.4411, 25.6475, 25.8536, 26.0591,
       26.2643, 26.4692, 26.6738, 26.8778, 27.0815, 27.2849, 27.488 ,
       27.6905, 27.893 , 28.0949])

#Store the data
sumarr = np.vstack([ds,omegas,xstable])

np.savetxt("Optimized_mech_freqs_for_n=2_equalizedat15_1nm.txt",sumarr)

#%%
ds = np.linspace(12,30,101) #Distances, in nm
xaxis = np.linspace(5,35,100001) #Distances at which the potential is evaluated.

Caspres_mechfreq = -4.101e-20*((xaxis*1e-9)**n3) #Casimir pressure for the mechanical frequency evaluation
Casforce = Caspres_mechfreq*np.pi*(0.5*11.3e-6)**2 #Casimir force based on the above and the drum parameters
Caspotential = np.cumsum(Casforce*(xaxis[-1] - xaxis[0])*1e-9/len(xaxis)) #Casimir potential based on the above.
Casplot = -Caspotential +Caspotential[-1] #Casimir potential when plotting

targetf = 10.01e6*np.ones(len(ds)) #Target frequency, which is the 'dressed frequency' that we measure.

omegas = 2*np.pi*40e6*np.ones(len(ds)) #Sets an initial guess for all the mechanical frequencies
#%%
#This cell needs to be run many times by hand, as I couldn't be bothered to write an optimization script for the frequencies.

#This line adjusts the mechanical frequencies (comment before running this cell for the first time, and uncomment after
#the script and paste the output from f_adjust instead of the array)


omegas = omegas - 2*np.pi*np.array([2.60817006e-01, -1.20479541e+02,  4.58747855e+01,  2.21825778e+01,
        9.32761193e+01, -7.10079877e+01, -3.54241570e+01,  4.00092712e-01,
        4.65030533e+01, -6.87247621e+01,  2.67235196e+01,  5.60825063e+01,
       -9.84953474e+00,  7.32393658e+01,  5.02958183e+01,  7.50505699e+01,
       -2.62987175e+01,  1.29516715e+01,  7.89501980e-01, -1.51839100e+01,
       -2.41416202e+01,  4.82442934e+01,  1.13929675e-01,  1.83861619e-01,
       -3.20923124e+01, -2.73490641e+01,  8.04145810e+01, -2.18611626e+01,
       -2.87824822e+01, -7.22430695e+00, -2.83816894e+01, -3.06438693e+01,
        2.06160268e+00,  1.23483204e+00,  5.92816709e-01, -1.76550801e+01,
        3.07832428e+01, -7.98064890e+00,  3.81006091e+01, -5.03533246e+01,
       -2.81499378e+01, -7.98584893e+01, -1.51677483e+01, -3.53322398e+00,
       -5.42766477e+01,  4.73615703e+00,  2.95255290e+01,  2.17206274e+00,
        3.26458462e+01,  1.35621273e+01,  1.19231131e+01,  1.54721310e+01,
        8.52424596e+00,  1.88012002e+01,  1.79169538e+01,  2.17856579e+01,
        2.28408276e+01,  2.18897997e+01,  2.41712290e+01,  1.31901329e+01,
        2.14295497e+01,  2.69395530e+01, -1.80776224e+01,  2.43093883e+01,
        3.22775217e+01,  1.89436801e+01,  2.36772719e+01,  2.43479417e+01,
        3.63438131e+01,  1.76126098e+00,  3.05048020e+01,  3.43344488e+01,
        3.13109347e+01,  3.28694195e+01,  6.51208852e+01,  3.45184018e+01,
        2.92345184e+01,  1.83256661e+01,  2.16694083e+01,  2.80526032e+01,
        3.42505991e+01,  4.26495603e+01,  4.06502245e+01,  4.06823507e+01,
        4.09343149e+01,  4.04231804e+01,  3.90111347e+01,  4.09866963e+01,
        3.78309981e+01,  3.84368260e+01,  3.48846861e+01,  3.89076345e+01,
        3.88027984e+01,  3.31535805e+01,  3.59695037e+01,  3.86970094e+01,
        3.88314708e+01,  3.70802721e+01,  3.54951458e+01,  3.69059442e+01,
        3.79871650e+01])

omegas = np.minimum(omegas,2*np.pi*5e7)

xstable = np.zeros(len(ds))
popts = np.zeros((3,len(ds)))
perrs = np.zeros((3,len(ds)))

#This fits a parabola to the calculate Casimir potential.
for idx,d in enumerate(ds):
    hbar = 1.0545718e-34
    kb = 1.38064852e-23
    omega_1 = omegas[idx]
    meff = 3.96e-14 #This is SI unit, which is kg. 
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
            
finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

#This is a scale factor, and you could replace it by a number like 10 or so.
f_adjust = (10*(finalf/omegas)**2)*(finalf - targetf)

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

ax3.set_xlabel('Initial separation (nm)')
ax3.set_ylabel('Frequency (Hz)')

finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

ax3.plot(ds,finalf, linestyle = 'None', marker = '.', label = 'Final frequency')
ax3.plot(ds,omegas/2/np.pi, linestyle = 'dotted', color = 'black', label = 'Start frequency')
ax3.plot(ds,targetf, linestyle = 'dashed', color = 'black', label = 'Target')
ax3.legend()
#ax3.set_ylim(0,4e7)

#Prints the output from f_adjust, for easy copying.
f_adjust

#%%

#After running this part of the script enough times, you should get the correct set of initial frequencies that result in the right 'dressed frequency'
omegas = np.array([2.18217266e+08, 2.11854397e+08, 2.05770005e+08, 1.99950452e+08,
       1.94380325e+08, 1.89046787e+08, 1.83936721e+08, 1.79037530e+08,
       1.74340421e+08, 1.69834437e+08, 1.65509072e+08, 1.61356863e+08,
       1.57368901e+08, 1.53537996e+08, 1.49856639e+08, 1.46318546e+08,
       1.42917126e+08, 1.39645701e+08, 1.36501031e+08, 1.33475548e+08,
       1.30565800e+08, 1.27766576e+08, 1.25074624e+08, 1.22484749e+08,
       1.19994543e+08, 1.17599630e+08, 1.15295888e+08, 1.13083080e+08,
       1.10955344e+08, 1.08911130e+08, 1.06948155e+08, 1.05062905e+08,
       1.03253132e+08, 1.01517201e+08, 9.98531124e+07, 9.82576273e+07,
       9.67287142e+07, 9.52654239e+07, 9.38642749e+07, 9.25238642e+07,
       9.12432847e+07, 9.00181463e+07, 8.88498356e+07, 8.77341570e+07,
       8.66683060e+07, 8.56532456e+07, 8.46842082e+07, 8.37617817e+07,
       8.28825633e+07, 8.20450821e+07, 8.12480794e+07, 8.04891244e+07,
       7.97667975e+07, 7.90788744e+07, 7.84246170e+07, 7.78016844e+07,
       7.72090013e+07, 7.66448999e+07, 7.61079219e+07, 7.55961893e+07,
       7.51097205e+07, 7.46459751e+07, 7.42038120e+07, 7.37833923e+07,
       7.33825325e+07, 7.29998664e+07, 7.26352652e+07, 7.22872733e+07,
       7.19553344e+07, 7.16377628e+07, 7.13352976e+07, 7.10461053e+07,
       7.07694525e+07, 7.05049303e+07, 7.02517685e+07, 7.00098469e+07,
       6.97780136e+07, 6.95558013e+07, 6.93431440e+07, 6.91393444e+07,
       6.89439021e+07, 6.87562851e+07, 6.85762233e+07, 6.84034024e+07,
       6.82374416e+07, 6.80779959e+07, 6.79247556e+07, 6.77773035e+07,
       6.76356268e+07, 6.74992246e+07, 6.73679026e+07, 6.72414785e+07,
       6.71196907e+07, 6.70022994e+07, 6.68892191e+07, 6.67801033e+07,
       6.66748926e+07, 6.65734286e+07, 6.64754419e+07, 6.63808470e+07,
       6.62894425e+07])

xstable = np.array([ 9.191 ,  9.341 ,  9.4916,  9.6437,  9.7964,  9.9506, 10.1057,
       10.2617, 10.4192, 10.5779, 10.7378, 10.8989, 11.0618, 11.2259,
       11.3918, 11.5589, 11.7278, 11.8985, 12.0707, 12.2447, 12.4205,
       12.5981, 12.7778, 12.9593, 13.1426, 13.3283, 13.5155, 13.7048,
       13.8962, 14.0897, 14.2853, 14.4824, 14.6819, 14.8832, 15.0863,
       15.2912, 15.4979, 15.7064, 15.9164, 16.1279, 16.3412, 16.5551,
       16.7708, 16.9874, 17.2049, 17.4236, 17.6429, 17.8628, 18.083 ,
       18.3041, 18.5252, 18.7466, 18.968 , 19.1897, 19.4111, 19.6325,
       19.8536, 20.0744, 20.2949, 20.5151, 20.7347, 20.954 , 21.1724,
       21.3908, 21.6083, 21.8252, 22.0415, 22.2572, 22.4723, 22.6865,
       22.9004, 23.1134, 23.3258, 23.5376, 23.7488, 23.9591, 24.1688,
       24.3779, 24.5864, 24.7943, 25.0016, 25.2083, 25.4144, 25.6199,
       25.8248, 26.0291, 26.2328, 26.4362, 26.6387, 26.8409, 27.0425,
       27.2438, 27.4445, 27.6446, 27.8444, 28.0439, 28.2428, 28.4411,
       28.6391, 28.8368, 29.0342])

#Store the data
sumarr = np.vstack([ds,omegas,xstable])

np.savetxt("Optimized_mech_freqs_for_n=3_equalizedat15_1nm.txt",sumarr)


#%%
ds = np.linspace(12,30,101) #Distances, in nm
xaxis = np.linspace(5,35,100001) #Distances at which the potential is evaluated.

Caspres_mechfreq = -6.191e-28*((xaxis*1e-9)**n4) #Casimir pressure for the mechanical frequency evaluation
Casforce = Caspres_mechfreq*np.pi*(0.5*11.3e-6)**2 #Casimir force based on the above and the drum parameters
Caspotential = np.cumsum(Casforce*(xaxis[-1] - xaxis[0])*1e-9/len(xaxis)) #Casimir potential based on the above.
Casplot = -Caspotential +Caspotential[-1] #Casimir potential when plotting

targetf = 10.01e6*np.ones(len(ds)) #Target frequency, which is the 'dressed frequency' that we measure.

omegas = 2*np.pi*60e6*np.ones(len(ds)) #Sets an initial guess for all the mechanical frequencies
#%%
#This cell needs to be run many times by hand, as I couldn't be bothered to write an optimization script for the frequencies.

#This line adjusts the mechanical frequencies (comment before running this cell for the first time, and uncomment after
#the script and paste the output from f_adjust instead of the array)


omegas = omegas - 2*np.pi*np.array([1.08309112e+00,  6.67695041e+01,  4.25898088e-02,  6.74011944e+01,
        6.55188037e+01, -1.29339650e+01,  7.10656939e+01,  4.89550612e+01,
        3.87680699e+00, -3.00235229e+01, -3.71475706e+01, -1.38369107e+02,
       -7.85764051e+00,  2.68767250e-01, -3.85943447e+01, -1.50814448e-02,
        7.34204753e+01,  2.34198383e+01, -1.25501483e+01,  2.84900629e+01,
       -1.13732476e+01,  1.02434718e+00, -2.30775102e+01,  2.70221634e+01,
        2.18211805e+01,  1.10617626e+01,  4.61955089e+01, -8.34043387e-02,
        2.92447444e+01, -1.97952420e+01, -2.74910791e+00, -4.87692336e-01,
        3.82885016e+00, -8.16859732e-01,  1.15722429e+01,  3.66689572e+01,
        9.59796180e+00,  1.66692653e+01,  1.65650871e+01, -1.04200608e-01,
        2.27421408e-01, -2.18793908e+01, -6.16181091e-02, -2.16744476e+01,
       -2.29879236e+00,  1.54981347e+00, -7.75187734e-01, -5.18249380e-01,
        2.07273967e+00,  1.45247906e+00,  4.80450249e+00,  1.75270153e+00,
        1.92288696e+00,  1.38851858e+00,  1.57234319e+00, -3.99294212e+01,
        1.60627565e+00,  1.60946373e+00,  1.88447586e+00, -1.74542229e+01,
        2.25025116e+00, -4.16492038e-01,  2.37022091e+00,  2.24175322e+00,
        2.23529393e+00,  1.46886780e+00,  1.87515691e+00,  7.26053252e+00,
        1.37148591e+00,  9.03680341e+00,  2.76107725e+00,  1.51675302e+00,
        4.83765942e+00,  3.86043660e+00,  2.39288658e+00,  1.99264812e+00,
        1.64971019e+00,  2.16120313e+00,  2.12473219e+00,  2.19093844e+00,
        1.76229854e+00,  2.16151764e+00,  2.12073385e+00,  2.08358704e+00,
        2.09078973e+00,  2.05841161e+00,  1.98691744e+00,  1.95352656e+00,
        1.96793546e+00,  1.81797183e+00,  1.50967318e+00,  1.84986109e+00,
        1.82280585e+00,  1.62650077e+00,  1.46427909e+00,  1.75499880e+00,
        1.59474994e+00,  1.66643137e+00,  1.68948742e+00,  1.69487497e+00,
        1.67143931e+00])

omegas = np.minimum(omegas,2*np.pi*5e7)

xstable = np.zeros(len(ds))
popts = np.zeros((3,len(ds)))
perrs = np.zeros((3,len(ds)))

#This fits a parabola to the calculate Casimir potential.
for idx,d in enumerate(ds):
    hbar = 1.0545718e-34
    kb = 1.38064852e-23
    omega_1 = omegas[idx]
    meff = 3.96e-14 #This is SI unit, which is kg. 
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
            
finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

#This is a scale factor, and you could replace it by a number like 10 or so.
f_adjust = (5*(finalf/omegas)**2)*(finalf - targetf)

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

ax3.set_xlabel('Initial separation (nm)')
ax3.set_ylabel('Frequency (Hz)')

finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

ax3.plot(ds,finalf, linestyle = 'None', marker = '.', label = 'Final frequency')
ax3.plot(ds,omegas/2/np.pi, linestyle = 'dotted', color = 'black', label = 'Start frequency')
ax3.plot(ds,targetf, linestyle = 'dashed', color = 'black', label = 'Target')
ax3.legend()
#ax3.set_ylim(0,4e7)

#Prints the output from f_adjust, for easy copying.
f_adjust

#%%

#After running this part of the script enough times, you should get the correct set of initial frequencies that result in the right 'dressed frequency'
omegas = np.array([2.77481592e+08, 2.67367093e+08, 2.57765501e+08, 2.48644826e+08,
       2.39974388e+08, 2.31727087e+08, 2.23876107e+08, 2.16399181e+08,
       2.09274618e+08, 2.02481107e+08, 1.96000552e+08, 1.89813893e+08,
       1.83909354e+08, 1.78267804e+08, 1.72877070e+08, 1.67723637e+08,
       1.62795146e+08, 1.58081626e+08, 1.53572453e+08, 1.49256242e+08,
       1.45126799e+08, 1.41173417e+08, 1.37389368e+08, 1.33767154e+08,
       1.30300210e+08, 1.26982307e+08, 1.23808170e+08, 1.20771190e+08,
       1.17866725e+08, 1.15090319e+08, 1.12437237e+08, 1.09904655e+08,
       1.07485159e+08, 1.05179738e+08, 1.02980286e+08, 1.00887043e+08,
       9.88934667e+07, 9.69987764e+07, 9.51985207e+07, 9.34901744e+07,
       9.18693900e+07, 9.03337054e+07, 8.88810420e+07, 8.75054528e+07,
       8.62062010e+07, 8.49790688e+07, 8.38219372e+07, 8.27296794e+07,
       8.16995371e+07, 8.07294901e+07, 7.98157100e+07, 7.89552485e+07,
       7.81452074e+07, 7.73831457e+07, 7.66650413e+07, 7.59876042e+07,
       7.53529661e+07, 7.47536177e+07, 7.41893005e+07, 7.36566638e+07,
       7.31559259e+07, 7.26832819e+07, 7.22370061e+07, 7.18161806e+07,
       7.14187266e+07, 7.10431125e+07, 7.06880515e+07, 7.03516712e+07,
       7.00347301e+07, 6.97336603e+07, 6.94483263e+07, 6.91787584e+07,
       6.89224435e+07, 6.86793250e+07, 6.84482077e+07, 6.82291113e+07,
       6.80206587e+07, 6.78222917e+07, 6.76336574e+07, 6.74539642e+07,
       6.72829425e+07, 6.71196843e+07, 6.69641860e+07, 6.68158145e+07,
       6.66741911e+07, 6.65389540e+07, 6.64097669e+07, 6.62863167e+07,
       6.61681298e+07, 6.60552756e+07, 6.59471442e+07, 6.58435433e+07,
       6.57443587e+07, 6.56493768e+07, 6.55582586e+07, 6.54708407e+07,
       6.53870371e+07, 6.53065779e+07, 6.52292962e+07, 6.51550469e+07,
       6.50838005e+07])

xstable = np.array([ 9.6998,  9.8531, 10.0073, 10.1621, 10.3178, 10.4747, 10.6319,
       10.7903, 10.9499, 11.1107, 11.2724, 11.435 , 11.5994, 11.765 ,
       11.9321, 12.1004, 12.2702, 12.4418, 12.6152, 12.7901, 12.9671,
       13.1459, 13.3265, 13.5095, 13.694 , 13.8812, 14.0702, 14.2616,
       14.4551, 14.651 , 14.8493, 15.0497, 15.2525, 15.4574, 15.6644,
       15.8738, 16.0853, 16.2983, 16.5134, 16.7303, 16.9487, 17.1683,
       17.3894, 17.6114, 17.8346, 18.0584, 18.2825, 18.5072, 18.7322,
       18.9572, 19.1819, 19.4069, 19.6313, 19.8551, 20.0786, 20.3012,
       20.5235, 20.7449, 20.9654, 21.185 , 21.404 , 21.6218, 21.839 ,
       22.055 , 22.2701, 22.4843, 22.6976, 22.91  , 23.1212, 23.3318,
       23.5415, 23.75  , 23.9579, 24.1649, 24.3713, 24.5765, 24.7811,
       24.9851, 25.1882, 25.3907, 25.5923, 25.7936, 25.994 , 26.1938,
       26.393 , 26.5916, 26.7896, 26.987 , 27.1841, 27.3803, 27.5762,
       27.7718, 27.9668, 28.1612, 28.3553, 28.5491, 28.7423, 28.9352,
       29.1278, 29.3201, 29.5118])

#Store the data
sumarr = np.vstack([ds,omegas,xstable])

np.savetxt("Optimized_mech_freqs_for_n=4_equalizedat15_1nm.txt",sumarr)

#%%
'''
Now that we have computed the initial frequencies for the different values of n, we load all the data

'''

init_freq_n1 = np.loadtxt("Optimized_mech_freqs_for_n=3_193_equalizedat15_1nm.txt")
init_freq_n2 = np.loadtxt("Optimized_mech_freqs_for_n=2_equalizedat15_1nm.txt")
init_freq_n3 = np.loadtxt("Optimized_mech_freqs_for_n=3_equalizedat15_1nm.txt")
init_freq_n4 = np.loadtxt("Optimized_mech_freqs_for_n=4_equalizedat15_1nm.txt")


#%%

#To save a lot of time, we load only two datasets and plot only the data measured for the Cavity drive power of -20 dBm
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

#%%
#We also rely on the optomechanical parameters not being changed, so we can import the calibration results/simulated displacement amplitudes
MWpower_sims = np.linspace(-40,-12.5,12)
VNApower_sims = np.linspace(-40,-10,13)

simpowers_VNAsweep_20 = np.loadtxt('Fit_Casimir_model_v27_Simdata_-20dBm_MW_altfreq.txt')

#We load the power factor outputs.
PF_file = np.loadtxt("Power_factors_9_12e-11.txt")
xmotions = PF_file[0,:]
lpower2s = PF_file[1,:]
All_PF = PF_file[2:,:]


#%%
#Initialize Matlab engine. (takes some time)
eng = matlab.engine.start_matlab()
eng.cd(r'C:\LocalUserData\User-data\dejonm1\MatCont\MatCont7p4', nargout = 0)

#%%
#The amplitudes predicted by our simulations are (in nm) stored in simpowers_VNAsweep_20[:,12]
#However we don't quite know what effective driving force would result in our Casimir simulation reaching this (maximum) amplitude
#The function below iterates the Casimir simulation for the given (hardcoded) Casimir parameters until the maximum amplitude has the right value.

def SimCas_w_prescribedamplitude(d,ampl_opto,Pcas,n,init_freqs,Fext = 0.1):
    #Fext = 0.1 #Initial value for Fext (driving force)
    
    meff = 3.96e-14 #In kg, from COMSOL simulation
    gamma_m = 2*np.pi*140#2*np.pi*170 #in Hz, Mechanical damping/loss rate (without Casimir force)
    area = np.pi*(0.5*11.3e-6)**2 #Area of the drum, from Louise's design file and verified by AFM
    omega_m = np.interp(d,init_freqs[0,:],init_freqs[1,:]) #in Hz, Mechanical frequency (without Casimir force)
    xinit = np.interp(d,init_freqs[0,:],init_freqs[2,:])-d #Initial guess for offset value of top plate position due to Casimir force, in nm
    #n1 = -3.193 #Extracted from fit to Casimir_pressure_dist_v3.h5, 95% confidence interval is +-2.37283016e-03
    Caspres = -Pcas/((d*1e-9)**-n) #Fitted function to Casimir_pressure_dist_v3.h5, to cover any distance d (in nm).
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
    
    return cascurve,Fext

#%%

d1 = 18.0
cascurve_30dBm_n1,Fext_30dBm_n1 = SimCas_w_prescribedamplitude(d1,1.94*10**(VNApower_sims[4]/20),Pcas_1,n1,init_freq_n1,Fext = 0.1363)
omega_m30dBm_n1 = np.interp(d1,init_freq_n1[0,:],init_freq_n1[1,:])
#%%
d2 = 16.2
cascurve_30dBm_n2,Fext_30dBm_n2 = SimCas_w_prescribedamplitude(d2,1.94*10**(VNApower_sims[4]/20),Pcas_2,n2,init_freq_n2,Fext = 0.13445)
omega_m30dBm_n2 = np.interp(d2,init_freq_n2[0,:],init_freq_n2[1,:])

#%%
d3 = 17.8
cascurve_30dBm_n3,Fext_30dBm_n3 = SimCas_w_prescribedamplitude(d3,1.94*10**(VNApower_sims[4]/20),Pcas_3,n3,init_freq_n3,Fext = 0.13446)
omega_m30dBm_n3 = np.interp(d3,init_freq_n3[0,:],init_freq_n3[1,:])

#%%
d4 = 18.5
cascurve_30dBm_n4,Fext_30dBm_n4 = SimCas_w_prescribedamplitude(d4,1.94*10**(VNApower_sims[4]/20),Pcas_4,n4,init_freq_n4,Fext = 0.1363)
omega_m30dBm_n4 = np.interp(d4,init_freq_n4[0,:],init_freq_n4[1,:])

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.31 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*3.56 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*2.3 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = 'Data')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n1 = cascurve_30dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_30dBm_n1[1,:]/2/np.pi*omega_m30dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

comp_cascurve3_n2 = cascurve_30dBm_n2[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2 = (cascurve_30dBm_n2[1,:]/2/np.pi*omega_m30dBm_n2-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2 = 20*np.log10(np.sqrt(comp_cascurve3_n2**2 + noisefl**2))+voffset_n2

comp_cascurve3_n3 = cascurve_30dBm_n3[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3 = (cascurve_30dBm_n3[1,:]/2/np.pi*omega_m30dBm_n3-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3 = 20*np.log10(np.sqrt(comp_cascurve3_n3**2 + noisefl**2))+voffset_n3

comp_cascurve3_n4 = cascurve_30dBm_n4[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4 = (cascurve_30dBm_n4[1,:]/2/np.pi*omega_m30dBm_n4-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4 = 20*np.log10(np.sqrt(comp_cascurve3_n4**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
ax2.legend()

#%%

#d1 = 18.1
cascurve_40dBm_n1,Fext_40dBm_n1 = SimCas_w_prescribedamplitude(d1,1.94*10**(VNApower_sims[0]/20),Pcas_1,n1,init_freq_n1,Fext = 0.0425)
omega_m40dBm_n1 = np.interp(d1,init_freq_n1[0,:],init_freq_n1[1,:])

#%%
#d2 = 16.1
cascurve_40dBm_n2,Fext_40dBm_n2 = SimCas_w_prescribedamplitude(d2,1.94*10**(VNApower_sims[0]/20),Pcas_2,n2,init_freq_n2,Fext = 0.0425)
omega_m40dBm_n2 = np.interp(d2,init_freq_n2[0,:],init_freq_n2[1,:])

#%%
#d3 = 17.8
cascurve_40dBm_n3,Fext_40dBm_n3 = SimCas_w_prescribedamplitude(d3,1.94*10**(VNApower_sims[0]/20),Pcas_3,n3,init_freq_n3,Fext = 0.0425)
omega_m40dBm_n3 = np.interp(d3,init_freq_n3[0,:],init_freq_n3[1,:])

#%%
#d4 = 18.5
cascurve_40dBm_n4,Fext_40dBm_n4 = SimCas_w_prescribedamplitude(d4,1.94*10**(VNApower_sims[0]/20),Pcas_4,n4,init_freq_n4,Fext = 0.0425)
omega_m40dBm_n4 = np.interp(d4,init_freq_n4[0,:],init_freq_n4[1,:])

#%%
#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.31 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*3.56 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*2.3 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = 'Data')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

comp_cascurve3_n1 = cascurve_40dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_40dBm_n1[1,:]/2/np.pi*omega_m40dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

comp_cascurve3_n2 = cascurve_40dBm_n2[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2 = (cascurve_40dBm_n2[1,:]/2/np.pi*omega_m40dBm_n2-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2 = 20*np.log10(np.sqrt(comp_cascurve3_n2**2 + noisefl**2))+voffset_n2

comp_cascurve3_n3 = cascurve_40dBm_n3[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3 = (cascurve_40dBm_n3[1,:]/2/np.pi*omega_m40dBm_n3-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3 = 20*np.log10(np.sqrt(comp_cascurve3_n3**2 + noisefl**2))+voffset_n3

comp_cascurve3_n4 = cascurve_40dBm_n4[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4 = (cascurve_40dBm_n4[1,:]/2/np.pi*omega_m40dBm_n4-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4 = 20*np.log10(np.sqrt(comp_cascurve3_n4**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
ax2.legend()

#%%

#d1 = 18.1
cascurve_20dBm_n1,Fext_20dBm_n1 = SimCas_w_prescribedamplitude(d1,1.94*10**(VNApower_sims[8]/20),Pcas_1,n1,init_freq_n1,Fext = 0.42475)
omega_m20dBm_n1 = np.interp(d1,init_freq_n1[0,:],init_freq_n1[1,:])

#%%
#d2 = 16.1
cascurve_20dBm_n2,Fext_20dBm_n2 = SimCas_w_prescribedamplitude(d2,1.94*10**(VNApower_sims[8]/20),Pcas_2,n2,init_freq_n2,Fext = 0.42475)
omega_m20dBm_n2 = np.interp(d2,init_freq_n2[0,:],init_freq_n2[1,:])

#%%
#d3 = 17.8
cascurve_20dBm_n3,Fext_20dBm_n3 = SimCas_w_prescribedamplitude(d3,1.94*10**(VNApower_sims[8]/20),Pcas_3,n3,init_freq_n3,Fext = 0.42475)
omega_m20dBm_n3 = np.interp(d3,init_freq_n3[0,:],init_freq_n3[1,:])

#%%
#d4 = 18.5
cascurve_20dBm_n4,Fext_20dBm_n4 = SimCas_w_prescribedamplitude(d4,1.94*10**(VNApower_sims[8]/20),Pcas_4,n4,init_freq_n4,Fext = 0.42475)
omega_m20dBm_n4 = np.interp(d4,init_freq_n4[0,:],init_freq_n4[1,:])

#%%
#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.31 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*3.56 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*2.3 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = 'Data')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

comp_cascurve3_n1 = cascurve_20dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_20dBm_n1[1,:]/2/np.pi*omega_m20dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

comp_cascurve3_n2 = cascurve_20dBm_n2[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2 = (cascurve_20dBm_n2[1,:]/2/np.pi*omega_m20dBm_n2-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2 = 20*np.log10(np.sqrt(comp_cascurve3_n2**2 + noisefl**2))+voffset_n2

comp_cascurve3_n3 = cascurve_20dBm_n3[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3 = (cascurve_20dBm_n3[1,:]/2/np.pi*omega_m20dBm_n3-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3 = 20*np.log10(np.sqrt(comp_cascurve3_n3**2 + noisefl**2))+voffset_n3

comp_cascurve3_n4 = cascurve_20dBm_n4[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4 = (cascurve_20dBm_n4[1,:]/2/np.pi*omega_m20dBm_n4-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4 = 20*np.log10(np.sqrt(comp_cascurve3_n4**2 + noisefl**2))+voffset_n4

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-12000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
ax2.legend()

#%%

'''
Compare all data plots and all fit plots for different n values

'''

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.31 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*3.56 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*2.3 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n1 = cascurve_20dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_20dBm_n1[1,:]/2/np.pi*omega_m20dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n1 = cascurve_30dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_30dBm_n1[1,:]/2/np.pi*omega_m30dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid')

comp_cascurve3_n1 = cascurve_40dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_40dBm_n1[1,:]/2/np.pi*omega_m40dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 3.193, d = 18.0 nm, P_cas = 1.275e-21')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_3_193_d_18_0.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_3_193_d_18_0.svg', format = 'svg', bbox_inches = 'tight')

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.31 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*3.56 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*2.3 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n2 = cascurve_20dBm_n2[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2 = (cascurve_20dBm_n2[1,:]/2/np.pi*omega_m20dBm_n2-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2 = 20*np.log10(np.sqrt(comp_cascurve3_n2**2 + noisefl**2))+voffset_n2

ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n2 = cascurve_30dBm_n2[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2 = (cascurve_30dBm_n2[1,:]/2/np.pi*omega_m30dBm_n2-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2 = 20*np.log10(np.sqrt(comp_cascurve3_n2**2 + noisefl**2))+voffset_n2

ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid')

comp_cascurve3_n2 = cascurve_40dBm_n2[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2 = (cascurve_40dBm_n2[1,:]/2/np.pi*omega_m40dBm_n2-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2 = 20*np.log10(np.sqrt(comp_cascurve3_n2**2 + noisefl**2))+voffset_n2

ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 2, d = 16.2 nm, P_cas = 2.712e-12')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_2_d_16_2.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_2_d_16_2.svg', format = 'svg', bbox_inches = 'tight')

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*1.78 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*3.35 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*3.56 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*2.3 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n3 = cascurve_20dBm_n3[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3 = (cascurve_20dBm_n3[1,:]/2/np.pi*omega_m20dBm_n3-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3 = 20*np.log10(np.sqrt(comp_cascurve3_n3**2 + noisefl**2))+voffset_n3

ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n3 = cascurve_30dBm_n3[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3 = (cascurve_30dBm_n3[1,:]/2/np.pi*omega_m30dBm_n3-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3 = 20*np.log10(np.sqrt(comp_cascurve3_n3**2 + noisefl**2))+voffset_n3

ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid')

comp_cascurve3_n3 = cascurve_40dBm_n3[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3 = (cascurve_40dBm_n3[1,:]/2/np.pi*omega_m40dBm_n3-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3 = 20*np.log10(np.sqrt(comp_cascurve3_n3**2 + noisefl**2))+voffset_n3

ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 3, d = 17.8 nm, P_cas = 4.101e-20')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_3_d_17_8.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_3_d_17_8.svg', format = 'svg', bbox_inches = 'tight')

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*1.78 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*3.35 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*3.56 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*2.3 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n4 = cascurve_20dBm_n4[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4 = (cascurve_20dBm_n4[1,:]/2/np.pi*omega_m20dBm_n4-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4 = 20*np.log10(np.sqrt(comp_cascurve3_n4**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n4 = cascurve_30dBm_n4[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4 = (cascurve_30dBm_n4[1,:]/2/np.pi*omega_m30dBm_n4-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4 = 20*np.log10(np.sqrt(comp_cascurve3_n4**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid')

comp_cascurve3_n4 = cascurve_40dBm_n4[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4 = (cascurve_40dBm_n4[1,:]/2/np.pi*omega_m40dBm_n4-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4 = 20*np.log10(np.sqrt(comp_cascurve3_n4**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 4, d = 18.5 nm, P_cas = 6.191e-28')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_4_d_18_5.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_4_d_18_5.svg', format = 'svg', bbox_inches = 'tight')

#%%

'''

And now see if I still can do this for 18.0 nm as the equal-pressure-point for all curves

'''


#Amplitudes of the Casimir pressure for other Casimir exponents such that they all yield the same force at 18.0 nm exactly.
Pcas_2a = -2.203e-12
Pcas_3a = -3.965e-20
Pcas_4a = -7.134e-28

Pressure_Casimir_1 = Pcas_1 * (ds*1e-9)**n1
Pressure_Casimir_2a = Pcas_2a * (ds*1e-9)**n2
Pressure_Casimir_3a = Pcas_3a * (ds*1e-9)**n3
Pressure_Casimir_4a = Pcas_4a * (ds*1e-9)**n4

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Plate separation (nm)')
ax1.set_ylabel('Pressure (Pa)')

ax1.plot(ds,-Pressure_Casimir_1, label = 'n = 3.193', color = 'black')
ax1.plot(ds,-Pressure_Casimir_2a, label = 'n = 2') 
ax1.plot(ds,-Pressure_Casimir_3a, label = 'n = 3')
ax1.plot(ds,-Pressure_Casimir_4a, label = 'n = 4')

#ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.legend()

fig1.savefig('Plot_of_Casimir_pressure_various_n_equalized_at_18_0_nm.png', dpi = 300, bbox_inches = 'tight')
fig1.savefig('Plot_of_Casimir_pressure_various_n_equalized_at_18_0_nm.svg', format = 'svg', bbox_inches = 'tight')

#%%

def Precalculate_mechfreqs(ds,omegas_init,Pcas,n):
    xaxis = np.linspace(5,35,100001) #Distances at which the potential is evaluated.
    Caspres_mechfreq = Pcas*((xaxis*1e-9)**n) #Casimir pressure for the mechanical frequency evaluation
    Casforce = Caspres_mechfreq*np.pi*(0.5*11.3e-6)**2 #Casimir force based on the above and the drum parameters
    Caspotential = np.cumsum(Casforce*(xaxis[-1] - xaxis[0])*1e-9/len(xaxis)) #Casimir potential based on the above.
    Casplot = -Caspotential +Caspotential[-1] #Casimir potential when plotting

    targetf = 10.01e6*np.ones(len(ds)) #Target frequency, which is the 'dressed frequency' that we measure.
    
    omegas_temp = omegas_init
    f_adjust = np.zeros(len(omegas_init))
    nmax = 201
    n = 0
    maxdiff = 10000
    
    while maxdiff > 3000 and n < nmax:
        omegas_temp = omegas_temp - 2*np.pi*f_adjust
        xstable = np.zeros(len(ds))
        popts = np.zeros((3,len(ds)))
        perrs = np.zeros((3,len(ds)))
    
        #This fits a parabola to the calculate Casimir potential.
        for idx,d in enumerate(ds):
            omega_1 = omegas_temp[idx]
            meff = 3.96e-14 #This is SI unit, which is kg. 
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
                    
        finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi
        
        #This is a scale factor, and you could replace it by a number like 10 or so.
        if n < 60:
            f_adjust = (5*(finalf/omegas_temp)**2)*(finalf - targetf)
        elif n < 120:
            f_adjust = (1*(finalf/omegas_temp)**2)*(finalf - targetf)
        elif n < 180:
            f_adjust = (0.1*(finalf/omegas_temp)**2)*(finalf - targetf)
        else:
            f_adjust = (0.01*(finalf/omegas_temp)**2)*(finalf - targetf)
        
        maxdiff = np.max(np.abs(finalf - targetf))
        
        n = n + 1
        '''
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)

        ax3.set_xlabel('Initial separation (nm)')
        ax3.set_ylabel('Frequency (Hz)')

        #finalf = np.sqrt(2*popts[0,:]*1e3/meff)/2/np.pi

        ax3.plot(ds,finalf, linestyle = 'None', marker = '.', label = 'Final frequency')
        #ax3.plot(ds,omegas/2/np.pi, linestyle = 'dotted', color = 'black', label = 'Start frequency')
        ax3.plot(ds,targetf, linestyle = 'dashed', color = 'black', label = 'Target')
        #ax3.legend()
        plt.show()
        plt.close(fig3)
        '''
        #print(n)
        print(maxdiff)
        
    omegas_final = omegas_temp - 2*np.pi*f_adjust
    xstable_final = xstable
    
    return omegas_final,xstable_final

ds = np.linspace(12,30,101) #Distances, in nm
omegas_final,xstable_final = Precalculate_mechfreqs(ds,np.ones(ds.shape[0])*40e6*2*np.pi,Pcas_2a,n2)
#Store the data
init_freq_n2a = np.vstack([ds,omegas_final,xstable_final])

#%%

omegas_final,xstable_final = Precalculate_mechfreqs(ds,np.ones(ds.shape[0])*40e6*2*np.pi,Pcas_3a,n3)
#Store the data
init_freq_n3a = np.vstack([ds,omegas_final,xstable_final])

#%%

omegas_final,xstable_final = Precalculate_mechfreqs(ds,np.ones(ds.shape[0])*60e6*2*np.pi,Pcas_4a,n4)
#Store the data
init_freq_n4a = np.vstack([ds,omegas_final,xstable_final])

#%%
d2a = 15.3
cascurve_30dBm_n2a,Fext_30dBm_n2a = SimCas_w_prescribedamplitude(d2a,1.94*10**(VNApower_sims[4]/20),Pcas_2a,n2,init_freq_n2a,Fext = 0.13445)
omega_m30dBm_n2a = np.interp(d2a,init_freq_n2a[0,:],init_freq_n2a[1,:])

#%%
d3a = 17.5
cascurve_30dBm_n3a,Fext_30dBm_n3a = SimCas_w_prescribedamplitude(d3a,1.94*10**(VNApower_sims[4]/20),Pcas_3a,n3,init_freq_n3a,Fext = 0.13446)
omega_m30dBm_n3a = np.interp(d3a,init_freq_n3a[0,:],init_freq_n3a[1,:])


d4a = 18.8
cascurve_30dBm_n4a,Fext_30dBm_n4a = SimCas_w_prescribedamplitude(d4a,1.94*10**(VNApower_sims[4]/20),Pcas_4a,n4,init_freq_n4a,Fext = 0.1363)
omega_m30dBm_n4a = np.interp(d4a,init_freq_n4a[0,:],init_freq_n4a[1,:])

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.82 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*6.05 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*4.9 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = 'Data')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n1 = cascurve_30dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_30dBm_n1[1,:]/2/np.pi*omega_m30dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

comp_cascurve3_n2a = cascurve_30dBm_n2a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2a = (cascurve_30dBm_n2a[1,:]/2/np.pi*omega_m30dBm_n2a-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2a = 20*np.log10(np.sqrt(comp_cascurve3_n2a**2 + noisefl**2))+voffset_n2

comp_cascurve3_n3a = cascurve_30dBm_n3a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3a = (cascurve_30dBm_n3a[1,:]/2/np.pi*omega_m30dBm_n3a-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3a = 20*np.log10(np.sqrt(comp_cascurve3_n3a**2 + noisefl**2))+voffset_n3

comp_cascurve3_n4a = cascurve_30dBm_n4a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4a = (cascurve_30dBm_n4a[1,:]/2/np.pi*omega_m30dBm_n4a-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4a = 20*np.log10(np.sqrt(comp_cascurve3_n4a**2 + noisefl**2))+voffset_n4

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2a,cascurve_a3_n2a, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3a,cascurve_a3_n3a, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4a,cascurve_a3_n4a, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
ax2.legend()

#%%
#d2a = 15.3
cascurve_20dBm_n2a,Fext_20dBm_n2a = SimCas_w_prescribedamplitude(d2a,1.94*10**(VNApower_sims[8]/20),Pcas_2a,n2,init_freq_n2a,Fext = 0.42475)
omega_m20dBm_n2a = np.interp(d2a,init_freq_n2a[0,:],init_freq_n2a[1,:])

#d3a = 17.5
cascurve_20dBm_n3a,Fext_20dBm_n3a = SimCas_w_prescribedamplitude(d3a,1.94*10**(VNApower_sims[8]/20),Pcas_3a,n3,init_freq_n3a,Fext = 0.42475)
omega_m20dBm_n3a = np.interp(d3a,init_freq_n3a[0,:],init_freq_n3a[1,:])

#d4a = 18.8
cascurve_20dBm_n4a,Fext_20dBm_n4a = SimCas_w_prescribedamplitude(d4a,1.94*10**(VNApower_sims[8]/20),Pcas_4a,n4,init_freq_n4a,Fext = 0.42475)
omega_m20dBm_n4a = np.interp(d4a,init_freq_n4a[0,:],init_freq_n4a[1,:])

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.82 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*6.05 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*4.9 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = 'Data')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

comp_cascurve3_n1 = cascurve_20dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_20dBm_n1[1,:]/2/np.pi*omega_m20dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

comp_cascurve3_n2a = cascurve_20dBm_n2a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2a = (cascurve_20dBm_n2a[1,:]/2/np.pi*omega_m20dBm_n2a-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2a = 20*np.log10(np.sqrt(comp_cascurve3_n2a**2 + noisefl**2))+voffset_n2

comp_cascurve3_n3a = cascurve_20dBm_n3a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3a = (cascurve_20dBm_n3a[1,:]/2/np.pi*omega_m20dBm_n3a-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3a = 20*np.log10(np.sqrt(comp_cascurve3_n3a**2 + noisefl**2))+voffset_n3

comp_cascurve3_n4a = cascurve_20dBm_n4a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4a = (cascurve_20dBm_n4a[1,:]/2/np.pi*omega_m20dBm_n4a-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4a = 20*np.log10(np.sqrt(comp_cascurve3_n4a**2 + noisefl**2))+voffset_n4

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2a,cascurve_a3_n2a, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3a,cascurve_a3_n3a, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4a,cascurve_a3_n4a, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
ax2.legend()

#%%
#d2a = 15.3
cascurve_40dBm_n2a,Fext_40dBm_n2a = SimCas_w_prescribedamplitude(d2a,1.94*10**(VNApower_sims[0]/20),Pcas_2a,n2,init_freq_n2a,Fext = 0.04)
omega_m40dBm_n2a = np.interp(d2a,init_freq_n2a[0,:],init_freq_n2a[1,:])

#d3a = 17.5
cascurve_40dBm_n3a,Fext_40dBm_n3a = SimCas_w_prescribedamplitude(d3a,1.94*10**(VNApower_sims[0]/20),Pcas_3a,n3,init_freq_n3a,Fext = 0.04)
omega_m40dBm_n3a = np.interp(d3a,init_freq_n3a[0,:],init_freq_n3a[1,:])

#d4a = 18.8
cascurve_40dBm_n4a,Fext_40dBm_n4a = SimCas_w_prescribedamplitude(d4a,1.94*10**(VNApower_sims[0]/20),Pcas_4a,n4,init_freq_n4a,Fext = 0.04)
omega_m40dBm_n4a = np.interp(d4a,init_freq_n4a[0,:],init_freq_n4a[1,:])

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.82 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*6.05 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*4.9 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = 'Data')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

comp_cascurve3_n1 = cascurve_40dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_40dBm_n1[1,:]/2/np.pi*omega_m40dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

comp_cascurve3_n2a = cascurve_40dBm_n2a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2a = (cascurve_40dBm_n2a[1,:]/2/np.pi*omega_m40dBm_n2a-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2a = 20*np.log10(np.sqrt(comp_cascurve3_n2a**2 + noisefl**2))+voffset_n2

comp_cascurve3_n3a = cascurve_40dBm_n3a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3a = (cascurve_40dBm_n3a[1,:]/2/np.pi*omega_m40dBm_n3a-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3a = 20*np.log10(np.sqrt(comp_cascurve3_n3a**2 + noisefl**2))+voffset_n3

comp_cascurve3_n4a = cascurve_40dBm_n4a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4a = (cascurve_40dBm_n4a[1,:]/2/np.pi*omega_m40dBm_n4a-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4a = 20*np.log10(np.sqrt(comp_cascurve3_n4a**2 + noisefl**2))+voffset_n4

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2a,cascurve_a3_n2a, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3a,cascurve_a3_n3a, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4a,cascurve_a3_n4a, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
ax2.legend()


#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.82 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*6.05 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*4.9 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n2a = cascurve_20dBm_n2a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2a = (cascurve_20dBm_n2a[1,:]/2/np.pi*omega_m20dBm_n2a-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2a = 20*np.log10(np.sqrt(comp_cascurve3_n2a**2 + noisefl**2))+voffset_n2

ax2.plot(cascurve_f3_n2a,cascurve_a3_n2a, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n2a = cascurve_30dBm_n2a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2a = (cascurve_30dBm_n2a[1,:]/2/np.pi*omega_m30dBm_n2a-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2a = 20*np.log10(np.sqrt(comp_cascurve3_n2a**2 + noisefl**2))+voffset_n2

ax2.plot(cascurve_f3_n2a,cascurve_a3_n2a, color = 'black', linestyle = 'solid')

comp_cascurve3_n2a = cascurve_40dBm_n2a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2a = (cascurve_40dBm_n2a[1,:]/2/np.pi*omega_m40dBm_n2a-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2a = 20*np.log10(np.sqrt(comp_cascurve3_n2a**2 + noisefl**2))+voffset_n2

ax2.plot(cascurve_f3_n2a,cascurve_a3_n2a, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 2, d = 15.3 nm, P_cas = 2.203e-12')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_2_d_15_3.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_2_d_15_3.svg', format = 'svg', bbox_inches = 'tight')

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.82 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*6.05 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*4.9 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n3a = cascurve_20dBm_n3a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3a = (cascurve_20dBm_n3a[1,:]/2/np.pi*omega_m20dBm_n3a-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3a = 20*np.log10(np.sqrt(comp_cascurve3_n3a**2 + noisefl**2))+voffset_n3

ax2.plot(cascurve_f3_n3a,cascurve_a3_n3a, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n3a = cascurve_30dBm_n3a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3a = (cascurve_30dBm_n3a[1,:]/2/np.pi*omega_m30dBm_n3a-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3a = 20*np.log10(np.sqrt(comp_cascurve3_n3a**2 + noisefl**2))+voffset_n3

ax2.plot(cascurve_f3_n3a,cascurve_a3_n3a, color = 'black', linestyle = 'solid')

comp_cascurve3_n3a = cascurve_40dBm_n3a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3a = (cascurve_40dBm_n3a[1,:]/2/np.pi*omega_m40dBm_n3a-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3a = 20*np.log10(np.sqrt(comp_cascurve3_n3a**2 + noisefl**2))+voffset_n3

ax2.plot(cascurve_f3_n3a,cascurve_a3_n3a, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 3, d = 17.5 nm, P_cas = 3.965e-20')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_3_d_17_5.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_3_d_17_5.svg', format = 'svg', bbox_inches = 'tight')

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*4.82 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*6.05 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*4.88 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n4a = cascurve_20dBm_n4a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4a = (cascurve_20dBm_n4a[1,:]/2/np.pi*omega_m20dBm_n4a-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4a = 20*np.log10(np.sqrt(comp_cascurve3_n4a**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n4a,cascurve_a3_n4a, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n4a = cascurve_30dBm_n4a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4a = (cascurve_30dBm_n4a[1,:]/2/np.pi*omega_m30dBm_n4a-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4a = 20*np.log10(np.sqrt(comp_cascurve3_n4a**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n4a,cascurve_a3_n4a, color = 'black', linestyle = 'solid')

comp_cascurve3_n4a = cascurve_40dBm_n4a[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4a = (cascurve_40dBm_n4a[1,:]/2/np.pi*omega_m40dBm_n4a-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4a = 20*np.log10(np.sqrt(comp_cascurve3_n4a**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n4a,cascurve_a3_n4a, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 4, d = 18.8 nm, P_cas = 7.134e-28')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_4_d_18_8.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_4_d_18_8.svg', format = 'svg', bbox_inches = 'tight')

#%%

'''

And now see if I still can do this when all the pressures are massively different (they still intersect at 15.1 nm nm)

'''


#Amplitudes of the Casimir pressure for other Casimir exponents such that they all yield the same force at 18.0 nm exactly.
Pcas_2b = -2.717e-13
Pcas_3b = -4.101e-21
Pcas_4b = -6.191e-29

Pressure_Casimir_1 = Pcas_1 * (ds*1e-9)**n1
Pressure_Casimir_2b = Pcas_2b * (ds*1e-9)**n2
Pressure_Casimir_3b = Pcas_3b * (ds*1e-9)**n3
Pressure_Casimir_4b = Pcas_4b * (ds*1e-9)**n4

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Plate separation (nm)')
ax1.set_ylabel('Pressure (Pa)')

ax1.plot(ds,-Pressure_Casimir_1, label = 'n = 3.193', color = 'black')
ax1.plot(ds,-Pressure_Casimir_2b, label = 'n = 2') 
ax1.plot(ds,-Pressure_Casimir_3b, label = 'n = 3')
ax1.plot(ds,-Pressure_Casimir_4b, label = 'n = 4')

#ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.legend()

fig1.savefig('Plot_of_Casimir_pressure_various_n_different_at_15_1_nm.png', dpi = 300, bbox_inches = 'tight')
fig1.savefig('Plot_of_Casimir_pressure_various_n_different_at_15_1_nm.svg', format = 'svg', bbox_inches = 'tight')

#%%

ds = np.linspace(8,30,101) #Distances, in nm
omegas_final,xstable_final = Precalculate_mechfreqs(ds,np.ones(ds.shape[0])*30e6*2*np.pi,Pcas_2b,n2)
#Store the data
init_freq_n2b = np.vstack([ds,omegas_final,xstable_final])



omegas_final,xstable_final = Precalculate_mechfreqs(ds,np.ones(ds.shape[0])*30e6*2*np.pi,Pcas_3b,n3)
#Store the data
init_freq_n3b = np.vstack([ds,omegas_final,xstable_final])



omegas_final,xstable_final = Precalculate_mechfreqs(ds,np.ones(ds.shape[0])*40e6*2*np.pi,Pcas_4b,n4)
#Store the data
init_freq_n4b = np.vstack([ds,omegas_final,xstable_final])
#%%

d2b = 8.8
cascurve_30dBm_n2b,Fext_30dBm_n2b = SimCas_w_prescribedamplitude(d2b,1.94*10**(VNApower_sims[4]/20),Pcas_2b,n2,init_freq_n2b,Fext = 0.13445)
omega_m30dBm_n2b = np.interp(d2b,init_freq_n2b[0,:],init_freq_n2b[1,:])

#%%
d3b = 10.9
cascurve_30dBm_n3b,Fext_30dBm_n3b = SimCas_w_prescribedamplitude(d3b,1.94*10**(VNApower_sims[4]/20),Pcas_3b,n3,init_freq_n3b,Fext = 0.13446)
omega_m30dBm_n3b = np.interp(d3b,init_freq_n3b[0,:],init_freq_n3b[1,:])


d4b = 12.4
cascurve_30dBm_n4b,Fext_30dBm_n4b = SimCas_w_prescribedamplitude(d4b,1.94*10**(VNApower_sims[4]/20),Pcas_4b,n4,init_freq_n4b,Fext = 0.1363)
omega_m30dBm_n4b = np.interp(d4b,init_freq_n4b[0,:],init_freq_n4b[1,:])

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*11.50 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*8.2 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*1.0 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = 'Data')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n1 = cascurve_30dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_30dBm_n1[1,:]/2/np.pi*omega_m30dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

comp_cascurve3_n2b = cascurve_30dBm_n2b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2b = (cascurve_30dBm_n2b[1,:]/2/np.pi*omega_m30dBm_n2b-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2b = 20*np.log10(np.sqrt(comp_cascurve3_n2b**2 + noisefl**2))+voffset_n2

comp_cascurve3_n3b = cascurve_30dBm_n3b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3b = (cascurve_30dBm_n3b[1,:]/2/np.pi*omega_m30dBm_n3b-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3b = 20*np.log10(np.sqrt(comp_cascurve3_n3b**2 + noisefl**2))+voffset_n3

comp_cascurve3_n4b = cascurve_30dBm_n4b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4b = (cascurve_30dBm_n4b[1,:]/2/np.pi*omega_m30dBm_n4b-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4b = 20*np.log10(np.sqrt(comp_cascurve3_n4b**2 + noisefl**2))+voffset_n4

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2b,cascurve_a3_n2b, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3b,cascurve_a3_n3b, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4b,cascurve_a3_n4b, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
ax2.legend()

#%%

#d2b = 8.8
cascurve_20dBm_n2b,Fext_20dBm_n2b = SimCas_w_prescribedamplitude(d2b,1.94*10**(VNApower_sims[8]/20),Pcas_2b,n2,init_freq_n2b,Fext = 0.42475)
omega_m20dBm_n2b = np.interp(d2b,init_freq_n2b[0,:],init_freq_n2b[1,:])

#d3b = 10.9
cascurve_20dBm_n3b,Fext_20dBm_n3b = SimCas_w_prescribedamplitude(d3b,1.94*10**(VNApower_sims[8]/20),Pcas_3b,n3,init_freq_n3b,Fext = 0.42475)
omega_m20dBm_n3b = np.interp(d3b,init_freq_n3b[0,:],init_freq_n3b[1,:])

#d4b = 12.4
cascurve_20dBm_n4b,Fext_20dBm_n4b = SimCas_w_prescribedamplitude(d4b,1.94*10**(VNApower_sims[8]/20),Pcas_4b,n4,init_freq_n4b,Fext = 0.42475)
omega_m20dBm_n4b = np.interp(d4b,init_freq_n4b[0,:],init_freq_n4b[1,:])

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*11.50 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*8.2 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*1.0 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = 'Data')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

comp_cascurve3_n1 = cascurve_20dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_20dBm_n1[1,:]/2/np.pi*omega_m20dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

comp_cascurve3_n2b = cascurve_20dBm_n2b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2b = (cascurve_20dBm_n2b[1,:]/2/np.pi*omega_m20dBm_n2b-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2b = 20*np.log10(np.sqrt(comp_cascurve3_n2b**2 + noisefl**2))+voffset_n2

comp_cascurve3_n3b = cascurve_20dBm_n3b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3b = (cascurve_20dBm_n3b[1,:]/2/np.pi*omega_m20dBm_n3b-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3b = 20*np.log10(np.sqrt(comp_cascurve3_n3b**2 + noisefl**2))+voffset_n3

comp_cascurve3_n4b = cascurve_20dBm_n4b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4b = (cascurve_20dBm_n4b[1,:]/2/np.pi*omega_m20dBm_n4b-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4b = 20*np.log10(np.sqrt(comp_cascurve3_n4b**2 + noisefl**2))+voffset_n4

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2b,cascurve_a3_n2b, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3b,cascurve_a3_n3b, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4b,cascurve_a3_n4b, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
ax2.legend()

#%%

#d2b = 8.8
cascurve_40dBm_n2b,Fext_40dBm_n2b = SimCas_w_prescribedamplitude(d2b,1.94*10**(VNApower_sims[0]/20),Pcas_2b,n2,init_freq_n2b,Fext = 0.04)
omega_m40dBm_n2b = np.interp(d2b,init_freq_n2b[0,:],init_freq_n2b[1,:])

#d3b = 10.9
cascurve_40dBm_n3b,Fext_40dBm_n3b = SimCas_w_prescribedamplitude(d3b,1.94*10**(VNApower_sims[0]/20),Pcas_3b,n3,init_freq_n3b,Fext = 0.04)
omega_m40dBm_n3b = np.interp(d3b,init_freq_n3b[0,:],init_freq_n3b[1,:])

#d4b = 12.4
cascurve_40dBm_n4b,Fext_40dBm_n4b = SimCas_w_prescribedamplitude(d4b,1.94*10**(VNApower_sims[0]/20),Pcas_4b,n4,init_freq_n4b,Fext = 0.04)
omega_m40dBm_n4b = np.interp(d4b,init_freq_n4b[0,:],init_freq_n4b[1,:])

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*11.50 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*8.2 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*1.0 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1
freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = 'Data')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

comp_cascurve3_n1 = cascurve_40dBm_n1[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n1 = (cascurve_40dBm_n1[1,:]/2/np.pi*omega_m40dBm_n1-omega_1-hoffset_n1*1e3)/2/np.pi
cascurve_a3_n1 = 20*np.log10(np.sqrt(comp_cascurve3_n1**2 + noisefl**2))+voffset_n1

comp_cascurve3_n2b = cascurve_40dBm_n2b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2b = (cascurve_40dBm_n2b[1,:]/2/np.pi*omega_m40dBm_n2b-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2b = 20*np.log10(np.sqrt(comp_cascurve3_n2b**2 + noisefl**2))+voffset_n2

comp_cascurve3_n3b = cascurve_40dBm_n3b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3b = (cascurve_40dBm_n3b[1,:]/2/np.pi*omega_m40dBm_n3b-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3b = 20*np.log10(np.sqrt(comp_cascurve3_n3b**2 + noisefl**2))+voffset_n3

comp_cascurve3_n4b = cascurve_40dBm_n4b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4b = (cascurve_40dBm_n4b[1,:]/2/np.pi*omega_m40dBm_n4b-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4b = 20*np.log10(np.sqrt(comp_cascurve3_n4b**2 + noisefl**2))+voffset_n4

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2b,cascurve_a3_n2b, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3b,cascurve_a3_n3b, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4b,cascurve_a3_n4b, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-5000,5000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,0,4000])
ax2.set_xticklabels(['-4000','$+\omega_m$','+4000'])
ax2.legend()

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*11.50 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*8.2 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*1.0 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n2b = cascurve_20dBm_n2b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2b = (cascurve_20dBm_n2b[1,:]/2/np.pi*omega_m20dBm_n2b-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2b = 20*np.log10(np.sqrt(comp_cascurve3_n2b**2 + noisefl**2))+voffset_n2

ax2.plot(cascurve_f3_n2b,cascurve_a3_n2b, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n2b = cascurve_30dBm_n2b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2b = (cascurve_30dBm_n2b[1,:]/2/np.pi*omega_m30dBm_n2b-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2b = 20*np.log10(np.sqrt(comp_cascurve3_n2b**2 + noisefl**2))+voffset_n2

ax2.plot(cascurve_f3_n2b,cascurve_a3_n2b, color = 'black', linestyle = 'solid')

comp_cascurve3_n2b = cascurve_40dBm_n2b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n2b = (cascurve_40dBm_n2b[1,:]/2/np.pi*omega_m40dBm_n2b-omega_1-hoffset_n2*1e3)/2/np.pi
cascurve_a3_n2b = 20*np.log10(np.sqrt(comp_cascurve3_n2b**2 + noisefl**2))+voffset_n2

ax2.plot(cascurve_f3_n2b,cascurve_a3_n2b, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 2, d = 8.8 nm, P_cas = 2.717e-13')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_2_d_8_8.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_2_d_8_8.svg', format = 'svg', bbox_inches = 'tight')

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*11.50 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*8.2 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*1.0 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n3b = cascurve_20dBm_n3b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3b = (cascurve_20dBm_n3b[1,:]/2/np.pi*omega_m20dBm_n3b-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3b = 20*np.log10(np.sqrt(comp_cascurve3_n3b**2 + noisefl**2))+voffset_n3

ax2.plot(cascurve_f3_n3b,cascurve_a3_n3b, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n3b = cascurve_30dBm_n3b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3b = (cascurve_30dBm_n3b[1,:]/2/np.pi*omega_m30dBm_n3b-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3b = 20*np.log10(np.sqrt(comp_cascurve3_n3b**2 + noisefl**2))+voffset_n3

ax2.plot(cascurve_f3_n3b,cascurve_a3_n3b, color = 'black', linestyle = 'solid')

comp_cascurve3_n3b = cascurve_40dBm_n3b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n3b = (cascurve_40dBm_n3b[1,:]/2/np.pi*omega_m40dBm_n3b-omega_1-hoffset_n3*1e3)/2/np.pi
cascurve_a3_n3b = 20*np.log10(np.sqrt(comp_cascurve3_n3b**2 + noisefl**2))+voffset_n3

ax2.plot(cascurve_f3_n3b,cascurve_a3_n3b, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 3, d = 10.9 nm, P_cas = 4.101e-21')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_3_d_10_9.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_3_d_10_9.svg', format = 'svg', bbox_inches = 'tight')

#%%

#Plot the data to verify
mwpower = 4
sbnum = 1

omega_1 = 2*np.pi*10.01e6

noisefl = 2e-4 #This is in linear scale
fshift = 500

hoffset_n1 = 2*np.pi*4.73 #in kHz, due to rescaling.
voffset_n1 = -11.3#16.5 #This is in log-scale

hoffset_n2 = 2*np.pi*11.50 #in kHz, due to rescaling.
voffset_n2 = -11.3#16.5 #This is in log-scale

hoffset_n3 = 2*np.pi*8.2 #in kHz, due to rescaling.
voffset_n3 = -11.3#16.5 #This is in log-scale

hoffset_n4 = 2*np.pi*1.0 #in kHz, due to rescaling.
voffset_n4 = -11.3#16.5 #This is in log-scale

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

sbnum = 1 #blue 1

ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Power (dBm)')

freqs1 = filtfreqs2[:,sbnum,8,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,8,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,8,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C2', label = '-20 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C2')

freqs1 = filtfreqs2[:,sbnum,4,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,4,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,4,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C1', label = '-30 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C1')

freqs1 = filtfreqs2[:,sbnum,0,mwpower][1:-1]-fcenters2[0,sbnum]-fshift
data1 = filtdata_up2[:,sbnum,0,mwpower][1:-1]
data2 = filtdata_down2[:,sbnum,0,mwpower][1:-1]

ax2.plot(freqs1,data1, marker = '.', linestyle = 'None', color = 'C0', label = '-40 dBm')
ax2.plot(freqs1,data2, marker = '.', linestyle = 'None', color = 'C0')

comp_cascurve3_n4b = cascurve_20dBm_n4b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4b = (cascurve_20dBm_n4b[1,:]/2/np.pi*omega_m20dBm_n4b-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4b = 20*np.log10(np.sqrt(comp_cascurve3_n4b**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n4b,cascurve_a3_n4b, color = 'black', linestyle = 'solid', label = 'Fit')

comp_cascurve3_n4b = cascurve_30dBm_n4b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4b = (cascurve_30dBm_n4b[1,:]/2/np.pi*omega_m30dBm_n4b-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4b = 20*np.log10(np.sqrt(comp_cascurve3_n4b**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n4b,cascurve_a3_n4b, color = 'black', linestyle = 'solid')

comp_cascurve3_n4b = cascurve_40dBm_n4b[0,:]#*np.sqrt(np.interp((cascurve_30dBm[0,:])*1e-9,xmotions/2,All_PF[1,:]))
cascurve_f3_n4b = (cascurve_40dBm_n4b[1,:]/2/np.pi*omega_m40dBm_n4b-omega_1-hoffset_n4*1e3)/2/np.pi
cascurve_a3_n4b = 20*np.log10(np.sqrt(comp_cascurve3_n4b**2 + noisefl**2))+voffset_n4

ax2.plot(cascurve_f3_n4b,cascurve_a3_n4b, color = 'black', linestyle = 'solid')

ax2.set_title('Fit n = 4, d = 12.4 nm, P_cas = 6.191e-29')

#ax2.plot(cascurve_f3_n1,cascurve_a3_n1, color = 'black', linestyle = 'solid', label = 'n = 3.193')
#ax2.plot(cascurve_f3_n2,cascurve_a3_n2, color = 'black', linestyle = 'solid', label = 'n = 2')
#ax2.plot(cascurve_f3_n3,cascurve_a3_n3, color = 'black', linestyle = 'solid', label = 'n = 3')
#ax2.plot(cascurve_f3_n4,cascurve_a3_n4, color = 'black', linestyle = 'solid', label = 'n = 4')


ax2.set_ylim(-84,-20)
ax2.set_xlim(-6000,3000)
#ax2.set_yticks([-80,-75,-70,-65,-60,-55,-50,-45,-40])
#ax2.set_yticklabels(['','','','','','','','',''])
ax2.set_xticks([-4000,-2000,0,2000])
ax2.set_xticklabels(['-4000','-2000','$+\omega_m$','+2000'])
ax2.legend()

fig2.savefig('Plot_of_Casimir_fits_n_4_d_12_4.png', dpi = 300, bbox_inches = 'tight')
fig2.savefig('Plot_of_Casimir_fits_n_4_d_12_4.svg', format = 'svg', bbox_inches = 'tight')