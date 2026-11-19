# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:19:31 2024

@author: dejonm1
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from matplotlib import gridspec
from scipy.signal import savgol_filter

#%%
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
sbnum = 1

fig1 = plt.figure(figsize = (6.5,2))

gs = gridspec.GridSpec(1, 7, width_ratios=[1,1,1,0.5,1,1,1]) 

ax1 = fig1.add_subplot(gs[0])
ax2 = fig1.add_subplot(gs[1])
ax3 = fig1.add_subplot(gs[2])
ax4 = fig1.add_subplot(gs[3])
ax5 = fig1.add_subplot(gs[4])
ax6 = fig1.add_subplot(gs[5])
ax7 = fig1.add_subplot(gs[6])

ax4.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Signal (dBm)')

mwpower = 2
sbnum = 4
ax1.plot(filtfreqs3[:,sbnum,6,mwpower],filtdata_up3[:,sbnum,6,mwpower], color = 'C0', marker = '.', alpha = 0.6, label = '-25 dBm+10', rasterized = True)
ax1.plot(filtfreqs3[:,sbnum,6,mwpower],filtdata_down3[:,sbnum,6,mwpower], color = 'C0', marker = '.', alpha = 0.6, rasterized = True)

sbnum = 2
ax2.plot(filtfreqs3[:,sbnum,6,mwpower],filtdata_up3[:,sbnum,6,mwpower], color = 'C0', marker = '.', alpha = 0.6, label = '-25 dBm+10', rasterized = True)
ax2.plot(filtfreqs3[:,sbnum,6,mwpower],filtdata_down3[:,sbnum,6,mwpower], color = 'C0', marker = '.', alpha = 0.6, rasterized = True)

sbnum = 0
ax3.plot(filtfreqs3[:,sbnum,6,mwpower][1:-1],filtdata_up3[:,sbnum,6,mwpower][1:-1], color = 'C0', marker = '.', alpha = 0.6, label = '-25 dBm+10', rasterized = True)
ax3.plot(filtfreqs3[:,sbnum,6,mwpower][1:-1],filtdata_down3[:,sbnum,6,mwpower][1:-1], color = 'C0', marker = '.', alpha = 0.6, rasterized = True)

ax4.axvspan((fcenters3[0,0]+fcenters3[0,1])/2-10000,(fcenters3[0,0]+fcenters3[0,1])/2+10000,color  = 'black', alpha = 0.2)

sbnum = 1
ax5.plot(filtfreqs3[:,sbnum,6,mwpower][1:-1],filtdata_up3[:,sbnum,6,mwpower][1:-1], color = 'C0', marker = '.', alpha = 0.6, label = '-25 dBm+10', rasterized = True)
ax5.plot(filtfreqs3[:,sbnum,6,mwpower][1:-1],filtdata_down3[:,sbnum,6,mwpower][1:-1], color = 'C0', marker = '.', alpha = 0.6, rasterized = True)

sbnum = 3
ax6.plot(filtfreqs3[:,sbnum,6,mwpower],filtdata_up3[:,sbnum,6,mwpower], color = 'C0', marker = '.', alpha = 0.6, label = '-25 dBm+10', rasterized = True)
ax6.plot(filtfreqs3[:,sbnum,6,mwpower],filtdata_down3[:,sbnum,6,mwpower], color = 'C0', marker = '.', alpha = 0.6, rasterized = True)

sbnum = 5
ax7.plot(filtfreqs3[:,sbnum,6,mwpower],filtdata_up3[:,sbnum,6,mwpower], color = 'C0', marker = '.', alpha = 0.6, label = '-25 dBm+10', rasterized = True)
ax7.plot(filtfreqs3[:,sbnum,6,mwpower],filtdata_down3[:,sbnum,6,mwpower], color = 'C0', marker = '.', alpha = 0.6, rasterized = True)

ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax5.spines['left'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax6.spines['left'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax7.spines['left'].set_visible(False)

ax1.yaxis.tick_left()
ax2.set_yticks([])
ax3.set_yticks([])
ax4.set_yticks([])
ax5.set_yticks([])
ax6.set_yticks([])
ax7.set_yticks([])
#ax6.yaxis.tick_right()
#ax1.tick_params(labelright='off')
#ax2.tick_params(labelright='off')

ax1.set_xticks([fcenters3[0,4]])
ax1.set_xticklabels(['$-3\omega_m$'])

ax2.set_xticks([fcenters3[0,2]])
ax2.set_xticklabels(['$-2\omega_m$'])

ax3.set_xticks([fcenters3[0,0]])
ax3.set_xticklabels(['$-\omega_m$'])

ax4.set_xticks([(fcenters3[0,0]+fcenters3[0,1])/2])
ax4.set_xticklabels(['$\omega_c$'])

ax5.set_xticks([fcenters3[0,1]])
ax5.set_xticklabels(['$+\omega_m$'])

ax6.set_xticks([fcenters3[0,3]])
ax6.set_xticklabels(['$+2\omega_m$'])

ax7.set_xticks([fcenters3[0,5]])
ax7.set_xticklabels(['$+3\omega_m$'])

fig1.subplots_adjust(wspace=0.05, hspace=0)

ax1.set_ylim((-84,-10))
ax2.set_ylim((-84,-10))
ax3.set_ylim((-84,-10))
ax4.set_ylim((-84,-10))
ax5.set_ylim((-84,-10))
ax6.set_ylim((-84,-10))
ax7.set_ylim((-84,-10))

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
ax2.plot((-d, +d), (-d, +d), **kwargs)

kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
ax2.plot((1-d, 1+d), (-d, +d), **kwargs)
ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
ax3.plot((-d, +d), (1-d, 1+d), **kwargs)
ax3.plot((-d, +d), (-d, +d), **kwargs)

kwargs = dict(transform=ax3.transAxes, color='k', clip_on=False)
ax3.plot((1-d, 1+d), (-d, +d), **kwargs)
ax3.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

kwargs.update(transform=ax4.transAxes)  # switch to the bottom axes
ax4.plot((-d, +d), (1-d, 1+d), **kwargs)
ax4.plot((-d, +d), (-d, +d), **kwargs)

kwargs = dict(transform=ax4.transAxes, color='k', clip_on=False)
ax4.plot((1-d, 1+d), (-d, +d), **kwargs)
ax4.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

kwargs.update(transform=ax5.transAxes)  # switch to the bottom axes
ax5.plot((-d, +d), (1-d, 1+d), **kwargs)
ax5.plot((-d, +d), (-d, +d), **kwargs)

kwargs = dict(transform=ax5.transAxes, color='k', clip_on=False)
ax5.plot((1-d, 1+d), (-d, +d), **kwargs)
ax5.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

kwargs.update(transform=ax6.transAxes)  # switch to the bottom axes
ax6.plot((-d, +d), (1-d, 1+d), **kwargs)
ax6.plot((-d, +d), (-d, +d), **kwargs)

kwargs = dict(transform=ax6.transAxes, color='k', clip_on=False)
ax6.plot((1-d, 1+d), (-d, +d), **kwargs)
ax6.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

kwargs.update(transform=ax7.transAxes)  # switch to the bottom axes
ax7.plot((-d, +d), (1-d, 1+d), **kwargs)
ax7.plot((-d, +d), (-d, +d), **kwargs)

#fig1.savefig("All_sidebands.png", dpi = 300, bbox_inches = 'tight')
#fig1.savefig("All_sidebands.svg", format = 'svg', bbox_inches = 'tight', dpi = 300)

#%%
#And now to make them overlap.

fig2 = plt.figure(figsize = (6.4,3))
ax2 = fig2.add_subplot(111)

ax2.set_xlabel("Scaled frequency (Hz)")
ax2.set_ylabel("Scaled signal (dB)")

cavfreq = (fcenters3[0,0]+fcenters3[0,1])/2

ROI1 = 200
red1_f = savgol_filter(-filtfreqs3[1:-1,0,6,mwpower]+fcenters3[0,0],5,1)
red1_up_d = savgol_filter(filtdata_up3[1:-1,0,6,mwpower],5,1)
red1_do_d = savgol_filter(filtdata_down3[1:-1,0,6,mwpower],5,1)

ax2.plot(red1_f[::3],red1_up_d[::3]-2.5, color = 'red', marker = '.', alpha = 0.8, label = '$-\omega_\mathrm{m}$')
ax2.plot(red1_f[ROI1::3],red1_do_d[ROI1::3]-2.5, color = 'red', marker = '.', alpha = 0.8)
ax2.plot(red1_f[0:ROI1:3],red1_do_d[0:ROI1:3]-2.5, color = 'red', marker = '.', alpha = 0.05)

#ax2.plot(-filtfreqs3[1:-1,0,6,mwpower]+fcenters3[0,0],filtdata_up3[1:-1,0,6,mwpower]-2.5, color = 'red', marker = '.', alpha = 0.4, rasterized = True)
#ax2.plot(-filtfreqs3[ROI1:-1,0,6,mwpower]+fcenters3[0,0],filtdata_down3[ROI1:-1,0,6,mwpower]-2.5, color = 'red', marker = '.', alpha =0.4, rasterized = True)
#ax2.plot(-filtfreqs3[0:ROI1,0,6,mwpower]+fcenters3[0,0],filtdata_down3[0:ROI1,0,6,mwpower]-2.5, color = 'red', marker = '.', alpha = 0.01, rasterized = True)

blue1_f = savgol_filter(filtfreqs3[1:-1,1,6,mwpower]-fcenters3[0,1],5,1)
blue1_up_d = savgol_filter(filtdata_up3[1:-1,1,6,mwpower],5,1)
blue1_do_d = savgol_filter(filtdata_down3[1:-1,1,6,mwpower],5,1)

ax2.plot(blue1_f[::3],blue1_up_d[::3], color = 'blue', marker = '.', alpha = 0.8, label = '$+\omega_\mathrm{m}$')
ax2.plot(blue1_f[::3],blue1_do_d[::3], color = 'blue', marker = '.', alpha = 0.8)

#ax2.plot(filtfreqs3[1:-1,1,6,mwpower]-fcenters3[0,1],filtdata_up3[1:-1,1,6,mwpower], color = 'blue', marker = '.', alpha = 0.4, rasterized = True)
#ax2.plot(filtfreqs3[1:-1,1,6,mwpower]-fcenters3[0,1],filtdata_down3[1:-1,1,6,mwpower], color = 'blue', marker = '.', alpha = 0.4, rasterized = True)

ROI2 = 160
red2_f = savgol_filter((-filtfreqs3[1:-1,2,6,mwpower]+fcenters3[0,2])*0.5,5,1)
red2_up_d = savgol_filter((filtdata_up3[1:-1,2,6,mwpower])*0.5,5,1)
red2_do_d = savgol_filter((filtdata_down3[1:-1,2,6,mwpower])*0.5,5,1)

ax2.plot(red2_f[::5],red2_up_d[::5]-1, color = 'orange', marker = 'x', alpha = 0.6, label = '$-2\omega_\mathrm{m}$')
ax2.plot(red2_f[ROI2::5],red2_do_d[ROI2::5]-1, color = 'orange', marker = 'x', alpha = 0.6)
ax2.plot(red2_f[0:ROI2:5],red2_do_d[0:ROI2:5]-1, color = 'orange', marker = 'x', alpha = 0.05)

#ax2.plot((-filtfreqs3[1:-1,2,6,mwpower]+fcenters3[0,2])*0.5,(filtdata_up3[1:-1,2,6,mwpower])*0.5-1, color = 'red', marker = 'x', alpha = 0.2, rasterized = True)
#ax2.plot((-filtfreqs3[ROI2:-1,2,6,mwpower]+fcenters3[0,2])*0.5,(filtdata_down3[ROI2:-1,2,6,mwpower])*0.5-1, color = 'red', marker = 'x', alpha = 0.2, rasterized = True)
#ax2.plot((-filtfreqs3[0:ROI2,2,6,mwpower]+fcenters3[0,2])*0.5,(filtdata_down3[0:ROI2,2,6,mwpower])*0.5-1, color = 'red', marker = 'x', alpha = 0.01, rasterized = True)

ROI3 = 320
blue2_f = savgol_filter((filtfreqs3[1:-1,3,6,mwpower]-fcenters3[0,3])*0.5,5,1)
blue2_up_d = savgol_filter((filtdata_up3[1:-1,3,6,mwpower])*0.5,5,1)
blue2_do_d = savgol_filter((filtdata_down3[1:-1,3,6,mwpower])*0.5,5,1)

ax2.plot(blue2_f[::5],blue2_up_d[::5]-0.5, color = 'teal', marker = 'x', alpha = 0.6, label = '$+2\omega_\mathrm{m}$')
ax2.plot(blue2_f[1:ROI3:5],blue2_do_d[1:ROI3:5]-0.5, color = 'teal', marker = 'x', alpha = 0.6)
ax2.plot(blue2_f[ROI3::5],blue2_do_d[ROI3::5]-0.5, color = 'teal', marker = 'x', alpha = 0.05)

#ax2.plot((filtfreqs3[1:-1,3,6,mwpower]-fcenters3[0,3])*0.5,(filtdata_up3[1:-1,3,6,mwpower])*0.5-0.5, color = 'blue', marker = 'x', alpha = 0.2, rasterized = True)
#ax2.plot((filtfreqs3[1:ROI3,3,6,mwpower]-fcenters3[0,3])*0.5,(filtdata_down3[1:ROI3,3,6,mwpower])*0.5-0.5, color = 'blue', marker = 'x', alpha = 0.2, rasterized = True)
#ax2.plot((filtfreqs3[ROI3:,3,6,mwpower]-fcenters3[0,3])*0.5,(filtdata_down3[ROI3:,3,6,mwpower])*0.5-0.5, color = 'blue', marker = 'x', alpha = 0.01, rasterized = True)

ROI4 = 200
red3_f = savgol_filter((-filtfreqs3[1:-1,4,6,mwpower]+fcenters3[0,4])*0.33,5,1)
red3_up_d = savgol_filter((filtdata_up3[1:-1,4,6,mwpower])*0.33,5,1)
red3_do_d = savgol_filter((filtdata_down3[1:-1,4,6,mwpower])*0.33,5,1)

ax2.plot(red3_f[::5],red3_up_d[::5]-6, color = 'gold', marker = 'v', alpha = 0.4, label = '$-3\omega_\mathrm{m}$')
ax2.plot(red3_f[ROI4::5],red3_do_d[ROI4::5]-6, color = 'gold', marker = 'v', alpha = 0.4)
ax2.plot(red3_f[0:ROI4:5],red3_do_d[0:ROI4:5]-6, color = 'gold', marker = 'v', alpha = 0.03)

#ax2.plot((-filtfreqs3[1:-1,4,6,mwpower]+fcenters3[0,4])*0.33,(filtdata_up3[1:-1,4,6,mwpower])*0.33-6, color = 'red', marker = 'v', alpha = 0.1, rasterized = True)
#ax2.plot((-filtfreqs3[ROI4:-1,4,6,mwpower]+fcenters3[0,4])*0.33,(filtdata_down3[ROI4:-1,4,6,mwpower])*0.33-6, color = 'red', marker = 'v', alpha = 0.1, rasterized = True)
#ax2.plot((-filtfreqs3[0:ROI4,4,6,mwpower]+fcenters3[0,4])*0.33,(filtdata_down3[0:ROI4,4,6,mwpower])*0.33-6, color = 'red', marker = 'v', alpha = 0.01, rasterized = True)

ROI5 = 300
blue3_f = savgol_filter((filtfreqs3[1:-1,5,6,mwpower]-fcenters3[0,5])*0.33,5,1)
blue3_up_d = savgol_filter((filtdata_up3[1:-1,5,6,mwpower])*0.33,5,1)
blue3_do_d = savgol_filter((filtdata_down3[1:-1,5,6,mwpower])*0.33,5,1)

ax2.plot(blue3_f[::5],blue3_up_d[::5]-5.5, color = 'cyan', marker = 'v', alpha = 0.3, label = '$+3\omega_\mathrm{m}$')
ax2.plot(blue3_f[1:ROI5:5],blue3_do_d[1:ROI5:5]-5.5, color = 'cyan', marker = 'v', alpha = 0.3)
ax2.plot(blue3_f[ROI5::5],blue3_do_d[ROI5::5]-5.5, color = 'cyan', marker = 'v', alpha = 0.03)

#ax2.plot((filtfreqs3[1:-1,5,6,mwpower]-fcenters3[0,5])*0.33,(filtdata_up3[1:-1,5,6,mwpower])*0.33-5.5, color = 'blue', marker = 'v', alpha = 0.1, rasterized = True)
#ax2.plot((filtfreqs3[1:ROI5,5,6,mwpower]-fcenters3[0,5])*0.33,(filtdata_down3[1:ROI5,5,6,mwpower])*0.33-5.5, color = 'blue', marker = 'v', alpha = 0.1, rasterized = True)
#ax2.plot((filtfreqs3[ROI5:,5,6,mwpower]-fcenters3[0,5])*0.33,(filtdata_down3[ROI5:,5,6,mwpower])*0.33-5.5, color = 'blue', marker = 'v', alpha = 0.01, rasterized = True)

ax2.set_xlim(-10000,3000)
ax2.set_ylim(-53,-20)

ax2.legend(ncol = 3, fontsize = 'small', loc = 'lower right')

ax2.set_xticks([-10000,-7500,-5000,-2500,0,2500])
ax2.set_xticklabels(['-10000','-7500','-5000','-2500','$\omega_\mathrm{m}$','+2500'])

fig2.savefig("Overlapping_sidebands_v2.png", dpi = 300, bbox_inches = 'tight')
fig2.savefig("Overlapping sidebands_v2.svg", format = 'svg', bbox_inches = 'tight')
