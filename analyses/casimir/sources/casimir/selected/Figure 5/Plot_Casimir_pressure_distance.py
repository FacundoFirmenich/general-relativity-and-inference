# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:16:26 2024

@author: dejonm1
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Data = pd.read_hdf('Casimir_pressure_dist_v3.h5')
#%%
Casmodel = -155.22*147000/((Data.iloc[:,0])**3.2)

fig1 = plt.figure(figsize = (6,4.8))
ax1 = fig1.add_subplot(111)

ax1.set_xlabel('Distance (nm)')
ax1.set_ylabel('Pressure (kPa)')

ax1.plot(Data.iloc[:,0],np.pi*Data.iloc[:,1]/1e3, marker = 'o', alpha = 0.6, label = 'Lifshitz formula')
ax1.plot(Data.iloc[:,0],np.pi*Casmodel/1e3, color = 'black', label = 'Fit $\propto d^{-3.2}$')

ax1.legend(fontsize = 'small',frameon = False)
ax1.set_xlim((9,30))

fig1.savefig('Plot_casimir_pressure_full_v2.png', dpi = 300, bbox_inches = 'tight')
fig1.savefig('Plot_casimir_pressure_full_v2.svg', format = 'svg', bbox_inches = 'tight')