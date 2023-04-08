# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:03:27 2020

@author: 20180239
"""

import numpy as np
import matplotlib.pyplot as plt

 
 
rad = np.pi/180
 
c3h6 = np.array([1.00, 1.00, 1.00, 1.00])
theata = np.linspace(15 * rad, 175 * rad, 4)
 
       
#Figure Setting 
fig, ax = plt.subplots(figsize=(16,9), dpi=72, nrows=1, ncols=1, subplot_kw={'projection':'polar'})      

 
#Ploting
ax.scatter(theata, c3h6, 
           s=400*np.abs(c3h6),
           label = 'C$_3$H$_6$',
           c = 'blue',
           alpha = 0.50,
           edgecolor = 'k')
 



#Set range
ax.set_thetamin(-10)
ax.set_thetamax(190)
#ax.set_rmin(-1)
#ax.set_rmax(5.5)
ax.tick_params(labelsize = 20)
ax.set_xticks([0*rad, 22.5*rad, 67.5*rad, 112.5*rad,  157.5*rad, 180*rad])
ax.set_xticklabels(['', 'Pt$_8$', 'Pt$_7$Sn', 'Pt$_6$Sn$_2$', 'Pt$_4$Sn$_4$', ''])
ax.axvline(x=40*rad, linewidth=1.50, linestyle='-', color='#393E46', alpha=0.50)
ax.axvline(x=90*rad, linewidth=1.50, linestyle='-', color='#393E46', alpha=0.50)           
ax.axvline(x=140*rad, linewidth=1.50, linestyle='-', color='#393E46', alpha=0.50) 
ax.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.) 
#ax.grid(linestyle='--', color='#393E46', alpha=1.00)


 