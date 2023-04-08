#!/usr/bin/env python3
# coding=UTF-8

"""
Created on Fri May 10 09:26:06 2019

@author: 20180239
"""
import sys
import numpy as np
import pandas as pd
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
#matplotlib.use('Agg')


#Initialization
filename='COHPCAR.lobster'
imgname='cohp.png'
erange = [-15.00, 2.00]
spin = 'up'
cc1 = 'red' #'red' or 'cyan'
cc2 = 'pink' # 'pink' or 'lightblue'

#Set Data 
argulen=len(sys.argv)
if argulen > 1:
    for arg in sys.argv[1:]:
        if arg.startswith(('fname=', 'filename=')):
            arg = arg.split('=')
            filename = arg[-1]
        elif arg.startswith(('imgname=', 'imagename=')):
            arg = arg.split('=')
            imgname = arg[-1]
        elif arg.startswith('spin='):
            arg = arg.split('=')
            spin = arg[-1]     
        elif arg.startswith('erange='):
            arg = arg.split('=')[-1].strip('(').strip(')').split(',')
            erange = float(arg[0]) , float(arg[-1]) 
else:
    print ("An exception keyword occurred") 
    #print ('Try syntax: pypltpcohp.py [options]')
    print ("Try syntax: pypltpcohp.py fname=COHPCAR.lobster imgname=COHP.png erange=-8,2")
    raise SyntaxError('invalid syntax')
            
#skip lines
txt = open(filename, 'r')
skipwords= ('ESCALE', '    ', 'Average', 'No.')
skiprows = 0
for line in txt:
    if line.startswith(skipwords): 
        skiprows += 1
    else:  
        break
txt.close()

#Read data
usecols = [0,1,3]
df = pd.read_csv(filename, skiprows=skiprows, header=None, usecols=usecols, delim_whitespace=True)
df.columns = ['energy', 'up', 'down' ]

#Set Figure size
fig, ax = plt.subplots(nrows=1, ncols=1,  sharey=True, figsize=(3,6), dpi=200)
fig.subplots_adjust(wspace=0.05)

#Set plot data
if spin == 'up':
    x = -1*df.up 
else:
    x = -1*df.down
y = df.energy


#Plot
#emask = (y1 >= erange[0]) & (y1 <= erange[1])
xmax = np.max([abs(np.min(x)), abs(np.max(x))])
xrange = [-xmax*1.10, xmax*1.10]
yrange = erange
ax.plot(x, y, alpha=1.00, label=None, color='gray', linewidth=0.50)
ax.fill_betweenx(y, x, -0.000,  where=y <= 0, interpolate=True, color=cc1, alpha=0.80) 
ax.fill_betweenx(y, x, +0.000,  where=y >= 0, interpolate=True, color=cc2, alpha=0.80) 

#set range 
ax.set_xlim(xrange[0], xrange[1])
ax.set_ylim(yrange[0], yrange[1])

#Set annotate
text_Ef = r'$\epsilon_F$' 
xycoords = 'axes fraction'
bbox=dict(boxstyle='square', pad=0.5, fc='w', ec='none', alpha=1.00)    
#ax.annotate(text_Ef, xy=(0.750, 0.785),  xycoords=xycoords, color= 'k',bbox=bbox, fontsize=20)
ax.annotate('anti-bonding', xy=(0.02, 0.010), xycoords='axes fraction', color= 'k',fontsize=10)
ax.annotate('bonding', xy=(0.70, 0.010), xycoords='axes fraction', color= 'k',fontsize=10)
ax.set_xlabel('-pCOHP', fontsize=15)
ax.set_ylabel('Energy (eV)', fontsize=15)
#ax.set_xticks(ticks=[])

#set vertical  and horizontal lines
ax.axvline(x=0.00, ymin=0.00, ymax=1.00, linewidth=1.20, linestyle='-', color='black')
ax.axhline(y=0.00, xmin=0.00, xmax=1.00, linewidth=1.20, linestyle='-', color='black') 

#Save image
fig.savefig(imgname)
