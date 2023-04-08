#!/usr/bin/env python3 
# -*- coding: utf-8 -*-

"""
Created on Tue Dec 11 09:56:16 2018
@author: mwchang
syntax 1: setidm.py [POSCAR/CONTCAR]
"""

import sys
import subprocess
import numpy as np
import vasp_modules.vasp_io2 as vio

if len(sys.argv) == 1:
    carfile = 'POSCAR'
else:
    carfile = sys.argv[1] 
    
vibinfo = vio.extra_vibinfo('OUTCAR') 
posmodes = []; negmodes = []
for line in vibinfo:
    if 'f  =' in line:
        line = line.replace(' f  =', 'f =')
        mode = line.split('=')[0].strip()
        posmodes.append(mode)
        print(line)
    elif ' f/i=' in line:
        line = line.replace(' f/i=','i =')
        mode = line.split('=')[0].strip()
        negmodes.append(mode)
        print(line)
        
while True:
    defmode = negmodes[-1]
    mode = input('Please select a vibrational mode (default: %s): ' %(defmode)).strip()
    if mode == '':
        mode = defmode
        
    if mode not in posmodes and mode not in negmodes:
        print ('The %s mode is not a reasonable mode.' %(mode))  
        continue
    else:
        print ('The %s mode has been selected.' %(mode))  
        break
    
molecule = vio.read_poscar(carfile)
molecule.write(filename='POSCAR.vasp', format='vasp')
dynmat = vio.get_dymatrix(vibinfo)   
eigvects = np.array(dynmat[mode])
f = open('POSCAR.vasp', 'a')
f.write(" ! here we define trial unstable direction:\n")
[f.write(" %18.15f %18.15f %18.15f\n" %(eigvect[0], eigvect[1], eigvect[2])) for eigvect in eigvects]
f.close()
print('POSCAR.vasp for improved dimer method has been generated')      
subprocess.call("sed -i 's/IBRION.*/IBRION = 44/' INCAR", shell=True) 
subprocess.call("sed -i 's/POTIM.*/POTIM = 0.50/' INCAR", shell=True) 
subprocess.call("sed -i 's/NSW.*/NSW = 100/' INCAR", shell=True) 
