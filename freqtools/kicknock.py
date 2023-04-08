#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:56:16 2018
@author: mwchang
syntax 1: kicknock.py 
syntax 2: kicknock.py OUTCAR CONTCAR
"""

import sys
import numpy as np
import vasp_modules.vasp_io2 as vio

if len(sys.argv) == 1:
    outcar = 'OUTCAR'
    contcar = 'CONTCAR'
elif len(sys.argv) == 3:
    outcar = sys.argv[1]
    contcar = sys.argv[2]
else:
    raise SyntaxError('please use command: kicknock.py or kicknock.py "OUTCAR" "CONTCAR"')

vibinfo = vio.extra_vibinfo(outcar) 
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
    
while True:        
    factor = input('Please assign a factor (default: 0.15): ').strip()
    if factor == '':
        factor = 0.15  
    try:
        factor = float(factor)
        print ('%s has been used.' %(factor))
        break
    except ValueError:
        print ('The factor should be a number')
        continue

molecule = vio.read_poscar(contcar)
dynmat = vio.get_dymatrix(vibinfo)
shfmat = factor * np.array(dynmat[mode])
newpos = molecule.positions + shfmat
molecule.set_atomic_positions(newpos)
molecule.write(filename='POSCAR.vasp', format='vasp')
print('Coordinates have been saved as POSCAR.vasp')

   
