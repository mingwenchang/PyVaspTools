#!/usr/bin/env python3
# coding=UTF-8

"""
NAME
        kicknock.py - eliminate imaginary frequencies
        
SYNTAX
        kicknock.py OUTCAR CONTCAR

DESCRIPTION
        This script is used to eliminate imaginary frequencies for vasp calculations.

        It reads the dynamical matrix from OUTCAR (generted by a freq calcualtion) and 
        then adds eigenvalues of the dynamical matrix to coordinations of CONTCAR
 
        output :
            new.vasp
 
            
DEVELOPER: 
    
    Dr. Ming-Wen Chang
    E-mail: ming.wen.c@gmail.com

"""

import sys
import numpy as np
import vasp_modules.vasp_io2 as vio


if __name__ == "__main__":
    if len(sys.argv) == 3:
        outcar = sys.argv[1]
        contcar = sys.argv[2]
    else:
        raise SyntaxError('Try Syntax: kicknock.py OUTCAR CONTCAR')

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
        factor = input('Please assign a factor (default: 0.30): ').strip()
        if factor == '':
            factor = 0.30  
        try:
            factor = float(factor)
            print ('%s has been used.' %(factor))
            break
        except ValueError:
            print ('The factor should be a int/float')
            continue

    molecule = vio.read_poscar(contcar)
    dynmat = vio.get_dymatrix(vibinfo)
    shfmat = factor * np.array(dynmat[mode])
    newpos = molecule.positions + shfmat
    molecule.set_atomic_positions(newpos)
    molecule.write(filename='new', format='vasp')
    print('New coordinates have been successfully generated and saved as new.vasp')
    print('Please re-perform an optimization calculation using the new coordinates')

   
