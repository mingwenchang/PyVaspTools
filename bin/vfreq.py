#!/usr/bin/env python3
# coding=UTF-8



"""
NAME
       vfreq.py - extra vibrational frequency information from OUTCAR and 
                  calculate zero-point energy  
        
SYNTAX
        vfreq.py
 
DESCRIPTION
            extra vibrational frequency information from OUTCAR and 
            calculate zero-point energy  
        
        output :
            zpe.dat
 
            
DEVELOPER: 
    
    Dr. Ming-Wen Chang
    E-mail: ming.wen.c@gmail.com

"""

import sys
import vasp_modules.assister as ast
from vasp_modules.vasp_io2 import extra_vibinfo, get_freqs 
 

vibinfo = extra_vibinfo()
freqs = get_freqs(vibinfo, 'meV')
argulen=len(sys.argv)
  
if argulen == 1:
    zpe = 0.001 *0.5* sum(i for i in freqs if i > 0) #unit in eV
    freqinfo = [line for line in vibinfo if 'f' in line]
    with open("zpe.dat", "w") as txt:
        for line in freqinfo:
            print (line)
            txt.write(line + '\n')
        print ('Zero Point Energy: %s eV' %(zpe))  
        txt.write('Zero Point Energy: %s eV' %(zpe))  

else:
    raise SyntaxError('Syntax: vfreq.py')
    

