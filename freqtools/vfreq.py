#!/usr/bin/env python3
# coding=UTF-8

"""
Syntaxs
    (1) vfreq 
    (2) vfreq -a 

"""
import sys
from vasp_modules.vasp_io2 import extra_vibinfo, get_freqs 
 

vibinfo = extra_vibinfo()
freqs = get_freqs(vibinfo, 'meV')
argulen=len(sys.argv)
  
if argulen == 1:
    zpe = 0.001 *0.5* sum(i for i in freqs if i > 0) #unit in eV
    freqinfo = [line for line in vibinfo if 'f' in line]
    for line in freqinfo:
        print (line)
    print ('Zero Point Energy: %s eV' %(zpe))  
elif argulen == 2 and sys.argv[1].lower() == '-a':
    zpe = 0.001 *0.5* sum(i for i in freqs if i > 0) #unit in eV
    for line in vibinfo:
          print (line)
    print ('Zero Point Energy: %s eV' %(zpe))
else:
    raise SyntaxError('invalid syntax')
    
