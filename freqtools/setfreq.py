#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:56:16 2018
@author: mwchang
syntax 1: setfreq.py 
syntax 2: setfreq.py --lines=n  CONTCAR
"""

import sys
import subprocess
import numpy as np
import vasp_modules.vasp_io2 as vio

if len(sys.argv) == 1:
    n = None
    contcar = 'CONTCAR'
elif len(sys.argv) == 3:
    n = int(sys.argv[1].strip('-').strip('lines='))
    contcar = sys.argv[2]
else:
    raise SyntaxError('please use command: setfreq.py or setfreq.py --lines=n  "CONTCAR"')

molobj = vio.read_poscar(contcar)
ntotal = molobj.get_total_atoms()
molobj.constraints = np.array([['F', 'F', 'F']]*ntotal)
if type(n) is int:
    molobj.constraints[-n:] = np.array([['T', 'T', 'T']]*n)
molobj.write(filename='POSCAR', format='vasp')

#Modify INCAR
subprocess.call("sed -i 's/.*NSW.*/NSW = 1/' INCAR", shell=True) 
subprocess.call("sed -i 's/.*NWRITE.*/NWRITE = 3/' INCAR", shell=True) 
subprocess.call("sed -i 's/.*IBRION.*/IBRION = 5/' INCAR", shell=True) 
subprocess.call("sed -i 's/.*POTIM.*/POTIM = 0.01/' INCAR", shell=True) 
print('POSCAR and INCAR for the freq calculation have been generated.')