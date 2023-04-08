#!/usr/bin/env python3
# coding=UTF-8

"""
NAME
        setfreq.py - set vibrational frequency calculation for vasp
        
SYNTAX
        setfreq.py -n=int CONTCAR
           => set Selective dynamic T, T, T for LAST n atoms,
              and F, F, F for others. 

        setfreq.py -h=int CONTCAR 
            => set Selective dynamic T, T, T for FIRST n atoms,
               and F, F, F for others. 
                
        setfreq.py -z=float CONTCAR 
            => set set Selective dynamic T, T, T for atoms with 
               z-coordinate > z, and F, F, F for others.
               
        setfreq.py -s=int1, int2, int3,... CONTCAR
            => set set Selective dynamic T, T, T for atomi, atomj, atomk, ...
               and F, F, F for others.
 
DESCRIPTION
        Re-set selective dynamic on atoms in POSCAR or CONTCAR and 
        keywords in INCAR to IBRION = 5, POTIM = 0.01 and NSW = 1
        
        output :
            INCAR and POSCAR
 
            
DEVELOPER: 
    
    Dr. Ming-Wen Chang
    E-mail: ming.wen.c@gmail.com

"""

import sys
import subprocess
import numpy as np
import vasp_modules.vasp_io2 as vio


if __name__ == "__main__":
    
    z = None #relaxing atoms with z > zmin
    n = None #relaxing last n atoms 
    h = None #relaxing start n atoms
    s = None   #relaxing specific atoms 
    contcar = 'CONTCAR'
    argulen=len(sys.argv)
    
    for arg in sys.argv[1:]:
        path = arg
        arg = arg.lower()
        if 'contcar' in arg:
            contcar = path
        elif 'poscar' in arg:
            contcar = path
        elif arg.startswith('-z='): 
            arg = arg.strip().split('=')
            z = float(arg[-1]) 
        elif arg.startswith('-n='):
            arg = arg.strip().split('=')
            n = int(arg[-1])
        elif arg.startswith('-h='):
            arg = arg.strip().split('=')
            h = int(arg[-1])
        elif arg.startswith('-s='):
            arg = arg.strip().split('=')[-1].split(',')
            s = [int(i) for i in arg]
        else:
            print ("An exception keyword occurred")
            print ('Try syntax: setfreq.py -n=int CONTCAR')
            print ('Try syntax: setfreq.py -h=int CONTCAR')
            print ('Try syntax: setfreq.py -z=float CONTCAR')
            print ('Try syntax: setfreq.py -s=int1, int2, int3,... CONTCAR')
            raise SyntaxError('invalid syntax') 
            
    molobj = vio.read_poscar(contcar)
    natoms = molobj.get_total_atoms()
    clength = molobj.get_cell_lengths()[-1]
    molobj.constraints = np.array([['F', 'F', 'F']]*natoms)
    
    if z is not None:
        z *= clength
        for i, pi in enumerate(molobj.positions):
            if pi[-1] >= z:
                molobj.constraints[i] = np.array(['T', 'T', 'T']) 
    elif s is not None:
        for i in s:
            molobj.constraints[i-1] = np.array(['T', 'T', 'T']) 
    elif n is not None:
        molobj.constraints[-n:] = np.array([['T', 'T', 'T']]*n)
    elif h is not None:
        molobj.constraints[0:h] = np.array([['T', 'T', 'T']]*h)
    else:
        raise SyntaxError('invalid syntax') 
    
    molobj.write(filename='POSCAR', format='vasp')
    
    #Modify INCAR
    subprocess.call("sed -i 's/.*NSW.*/NSW = 1/' INCAR", shell=True) 
    subprocess.call("sed -i 's/.*NWRITE.*/NWRITE = 3/' INCAR", shell=True) 
    subprocess.call("sed -i 's/.*IBRION.*/IBRION = 5/' INCAR", shell=True) 
    subprocess.call("sed -i 's/.*POTIM.*/POTIM = 0.01/' INCAR", shell=True) 
    print('POSCAR and INCAR are ready for the freq calculation.')
