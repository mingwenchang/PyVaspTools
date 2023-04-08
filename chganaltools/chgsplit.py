#!/usr/bin/env python3
# coding=UTF-8
"""
Created on Thu Sep 12 09:44:40 2019
Name:chgsplit.py
Developer: Ming-Wen Chang
E-mail: ming.wen.c@gmail.com

Read CHGCAR files and split it into CHGCAR and CHGCAR_MAG
syntax: chgdiff.py CHGCAR 
output: CHGCAR_TOT and CHGCAR_MAG
"""

import sys
from ase.calculators.vasp import VaspChargeDensity 

def write_chgcar(chgobj, fname='CHGCAR_NEW', vasp5=True):
    print ("Writing charge density data to %s file..." %(fname))
    if vasp5:
        import ase.io.vasp as aiv
        chg = chgobj.chg[-1]
        atomsobj = chgobj.atoms[-1]
        vol = atomsobj.get_volume()
        f = open(fname, 'w')
        aiv.write_vasp(f, atomsobj, direct=True, long_format=False,vasp5=True)
        f.write('\n')
        #Write ngrids 
        for ngrid in chg.shape:
            f.write(' %4i' %ngrid)
        #Write charge density data 
        f.write('\n')
        chgobj._write_chg(f, chg, vol, format='chgcar')
        f.close()
    else: #for vasp 4
        chgobj.write(filename=fname, format='chgcar')
        
    print ('Charge density data have been saved to %s file.' %(fname))
    

if __name__ == "__main__":
    nargs  = len(sys.argv)
    if nargs < 2:
        print ('You have to specify a name at least on command line')
        print ('Syntax: chgsplit.py CHGCAR ')
        print ('Output: CHGCAR_TOT and CHGCAR_MAG')
        exit()
    
    name = sys.argv[1]
    print ("Reading charge density data from %s file...." %(name))
    chgobj = VaspChargeDensity(filename=name)
    ispin = chgobj.is_spin_polarized()
    if ispin:
        write_chgcar(chgobj, fname='CHGCAR_TOT', vasp5=True)
        chgobj.chg[-1]  = chgobj.chgdiff[-1]   
        write_chgcar(chgobj, fname='CHGCAR_MAG', vasp5=True)
    else:
        raise ValueError ('No spin-polarized calculation was performed in %s file' %(name))
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    