#!/usr/bin/env python3
# coding=UTF-8
"""
Created on Thu Sep 12 09:44:40 2019
Name:chgdiff.py
Developer: Ming-Wen Chang
E-mail: ming.wen.c@gmail.com

Reading CHGCAR files and calculating their difference
syntax: chgdiff.py CHGCAR1 CHGCAR2 [CHGCAR3 ...]
output: CHGCAR_DIFF ( = CHGCAR1 - CHGCAR2 - [CHGCAR3 - ....] )
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
    if nargs < 3:
        print ('You have to specify two names at least on command line')
        print ('Syntax: chgdiff.py CHGCAR1 CHGCAR2 [CHGCAR3 CHGCAR4 etc.] ')
        print ('Output: CHGCAR_DIFF ( = CHGCAR1 - CHGCAR2 - [CHGCAR3 - ....])')
        exit()
    
    names =sys.argv[1:]
    for idx, name in  enumerate(names):
        print ("Reading charge density data from %s file...." %(name))
        chgobj = VaspChargeDensity(filename=name)
        if idx == 0:
            print ('The reference charge density is taken from %s file' %(name))
            chgdiff = chgobj.chg[-1] 
            atomsobj = chgobj.atoms[-1]
            
        else:
            print ("Subtracting data from the reference charge density data .....")
            try:
                chgdiff -= chgobj.chg[-1]
            except (IOError, ValueError, IndexError):
                print("IOError ---> Probably %s file is empty" % (name))
                print ("ValueError or IndexError ---> Two sets of data are not with the same grid.")
   
    print ("Subtracting is done") 
    
    chgobj = VaspChargeDensity(filename=None)
    chgobj.atoms.append(atomsobj)
    chgobj.chg.append(chgdiff)
    write_chgcar(chgobj, fname='CHGCAR_DIFF', vasp5=True)
    
 
 
