#!/usr/bin/env python3
# coding=UTF-8
"""
Created on Thu Sep 12 09:44:40 2019
Name:chgave.py
Developer: Ming-Wen Chang
E-mail: ming.wen.c@gmail.com

Read CHGCAR/LOCPOT files and integrating/averaging density/potential 
along the x/y/z-direction
syntax: chgave.py -[x/y/z] CHGCAR/LOCPOT 
output: CHGCAR_[X/Y/Z] or LOCPOT_[X/Y/Z] 
"""

import sys
import numpy as np
from ase.calculators.vasp import VaspChargeDensity 

if __name__ == "__main__":
    """
    Format of CHGCAR 
    
    In CHGCAR, total charge density is preserved as volumetric data with the
    fine FFT-grid of ngx by ngy by ngz. Each element in volumetric data 
    represents charge density in the voxel of dV = dx * dy * dz.  One can get 
    total charge in the simulation box by integrating all voxels.  
    (cf. https://cms.mpi.univie.ac.at/wiki/index.php/CHGCAR)
    
    Format of LOCPOT is similar to that of CHGCAR, but each element in 
    volumetric data represents total LOCAL potential (in eV) in the voxel of 
    dV = dx * dy * dz. An average potential of the simulation box can be 
    obtained by integrating all voxels then divide the value by total number of
    voxels. 
    (cf. https://cms.mpi.univie.ac.at/vasp/vasp/LOCPOT_file.html)
    
    The relations decribe below: 
    
    Lx = length of the simulation box in x-direction
    Ly = length of the simulation box in y-direction
    Lz = length of the simulation box in z-direction
    V = volume of the simulation box (Lx * Ly * Lz) 
    
    ngx = number of grids in x-direction
    ngy = number of grids in y-direction
    ngz = number of grids in z-direction 
    ngt = ngx * ngy * ngz (number of voxels in the simulation box) 
    
    dx = Lx/ngx 
    dy = Ly/ngy
    dz = Lz/ngz
    dV = dx * dy * dz #the volume of a voxel 
    
    rho_ijk = charge density in the voxel of ijk 
    rho = sum_ijk(rho_ijk) #total charge density in the simulation box
    
    rho_ijk *dV = number of electrons in the voxel of ijk 
    nelects = sum_ijk(rho_ijk * dV) #number of electrons in the simulation box
    
    
    eta_ijk = local potential in the voxel of ijk 
    avepot = sum_ijk(eta_ijk )/ngt  #average local potential in the simulation box
    
    """
    
    nargs  = len(sys.argv)
    if nargs < 3:
        print ('You have to specify a direction for integration and a reference file')
        print ('Syntax: chgave.py  -[x/y/z]  [CHGCAR/LOCPOT]')
        print ('e.g.: chgave.py  -z  CHGCAR')
        print ('Output: CHGCAR_Z')
        exit()
    
    #name ='CHGCAR_DIFF'
    #direction = 'Z'
    direction = sys.argv[1].strip('-').upper()
    name = sys.argv[2]
    
    print ('Reading %s file and analyzing volumetric data.....' %(name))
    chgobj = VaspChargeDensity(filename=name)
    rho = chgobj.chg[-1] #charge density 
    atoms =  chgobj.atoms[-1]
    vol = atoms.get_volume()
    cell = atoms.cell
    ngx, ngy, ngz = rho.shape
    Lx, Ly, Lz = np.sqrt(np.dot(cell, cell.T).diagonal()) 
    dx, dy, dz = Lx/ngx, Ly/ngy, Lz/ngz

    print ("Integrating bonding charge/local potential in planes perpendicular to the given direction")
    if direction == 'Z':
        df = dz
        dA = dx * dy
        ngt = ngx * ngy #number of grid points on xy plane  
        integ_rho = np.zeros(ngz)
        for idx in range(ngz):
            rho[:,:,idx] *= dA
            integ_rho[idx] = rho[:,:,idx].sum()
    elif direction == 'Y':
        df = dy
        dA = dx * dz
        ngt = ngx * ngz #number of grid points on xz plane  
        integ_rho = np.zeros(ngy) 
        for idx in range(ngy):
            rho[:,idx,:] *= dA
            integ_rho[idx] = rho[:,idx,:].sum()
    elif direction == 'X':
        df = dx
        dA = dy * dz
        ngt = ngy * ngz #number of grid points on yz plane
        integ_rho = np.zeros(ngx)
        for idx in range(ngx):
            rho[idx,:,:] *= dA
            integ_rho[idx] = rho[idx,:,:].sum()
    else:
        df = 0
        dV = dx * dy * dz
        ngt = ngx * ngy * ngz #number of grid points in the box
        rho[:,:,:] *= dV 
        integ_rho = rho.sum()
    
    print ("Averaging charge/local potential with respect to the planes perpendicular to the given direction" )    
    if 'LOCPOT' in name:
        tag1 = 'Distance(A)'
        tag2 = 'Potential(eV)' 
        integ_rho = (integ_rho * vol )/ (dA *  ngt)
    else:
        tag1 = 'Distance(A)'
        tag2 = 'rho(|e|/A^3)'
        integ_rho /= (dA * ngt)
  
    name = '%s_%s' %(name, direction)
    with open(name, 'w') as txt:
        print ('Writing averaged density/potential to %s file' %(name)) 
        txt.write("     %s    %s\n" %(tag1, tag2))
        for idx in range(len(integ_rho)):
            d = idx * df  
            txt.write("%16.8f %16.8f\n" % (d , integ_rho[idx]))
        print('Done !!')

