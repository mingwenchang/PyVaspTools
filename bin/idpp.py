#!/usr/bin/env python3
# coding=UTF-8
"""
Created on Mon Sep 24 13:42:41 2018

@author: 20180239
"""
import os, sys
from ase.neb import NEB
from ase.io import read, write


initial =  sys.argv[1] 
final   =  sys.argv[2] 
nimages =  int(sys.argv[3])  

IS = read(initial, format='vasp')
FS = read(final, format='vasp')

#create a list of images for interpolation
images = [IS]

for i in range(nimages):
    images.append(IS.copy())
    
images.append(FS)


neb = NEB(images)    
neb.interpolate('idpp')


for i in range(nimages+2):
    dirname = str(i).zfill(2)
    os.mkdir(dirname)
    poscar = '%s/%s' %(dirname, 'POSCAR')
    write(poscar, neb.images[i], format='vasp')
    
