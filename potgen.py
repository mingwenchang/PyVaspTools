#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Modified on Fri Sep 19 10:57:09 2019
Name: potgen (for generating POTCAR )
Developer: Ming-Wen Chang
E-mail: ming.wen.c@gmail.com
"""
import os, sys, subprocess
import vasp_modules.vasp_io2 as vio


root='/home/mwchang/pkg/pot' #path to potcar dirs
pseudos=('pbe', 'pbe.52', 'lda', 'lda.52', 'gga', 'ugga', 'ulda')   

syntax = """
    (1) potgen
    (2) potgen POSCAR
    (3) potgen [pbe/pbe.52/lda/lda.52/gga/ugga/ulda] POSCAR   
    (4) potgen [pbe/pbe.52/lda/lda.52/gga/ugga/ulda] [elements]
"""
def assign_pseudo():
    print ('Which pseudopotential do you want to use?')
    print ('\n[1]: paw_PBE\n[2]: paw_GGA\n[3]: paw_LDA\n[4]: ussp_gga\n[5]: ussp_lda\n')
    print ('or press CTRL + C keys to leave this script')
    select = input('Please enter a number (default: [1]: paw_PBE): ').strip()
   
    if select == '':
        pseudo = 'pbe'
    elif select == '1':
        pseudo = 'pbe'
    elif select == '2':
        pseudo = 'gga'
    elif select == '3':
       pseudo = 'lda'
    elif select == '4':
       pseudo = 'ugga'        
    elif select == '5':
       pseudo = 'ulda'  
    else:
        raise IOError('Please enter a correct number!!!')

    return pseudo

def check_pseudo(pseudo):
    if pseudo not in pseudos:
        raise IOError('Please enter a correct pseudopotential name!!!')

def assign_elements():
    print ('Please enter atomic species. e.g. H, C, O, Ce')
    print ('If you want to use other special pseudopotential, for example: a harder pseudopotential,') 
    print ('p or s semi-core states treated as valence states, etc.')
    print ('just put an extension like _pv, _sv and _d after the element name. e.g. H, V_sv, Ru_pv')

    elements = input('Atomic species: ').strip().split()
    
    if len(elements) != 0:
        elements = [element.strip(',') for element in elements]
    else:
        raise IOError('Plese enter atomic species!!!')
   
    return elements

def get_elements(file='POSCAR'):
    poscar = vio.read_poscar(file)
    return poscar.atomtypes

#Main Program

if __name__ == "__main__":
    if os.path.exists('POTCAR'):
        os.remove('POTCAR')
    
    argulen=len(sys.argv)
    if argulen == 1:
        pseudo = assign_pseudo()
        elements = assign_elements()
    elif argulen == 2 and sys.argv[1] == 'POSCAR':
        pseudo = 'pbe'
        elements= get_elements('POSCAR')
    elif argulen == 3 and sys.argv[1].lower() in pseudos and sys.argv[2] == 'POSCAR':
        pseudo = sys.argv[1].lower().strip()
        elements= get_elements('POSCAR')
    elif argulen > 3 and sys.argv[1].lower() in pseudos:
        pseudo = sys.argv[1].lower().strip()
        elements = [element.strip(',') for element in sys.argv[2:]]
    else:
        raise SyntaxError('\nSyntaxs: %s' %(syntax))
        
    for element in elements:
        potcar ='%s/%s/%s/POTCAR' %(root, pseudo, element)
        cmd='cat %s >> POTCAR' %(potcar)
        
        if os.path.exists(potcar):
            subprocess.call(cmd, shell=True)
        else:
            raise IOError('The POTCAR of %s element does not exist!!' %potcar)
            
        
