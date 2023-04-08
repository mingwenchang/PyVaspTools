#!/usr/bin/python
# coding=UTF-8

import os, sys, subprocess

def ShowWelcome():
    welcome = (
    '********************************************************************************',
    '*                          _     _     _     _     _     _     _               *',
    '*                         / \   / \   / \   / \   / \   / \   / \              *',
    '*            Wellcome to ( T ) ( U ) ( / ) ( E ) ( H ) ( P ) ( C )             *',
    '*                         \_/   \_/   \_/   \_/   \_/   \_/   \_/              *',
    '*                                                                              *',
    '********************************************************************************',
    '*                                                                              *',
    '*                  Submission assister for VASP calculations                   *',
    '*                                                                              *',
    '*  Syntaxs                                                                     *',
    '*                                                                              *',
    '*    (1) vasprun                                                               *',
    '*    (2) vasprun "jobname"                                                     *',
    '*    (3) vasprun "jonname" -np "ncpu"                                          *',
    '*    (4) vasprun "jonname" -np "ncpu" -q "qname"                               *',
    '*    (5) vasprun -all "jobname" "ncpu" "qname" "version" "events"              *',
    '*                                                                              *',
    '*  Hello Users:                                                                *',
    '*    Are you struggling in submitting vasp job to TU/E HPC?                    *',
    '*    Take it easy and keep relax. This script will help                        *',
    '*    you submit vasp calculatins to TUE local cluster.                         *',
    '*    Please follow the guide below                                             *',
    '********************************************************************************',
    '*  If you have further question, please contact:                               *',
    '*                                                                              *',
    '*    Dr. Ming-Wen Chang                                                        *',
    '*    E-Mail: m.chang@tue.nl;  STW 3.43                                         *',
    '*                                                                              *',
    '********************************************************************************')

    for line in welcome:
        print (line)

def GetJobName():
    print ('Please enter your job name,')
    print ('If you want to leave, easily press "CTRL" and "C" keys')
    jobname=raw_input('Job name: ').strip()
	
    if jobname == '':
        jobname = 'vasp'

    return jobname

def AssignCPU():
    print ('How many CPUs do you want to use?') 
    print ('These are the recommended cpu numbers on TU/E HPC:')
    print ('8cpu', '16cpu', '32cpu', '48cpu', '64cpu', '96cpu', '128cpu', '192cpu')
    ncpu=raw_input('I want to use (default: 8cpu): ').lower().strip('cpu')
	
    if ncpu == '':
        ncpu=8
    else:
        ncpu=int(ncpu)

    return ncpu	

def AssignQueue():
    avblq=('all.q', 'smk.q')
    queue=raw_input('Which queue do you want to use? all.q or smk.q? (defalut: all.q): ').lower().strip()
	
    if queue == '' or queue not in avblq:
        queue='all.q'	
		
    return queue
	
def AssignVersion():
    avblver=('mpi', 'gammampi')
    print ('Whic vasp version do you want to use? mpi or gammampi?') 
    print ('PS: gammampi is twice as fast, but only for use in gamma point.')	
    version=raw_input('I want to use (defalut: mpi): ').lower().strip()
	
    if version == '':
        version = avblver[0] #mpi
    elif version[0] == 'm': 
        version = avblver[0] #mpi
    elif version[0] == 'g':
        version = avblver[1] #gammampi
    else:
        version = avblver[0] #mpi

    return version

def ReceiveNotice():

    print ('Do you want to receive notifications by mail when the job begins, ends and aborts?')
    ans=raw_input('yes or no (default: No): ').lower().strip()
	
    if ans == 'yes' or ans == 'y':
        argument='beas'
    else:
        argument='n'

    return argument		
	
def WritePBScript():

    values={'JOBNAME':jobname,
            'QUEUE':queue,
            'NCPU':ncpu,
            'VER':version,
            'ARG':argument}

    pbstmp ="""#!/bin/sh
#$ -N %(JOBNAME)s         
#$ -e %(JOBNAME)s.err
#$ -o %(JOBNAME)s.out
#$ -q %(QUEUE)s
#$ -m %(ARG)s
#$ -pe openmpi %(NCPU)s
#$ -cwd
 
#Import modules 
module purge
module load shared
module load intel/mkl/64
module load intel-mpi/64
module load vasp/5.3.5/%(VER)s

#mpirun 
time mpirun -np %(NCPU)s vasp.real
wait
	
"""
    scrname=jobname+'.run'
    pbscr = open(scrname, 'w')
    print >> pbscr, pbstmp %values
    pbscr.close()
    return scrname

#Main Program
#vaspinputs=('INCAR', 'POSCAR', 'POTCAR', 'KPOINTS')
#missings=[vaspinput for vaspinput in vaspinputs if not os.path.exists(vaspinput)]
#if len(missings) != 0:
#    raise IOError('%s files are not found.' %(missings))

argulen=len(sys.argv)
if argulen == 1:
    ShowWelcome()
    jobname = GetJobName()
    ncpu = AssignCPU()
    queue = AssignQueue()
    version = AssignVersion()
    argument = ReceiveNotice()
elif argulen == 2:
    jobname = sys.argv[1]
    ncpu = AssignCPU()
    queue = AssignQueue()
    version = AssignVersion()
    argument = ReceiveNotice()
elif argulen == 4:
    if sys.argv[2].lower() == '-np':
        jobname = sys.argv[1]
        ncpu = sys.argv[3].lower().strip('cpu')
        queue = 'all.q' 
        version = 'mpi' 
        argument = 'n' 
    else:
        raise SyntaxError('invalid syntax')
elif argulen == 6:
    if sys.argv[2].lower() == '-np' and sys.argv[4].lower() == '-q':
        jobname = sys.argv[1]
        ncpu = sys.argv[3].lower().strip('cpu')
        queue = sys.argv[5].lower().strip()
        version = 'mpi' 
        argument = 'n'
    else: 
        raise SyntaxError('invalid syntax')
elif argulen == 7:
    if sys.argv[1].lower() == '-all':
        jobname=sys.argv[2]
        ncpu = sys.argv[3].lower().strip('cpu')
        queue = sys.argv[4].lower().strip() 
        version = sys.argv[5].lower().strip('"').strip("'") 
        argument = sys.argv[6].lower().strip('"').strip("'")
    else: 
        raise SyntaxError('invalid syntax')
else:
    raise SyntaxError('invalid syntax')

pbsfile = WritePBScript()
subprocess.Popen('chmod +x %s' %(pbsfile), shell=True,stdout=subprocess.PIPE)
subprocess.Popen('qsub < %s'   %(pbsfile), shell=True,stdout=subprocess.PIPE)


    
