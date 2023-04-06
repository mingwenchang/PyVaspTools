#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Modified on Fri Sep 19 10:57:09 2019
Name: vasprun (assister for job submission )
Developer: Ming-Wen Chang
E-mail: ming.wen.c@gmail.com
"""

import sys, subprocess

#-----------------------------FOR CARTESIUS-----------------------------------#
#commands 
qsub = 'sbatch < ' #The command to submit a job
#qstat = 'squeue -j ' #The command to check status of jobs

#basic inforamtion  
module = 'vasp/5.3.5-avx'
account = 'ncteh231' 
cputype ='haswell'

#Available sources
avblpart = ['normal', 'short', 'fat']
avblcpus = ['12cpu', '24cpu', '48cpu', '96cpu', '196cpu', '240cpu']
avblhour = ['4HR', '8HR', '12HR', '24HR', '48HR', '72HR', '96HR', '120HR']
avblvasp = ['vasp.mpi', 'vasp.gamma']

pbstmp ="""#!/bin/sh
#SBATCH -J %(JOBNAME)s          #specify job name
#SBATCH -e %(JOBNAME)s.err      #specify standard error 
#SBATCH -o %(JOBNAME)s.out      #specify standard output 
#SBATCH -n %(NCPUS)s            #request slot range for parallel jobs
#SBATCH -C %(CPUTYPE)s          #bind job to haswell/ivy cpu
#SBATCH -t %(NHOUR)s:00:00      #Walltime
#SBATCH -A %(ACCOUNT)s          #Computational budgets account
#SBATCH -p %(PARTITION)s        #Request a specific partition: normal/fat/short

#Set vasp running environment
#(Import modules)
module load %(MODULE)s 

#Set mpirun and ncpus
#(Execute vasp)
time srun -n %(NCPUS)s %(VASP)s
wait

"""

pbsopt ="""#!/bin/sh
#SBATCH -J %(JOBNAME)s          #specify job name
#SBATCH -e %(JOBNAME)s.err      #specify standard error 
#SBATCH -o %(JOBNAME)s.out      #specify standard output 
#SBATCH -n %(NCPUS)s            #request slot range for parallel jobs
#SBATCH -C %(CPUTYPE)s          #bind job to haswell/ivy cpu
#SBATCH -t %(NHOUR)s:00:00      #Walltime
#SBATCH -A %(ACCOUNT)s          #Computational budgets account
#SBATCH -p %(PARTITION)s        #Request a specific partition: normal/fat/short

#Set vasp running environment
#(Import modules)
module load %(MODULE)s 

#Set mpirun and ncpus
#(Execute vasp)
sed -i 's/.*IBRION.*/IBRION = 2 /' INCAR
sed -i 's/.*EDIFF[^G].*/EDIFF = 1E-4 /' INCAR
time srun -n %(NCPUS)s %(VASP)s
wait

cp CONTCAR POSCAR
sed -i 's/.*ISTART.*/ISTART = 1 /' INCAR
sed -i 's/.*ICHARG.*/ICHARG = 1 /' INCAR
sed -i 's/.*IBRION.*/IBRION = 1 /' INCAR
sed -i 's/.*EDIFF[^G].*/EDIFF = 1E-6 /' INCAR
time srun -n %(NCPUS)s %(VASP)s
wait

"""
#-----------------------FOR CARTESIUS (END)-----------------------------------#


'''
#-----------------------FOR TUE LOCAL CLUSTER---------------------------------#

#commands 
qsub = 'qsub < ' #The command to submit a job

#basic inforamtion  
module = 'vasp/5.3.5'
account = 'unknown' 
cputype ='unknown'

#Available sources
avblpart = ['smk.q', 'all.q']
avblcpus = ['8cpu', '16cpu', '32cpu', '64cpu']
avblhour = ['00HR', '4HR', '8HR', '12HR', '24HR', '48HR', '72HR', '96HR', '120HR']
avblvasp = ['mpi', 'gammampi']

pbstmp="""#!/bin/sh
#$ -N %(JOBNAME)s             #specify job name  
#$ -e %(JOBNAME)s.err         #specify standard error 
#$ -o %(JOBNAME)s.out         #pecify standard output 
#$ -q %(PARTITION)s           #Request a specific partition
#$ -pe openmpi %(NCPUS)s      #request slot range for parallel jobs
#$ -l h_rt=%(NHOUR)s:00:00    #Walltime
#$ -m n
#$ -cwd
 
#Set vasp running environment
#(Import modules)
module purge
module load shared
module load intel/mkl/64
module load intel-mpi/64
module load %(MODULE)s/%(VASP)s

#Set mpirun and ncpus
#(Execute vasp)
time mpirun -np %(NCPUS)s vasp.real
wait

"""	
#-----------------------FOR TUE LOCAL CLUSTER (END)----------------------------#

'''

welcome = """
********************************************************************************
*                          _     _     _     _     _     _     _               *
*                         / \   / \   / \   / \   / \   / \   / \              *
*            Wellcome to ( S ) ( U ) ( R ) ( F ) ( H ) ( P ) ( C )             *
*                         \_/   \_/   \_/   \_/   \_/   \_/   \_/              *
*                                                                              *
********************************************************************************
*                                                                              *
*                        The assister for job submission                       *
*  Hello Users:                                                                *
*    Are you struggling and confused how to submite a job to HPC?              *
*    This script will help you do this. Just follows the guide below           *
*                                                                              *    
*  Syntaxs:                                                                    *
*    (1) vasprun                                                               *
*    (2) vasprun "jobname"                                                     *
*    (3) vasprun "jobname" "ncpu" "partition"  "walltime"  "version"           *
*                                                                              *
*  PS:                                                                         *
*    You can terminate this program at anytime,                                *                              
*    by presssing "CTRL" and "C" keys'                                         *
*                                                                              *    
********************************************************************************
*                                                                              *    
*  If you have further questions, please contact:                              *
*                                                                              *
*  Developer: Dr. Ming-Wen Chang                                               *
*  E-mail: ming.wen.c@gmail.com                                                *
*                                                                              *
********************************************************************************
"""

syntax = """                                                                   
    (1) vasprun                                                               
    (2) vasprun "jobname"                                                     
    (3) vasprun "jobname" "ncpu" "partition"  "walltime"  "version"           
"""


syntax2 = """                                                                     
     (1) vasprun -opt                                                              
     (2) vasprun -opt "jobname"                                                     
     (3) vasprun -opt "jobname" "ncpu" "partition"  "walltime"  "version"  
     (4) vasprun -mwc "jobname"  
     (5) vasprun -optmwc "jobname"  
"""
    
def assign_jobname():
    print ('To terminate this program: presssing "CTRL" and "C" keys')
    jobname=input('Please enter your job name (default: %s): ' %('job')).strip()
    print ('')
    if jobname == '':
        jobname = 'job'
    return jobname

def assign_ncpu(ncpu=None):
    if ncpu is None:
        print ('How many CPUs do you want to use?') 
        print ('Recommended CPUs:', ', '.join(avblcpus))
        ncpu=input('I want to use (default: %s): ' %(avblcpus[0]))
        print ('')
        
    ncpu = ncpu.lower().strip('cpu')
    if ncpu == '':
        ncpu = int(avblcpus[0].strip('cpu'))
    else:
        ncpu=int(ncpu)
    return ncpu	

def assign_partition(partition=None):
    if partition is None:
        print ('Which partition do you want to use?') 
        print ('Available partitions:', ', '.join(avblpart) )
        partition = input('I want to use (default: %s): ' %(avblpart[0]))
        print ('')
        
    partition = partition.lower().strip()
    if partition == '' or partition not in avblpart:
        partition = avblpart[0]		
    return partition

def assign_walltime(nhour=None):
    if nhour is None:
        print ('Please estimate how many hours the job will take to complete') 
        print ('For examples: 4HR, 8HR, 12HR, 24HR, 48HR, 72HR, 96HR, or 120HR(max)')
        nhour = input('I want to use (default: %s): ' %(avblhour[0]))
        print ('')
        
    nhour = nhour.lower().strip('hr')
    if partition == 'short':
        nhour = 1
    else:
        if nhour == '':
            nhour = 4
        else:
            nhour = int(nhour)
    return nhour	
	
def assign_version(version=None):
    if version is None:
        print ('Which vasp version do you want to use?') 
        print ('Available versions:', ', '.join(avblvasp))
        print ('PS: gamma version is only for use of kpoint=(1x1x1).')	
        version = input('I want to use (default: %s): ' %(avblvasp[0])).lower().strip()
        print ('')
        
    version = version.lower().strip()
    if version == '' or version not in avblvasp:
        version = avblvasp[0]	
    return version

if __name__ == "__main__":
    
    argv = sys.argv
    if len(argv) > 1:
        if argv[1] =='-opt':
            argv.pop(1)
            pbstmp = pbsopt
        elif argv[1] =='-mwc':
            argv.pop(1)
            jobname = argv[1]
            argv = ['vasprun', jobname, '96cpu', 'short', '1HR', 'vasp.gamma']
        elif argv[1] =='-optmwc':
            argv.pop(1)
            jobname = argv[1]
            argv = ['vasprun', jobname, '96cpu', 'short', '1HR', 'vasp.gamma']
            pbstmp = pbsopt
        #else:
        #    raise SyntaxError('\nSyntaxs: %s' %(syntax2))

    nargs = len(argv)
    if nargs == 1:
        print (welcome)
        jobname = assign_jobname()
        ncpu = assign_ncpu()
        partition = assign_partition()
        nhour = assign_walltime()
        version = assign_version()
    elif nargs == 2:
        print (welcome)
        jobname = argv[1]
        ncpu = assign_ncpu()
        partition = assign_partition()
        nhour = assign_walltime()
        version = assign_version()
    elif nargs == 6:
        try:
            jobname = argv[1]
            ncpu = assign_ncpu(argv[2])
            partition = assign_partition(argv[3])
            nhour = assign_walltime(argv[4])
            version = assign_version(argv[5])
        except:
            raise SyntaxError('\nSyntaxs: %s' %(syntax))

    else:
        raise SyntaxError('\nSyntaxs: %s' %(syntax))
 
    #Write pbs file 
    variables= {'JOBNAME':jobname,
                'NCPUS': ncpu,
                'PARTITION': partition,
                'NHOUR': nhour,
                'VASP':version,
                'CPUTYPE':cputype,
                'ACCOUNT': account,
                'MODULE': module}    
    
    print ('Generating a pbs script')  
    pbsfile = jobname + '.run' 
    with open(pbsfile, 'w') as txt:
        txt.write(pbstmp %(variables))
    
    subprocess.Popen('chmod +x %s' %(pbsfile), shell=True,stdout=subprocess.PIPE)
    subprocess.Popen(qsub + '%s'   %(pbsfile), shell=True,stdout=subprocess.PIPE)
    print ('Your job: %s has been submitted' %(jobname))
    
