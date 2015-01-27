'''
Created on Oct 9, 2014

@author: dtamayo
'''

import numpy as np
from subprocess import call

numiter=20
numcores=2

massmin = 5
massmax = 1000
Nmass = 10

mearth = 3.0024584e-6 # in solar masses

if Nmass == 1:
    deltamass=0
else:
    deltamass = (np.log(massmax)-np.log(massmin))/(Nmass-1)

for j in range(0,Nmass):
    logmass = np.log(massmin) + j*deltamass
    mass = pow(np.e,logmass)
    with open("m_{0:.1f}".format(mass), "w") as of:
        of.write("#!/bin/bash\n")
        of.write("#PBS -l nodes=1:ppn=4\n")
        of.write("#PBS -q workq\n")
        of.write("#PBS -r n\n")
        of.write("#PBS -l walltime=2:00:00\n")
        of.write("#PBS -N pygrowparallel\n")
        of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
        of.write("cd $PBS_O_WORKDIR\n")
        of.write("qsub python problem.py --mass={0:.3e}\n".format(mass*mearth))
    
        call("chmod u=rwx m_{0:.1f}".format(mass), shell=True)

with open("sunnyscript", "w") as of:
    of.write("#!/bin/bash\n")
    for j in range(0,Nmass):
        logmass = np.log(massmin) + j*deltamass
        mass = pow(np.e,logmass)
        of.write("./m_{0:.1f}\n".format(mass))
        
    call("chmod u=rwx sunnyscript", shell=True)