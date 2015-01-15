'''
Created on Oct 9, 2014

@author: dtamayo
'''
import numpy as np
from subprocess import call

delta = 4.5
numiter=20
numcores=20

massmin = 5
massmax = 1000
Nmass = 50
if Nmass == 1:
    deltamass=0
else:
    deltamass = (np.log(massmax)-np.log(massmin))/(Nmass-1)

mearth = 0.00314635457 # in mjup

itspercore = numiter/numcores
print("%s iterations on each core"%(itspercore))
print("%s runs per core"%(itspercore*Nmass))

of = open("mmasterbatch", "w")
of.write("#!/bin/bash\n")

for q in range(numcores):
    of.write("qsub mbatch%d\n"%(q))
of.close()
call("chmod u=rwx mmasterbatch", shell=True)

for q in range(numcores):
    of = open("mbatch%d"%(q), "w")
    of.write("#!/bin/bash\n")
    of.write("#PBS -l nodes=1:ppn=1\n")
    of.write("#PBS -q workq\n")
    of.write("#PBS -r n\n")
    of.write("#PBS -l walltime=47:00:00\n")
    of.write("#PBS -N mbatch%d\n"%(q))
    of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
    of.write("cd $PBS_O_WORKDIR\n")
    of.write("./mscript%d\n"%(q))
    of.close()
    call("chmod u=rwx mbatch%d"%(q), shell=True)

    of = open("mscript%d"%(q), "w")
    of.write("#!/bin/bash\n")
    
    for j in range(0,Nmass):
        logmass = np.log(massmin) + j*deltamass
        mass = pow(np.e,logmass)
        of.write("./nbody --mass=%0.2e --delta=%0.2e --it=%d\n"%(mass*mearth,delta,q))
    of.close()
    call("chmod u=rwx mscript%d"%(q), shell=True)
