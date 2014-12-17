'''
Created on Oct 9, 2014

@author: dtamayo
'''
import numpy as np
from subprocess import call

mass = 5.11e-2
numiter=1
numcores=1

delta = np.linspace(1,2,3)

itspercore = numiter/numcores
print("%s iterations on each core"%(itspercore))
print("%s runs per core"%(itspercore*len(delta)))

of = open("masterbatch", "w")
of.write("#!/bin/bash\n")

for q in range(numcores):
    of.write("./batch%d\n"%(q))
of.close()
call("chmod u=rwx masterbatch", shell=True)

for q in range(numcores):
    of = open("batch%d"%(q), "w")
    of.write("#!/bin/bash\n")
    of.write("#PBS -l nodes=1:ppn=1\n")
    of.write("#PBS -q workq\n")
    of.write("#PBS -r n\n")
    of.write("#PBS -l walltime=47:00:00\n")
    of.write("#PBS -N batch%d\n"%(q))
    of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
    of.write("cd $PBS_O_WORKDIR\n")
    of.write("./script%d\n"%(q))
    of.close()
    call("chmod u=rwx batch%d"%(q), shell=True)

    of = open("script%d"%(q), "w")
    of.write("#!/bin/bash\n")
    
    for j in delta:
        of.write("./nbody --mass=%0.2e --delta=%0.2e --it=%d\n"%(mass,j,q))
    of.close()
    call("chmod u=rwx script%d"%(q), shell=True)
