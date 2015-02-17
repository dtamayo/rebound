'''
Created on Oct 9, 2014

@author: dtamayo
'''

import numpy as np
from subprocess import call

j = 2.
numnRH = 30 
Delta = 0.05

nRHs = np.linspace(3.5, 4.5, numnRH, endpoint=True)

for nRH in nRHs:
    with open("nRH_{0:.3f}".format(nRH), "w") as of:
        of.write("#!/bin/bash\n")
        of.write("#PBS -l nodes=1:ppn=8\n")
        of.write("#PBS -q workq\n")
        of.write("#PBS -r n\n")
        of.write("#PBS -l walltime=1:00:00\n")
        of.write("#PBS -N Deltacrit\n")
        of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
        of.write("cd $PBS_O_WORKDIR\n")
        of.write("python problem.py -D {0:.3e} -n {1:.3e} -j {2:.2e}\n".format(Delta,nRH,j))
    
        call("chmod u=rwx nRH_{0:.3f}".format(nRH), shell=True)

with open("sunnyscript", "w") as of:
    of.write("#!/bin/bash\n")
    for nRH in nRHs:
        of.write("qsub nRH_{0:.3f}\n".format(nRH))
        
    call("chmod u=rwx sunnyscript", shell=True)
