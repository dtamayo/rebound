'''
Created on Oct 9, 2014

@author: dtamayo
'''

import numpy as np
from subprocess import call

j = 2.
nRH = 3.7
numDeltas = 101 # choose odd number to include delta = 0

daOveraRes = ((j+1)/j)**(2./3.) - 1.
Deltacrit = 10*(1.5*(daOveraRes/nRH)**3)/daOveraRes
Deltas = np.linspace(-3*Deltacrit, 3*Deltacrit, numDeltas, endpoint=True)

for Delta in Deltas:
    with open("Delta_{0:.5f}".format(Delta), "w") as of:
        of.write("#!/bin/bash\n")
        of.write("#PBS -l nodes=1:ppn=8\n")
        of.write("#PBS -q workq\n")
        of.write("#PBS -r n\n")
        of.write("#PBS -l walltime=1:00:00\n")
        of.write("#PBS -N Deltacrit\n")
        of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
        of.write("cd $PBS_O_WORKDIR\n")
        of.write("python problem.py -D {0:.5e} -n {1:.3e} -j {2:.2e}\n".format(Delta,nRH,j))
    
        call("chmod u=rwx Delta_{0:.5f}".format(Delta), shell=True)

with open("sunnyscript", "w") as of:
    of.write("#!/bin/bash\n")
    for Delta in Deltas:
        of.write("qsub Delta_{0:.5f}\n".format(Delta))
        
    call("chmod u=rwx sunnyscript", shell=True)
