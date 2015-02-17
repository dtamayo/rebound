'''
Created on Oct 9, 2014

@author: dtamayo
'''

import numpy as np
from subprocess import call

def runef(Delta):
    j = 2.
    numefs = 100

    efs = np.linspace(0.001, 0.1, numefs, endpoint=True)
    
    for ef in efs:
        with open("ef_{0:.5e}".format(ef), "w") as of:
            of.write("#!/bin/bash\n")
            of.write("#PBS -l nodes=1:ppn=8\n")
            of.write("#PBS -q workq\n")
            of.write("#PBS -r n\n")
            of.write("#PBS -l walltime=1:00:00\n")
            of.write("#PBS -N Deltacrit\n")
            of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
            of.write("cd $PBS_O_WORKDIR\n")
            of.write("python problem.py -D {0:.5e} -e {1:.5e} -j {2:.2e}\n".format(Delta,ef,j))
    
            call("chmod u=rwx ef_{0:.5e}".format(ef), shell=True)

    with open("sunnyscript", "w") as of:
        of.write("#!/bin/bash\n")
        for ef in efs:
            of.write("qsub ef_{0:.5e}\n".format(ef))
        
        call("chmod u=rwx sunnyscript", shell=True)

    call("./sunnyscript", shell=True)
    
if __name__ == '__main__':
    Delta = 0.01
    runef(Delta)