'''
Created on Oct 9, 2014

@author: dtamayo
'''
import numpy as np
from subprocess import call

daOvera = 0.45
a3 = 66
numiter=20

starmass = 0.55 #solar masses
msun = 1048. #jupiter masses.

delta = np.linspace(1,6,51)
q = np.linspace(1,10,19)
numcores=len(q)

of = open("masterbatch", "w")
of.write("#!/bin/bash\n")

for m in range(len(q)):
    of.write("qsub a3_%s_da_%s_q%s\n"%(a3, daOvera, q[m]))
of.close()
call("chmod u=rwx masterbatch", shell=True)

for m in range(len(q)):
    of = open("a3_%s_da_%s_q%s"%(a3, daOvera, q[m]), "w")
    of.write("#!/bin/bash\n")
    of.write("#PBS -l nodes=1:ppn=1\n")
    of.write("#PBS -q workq\n")
    of.write("#PBS -r n\n")
    of.write("#PBS -l walltime=47:55:00\n")
    of.write("#PBS -N a3_%s_da_%s_m%s\n"%(a3, daOvera, q[m]))
    of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
    of.write("cd $PBS_O_WORKDIR\n")
    of.write("./scripta3_%s_da_%s_q%s\n"%(a3, daOvera, q[m]))
    of.close()
    call("chmod u=rwx a3_%s_da_%s_q%s"%(a3, daOvera, q[m]), shell=True)

    of = open("scripta3_%s_da_%s_q%s"%(a3, daOvera, q[m]), "w")
    of.write("#!/bin/bash\n")
    
    for j in delta:
        for k in range(numiter):
            mass = starmass*msun*1.5*(daOvera/j)**3
            of.write("./nbody --a3=%0.2e --daOvera=%0.2e --q=%s --mass=%0.2e --it=%d\n"%(a3,daOvera,q[m],mass,k))
    of.close()
    call("chmod u=rwx scripta3_%s_da_%s_q%s"%(a3, daOvera, q[m]), shell=True)
