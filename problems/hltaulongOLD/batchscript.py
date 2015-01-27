'''
Created on Oct 9, 2014

@author: dtamayo
'''

import numpy as np
from subprocess import call

tauis= [-1] # use -1 to use same as taue

numiter=20
numcores=2

massmin = 5
massmax = 1000
Nmass = 10

tauemin = 1e4
tauemax = 1e8
Ntaues = 5

k = 100
starmass = 0.55
sigmae = 0.01

if Nmass == 1:
    deltamass=0
else:
    deltamass = (np.log(massmax)-np.log(massmin))/(Nmass-1)
    
if Ntaues == 1:
    deltataue=0
else:
    deltataue = (np.log(tauemax) - np.log(tauemin))/(Ntaues-1)

mearth = 0.00314635457 # in mjup

itspercore = numiter/numcores
print("%s iterations on each core"%(itspercore))
print("%s runs per core"%(itspercore*Nmass*Ntaues))


of = open("masterbatch", "w")
of.write("#!/bin/bash\n")

for taui in tauis:
    for q in range(numcores):
        of.write("qsub taui=%0.2e_batch%d\n"%(taui,q))

of.close()
call("chmod u=rwx masterbatch", shell=True)

for taui in tauis:
    for q in range(numcores):
        of = open("taui=%0.2e_batch%d"%(taui,q), "w")
        of.write("#!/bin/bash\n")
        of.write("#PBS -l nodes=1:ppn=1\n")
        of.write("#PBS -q workq\n")
        of.write("#PBS -r n\n")
        of.write("#PBS -l walltime=18:00:00\n")
        of.write("#PBS -N taui=%0.2e_%d\n"%(taui,q))
        of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
        of.write("cd $PBS_O_WORKDIR\n")
        of.write("./taui=%0.2e_%d\n"%(taui,q))
        of.close()
        call("chmod u=rwx taui=%0.2e_batch%d"%(taui,q), shell=True)


    for q in range(numcores):
        of = open("taui=%0.2e_%d"%(taui,q), "w")
        of.write("#!/bin/bash\n")
    
        for j in range(0,Nmass):
            for i in range(Ntaues):
                for k in range(q*itspercore,(q+1)*itspercore):
                    logmass = np.log(massmin) + j*deltamass
                    mass = pow(np.e,logmass)
                    logtaue = np.log(tauemin) + i*deltataue
                    taue = pow(np.e,logtaue)
                    if taui == -1:
                        of.write("./nbody --mass=%0.2e --taue=%0.e --sigmae=%0.e --it=%d\n"%(mass*mearth,taue,sigmae,k))
                    else:
                        of.write("./nbody --mass=%0.2e --taue=%0.e --taui=%0.e --sigmae=%0.e --it=%d\n"%(mass*mearth,taue,taui,sigmae,k))
                
        of.close()
        call("chmod u=rwx taui=%0.2e_%d"%(taui,q), shell=True)