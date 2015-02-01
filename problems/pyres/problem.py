import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy as np
import math
import os
import pytools
import random
import matplotlib.pyplot as plt
import pickle
from interruptible_pool import InterruptiblePool
import time

start_time = time.time()

def check_jumps(a, particles):
    com = particles[0]
    for i in range(1,rebound.get_N()):
        o = pytools.p2orbit(particles[i],com)

        if o.a*o.e > dr_thresh:
            print("Planet {} had a*e > {} AU\n".format(i,dr_thresh))
            return True

        if a[i] == 0.:
            a[i] = o.a
            
        if math.fabs(o.a - a[i]) > dr_thresh:
            print("Planet {} had semimajor axis jump by more than {} AU".format(i, dr_thresh))
            return True
        
        a[i] = o.a
        com = pytools.get_center_of_mass(com, particles[i]) 
    return

def integrate(args):    
    mthresh, taue, q = args
    with open('/Users/dtamayo/Desktop/starters/m_{0:.2e}_taue_{1:.2e}_{2:1d}.pickle'.format(mthresh,taue,q), 'r') as f:
        data = pickle.load(f)
        
    atrack = [0.,0.,0.,0.,0.,0.]
    rebound.reset()
    
    outputdelta=10.
    rebound.set_G(4.*math.pi**2)  
    tmax = 1.e2
    
    prevparticles = data[0]
    tauas = data[1]
    taues = data[2]
    tauis = data[3]
    
    for p in prevparticles: # data[0] is a list of the particles
        rebound.particle_add(p)

    N = rebound.get_N()  
    rebound.init_damping_forces()
    rebound.add_migration(tauas)
    rebound.add_e_damping(taues)
    rebound.add_i_damping(tauis)
    
    particles = rebound.particles_get()

    e = []
    t = []
    a = []
    P = []
    
    for q in range(N-1):
        e.append([])
        t.append([])   
        a.append([])
        P.append([])     
        
    last_t = -1e6
    N_output = 0

    while rebound.get_t()<tmax:
        _t = rebound.get_t()
        if _t - last_t > outputdelta:
            com = particles[0]
            for i in range(1,N):
                o = pytools.p2orbit(particles[i],com)
                t[i-1].append(_t)
                e[i-1].append(o.e)
                a[i-1].append(o.a)
                P[i-1].append(o.P)
                pytools.get_center_of_mass(com, particles[i])    
            N_output += 1
            last_t = _t
        rebound.step()
        breakflag = check_jumps(atrack,particles)
        if breakflag is True:
            break
    
    return (t,e,a,P)


dr_thresh = 5.        
taue=5.e5
mthresh=2.5e-4
n_restarts = 3

massmin = 5
massmax = 1000
Nmass = 10

mearth = 3.0024584e-6 # in solar masses

if Nmass == 1:
    deltamass=0
else:
    deltamass = (np.log(massmax)-np.log(massmin))/(Nmass-1)

masses = []
for j in range(0,Nmass):
    logmass = np.log(massmin) + j*deltamass
    masses.append(pow(np.e,logmass)*mearth)
    
'''args = []
for q in range(n_restarts):
    args.append((mthresh,taue,q))

pool = InterruptiblePool()
vars = pool.map(integrate,args)'''

t,e,a,P = integrate((masses[3],1.e7,3))

print(np.asarray(t).shape)

'''t0 = vars[0][0][0]
t1 = vars[1][0][0]
t2 = vars[2][0][0]

e0 = vars[0][1][2] # 0th shadow run, take the e from (t,e) and particle index 2 (middle one, most unstable)
e1 = vars[1][1][2] # 1st shadow run...
e2 = vars[2][1][2]'''

start = 0
end = max([len(t[i]) for i in range(5)])
t0 = t[0]

Pratio1 = [P[4][q]/P[3][q] for q in range(min(len(P[4]),len(P[3])))]
Pratio2 = [P[3][q]/P[2][q] for q in range(min(len(P[3]),len(P[2])))]

fig,ax = plt.subplots(3)

for q in range(2,5):
    ax[0].plot(t0[start:end],e[q][start:end], '.')
    
ax[1].plot(t0[start:end],Pratio1, '.')
ax[1].plot(t0[start:end],Pratio2, '.')

for q in range(5):
    ax[2].plot(t0[start:end],a[q][start:end], '.')


print("Took {0} seconds".format(time.time() - start_time))
plt.show()