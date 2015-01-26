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
    with open('m_{:.2e}_taue_{:.2e}_{:2d}.bin'.format(mthresh,taue,q), 'rb') as f:
        data = pickle.load(f)
        
    a = [0.,0.,0.,0.]
    rebound.reset()
    
    outputdelta=10.
    rebound.set_G(4.*math.pi**2)  
    tmax = 1.e6
    
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
    
    for q in range(N-1):
        e.append([])
        t.append([])        
        
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
                pytools.get_center_of_mass(com, particles[i])    
            N_output += 1
            last_t = _t
        rebound.step()
        breakflag = check_jumps(a,particles)
        if breakflag is True:
            break
    
    return (t,e)


dr_thresh = 5.        
taue=5.e5
mthresh=2.5e-4
n_restarts = 3

args = []
for q in range(n_restarts):
    args.append((mthresh,taue,q))

pool = InterruptiblePool()
vars = pool.map(integrate,args)

t0 = vars[0][0][0]
t1 = vars[1][0][0]
t2 = vars[2][0][0]

e0 = vars[0][1][2] # 0th shadow run, take the e from (t,e) and particle index 2 (middle one, most unstable)
e1 = vars[1][1][2] # 1st shadow run...
e2 = vars[2][1][2]

start = 0
end0 = len(t0)
end1 = len(t1)
end2 = len(t2)

fig,ax = plt.subplots(3)
ax[0].plot(t0[start:end0],e0[start:end0], '.')
ax[1].plot(t1[start:end1],e1[start:end1], '.')
ax[2].plot(t2[start:end2],e2[start:end2], '.')


print("Took {0} seconds".format(time.time() - start_time))
plt.show()