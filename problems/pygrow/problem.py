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

import time
start_time = time.time()

def make_restarts():
    #make the first one the original simulation
    ps = []
    for q in range(N):
        ps.append(particles[q])
        if q == 1:
            print(particles[q].x)
    with open('m_{:.2e}_taue_{:.2e}_{:2d}.bin'.format(mthresh,taues[1],0), 'wb') as f:
        pickle.dump((ps,tauas,taues,tauis),f)
    #now perturb n_restarts-1 times
    for i in range(n_restarts-1):
        com = particles[0] # just use star as com.  We're then jiggling the heliocentric f (diff between heliocentric and jacobi shouldn't matter)
        ps = []
        ps.append(particles[0])
        for q in range(1,N):
            o = pytools.p2orbit(particles[q],com)
            o.f += np.random.uniform(-delta,delta)
            p = pytools.kepler_particle(particles[q].m,com,o.a,o.f,o.e,o.omega,o.inc,o.Omega)
            ps.append(p)
            if q == 1:
                print(p.x)
        with open('m_{:.2e}_taue_{:.2e}_{:2d}.bin'.format(mthresh,taues[1],i+1), 'wb') as f:
            pickle.dump((ps,tauas,taues,tauis),f)
    
    
mthresh = 3.e-4 #solar masses
n_restarts = 3
delta = 2e-2

outputdelta=20.
rebound.set_G(4.*math.pi**2)  

tmax = 2.e7

with open('starter.pickle', 'rb') as f:
    data = pickle.load(f)

prevparticles = data[0]
tauas = data[1]
taues = data[2]
tauis = data[3]

for p in prevparticles: # data[0] is a list of the particles
    rebound.particle_add(p)

rebound.add_migration(tauas)
rebound.add_e_damping(taues)
rebound.add_i_damping(tauis)

particles = rebound.particles_get()

N = rebound.get_N()

Om = []
a = []
e = []
w = []
ml = []
t = []
P = []
mass = []

for q in range(N-1):
    Om.append([])
    a.append([])
    e.append([])
    w.append([])
    ml.append([])
    t.append([])
    P.append([])
    mass.append([])

last_t = -1 # to keep track of output
N_output = 0

while rebound.get_t()<tmax:
    _t = rebound.get_t()
    if _t - last_t > outputdelta:
        com = particles[0]
        #dt = 0 if last_t == -1e6 else _t - last_t
        for i in range(1,N):
            o = pytools.p2orbit(particles[i],com)#particles[0])
            t[i-1].append(_t)
            a[i-1].append(o.a)
            e[i-1].append(o.e)
            Om[i-1].append(o.Omega)
            w[i-1].append(o.omega)
            ml[i-1].append(o.l)
            P[i-1].append(o.P)
            mass[i-1].append(particles[i].m)
            
            if i==N-1 and particles[i].m < mthresh:
                tlib = 0.078*(particles[i].m/0.55)**(-2./3.)*o.P
                deltaM = particles[i].m * (_t - last_t) / 10. / tlib
                
            com = pytools.get_center_of_mass(com, particles[i])    
        N_output += 1
        last_t = _t
        for i in range(1,N):
            particles[i].m += deltaM
    
        if particles[3].m > mthresh:
            make_restarts()
            break

    rebound.step()

print("time at end = {}".format(_t))
print("mass at end = {}".format(particles[3].m))
print("N planets = {}".format(rebound.get_N()))
print(N_output)
if (mthresh > particles[3].m):
    print("Didn't reach specified mthresh")
    
phi = pytools.resarg(4,-3,-1,0, ml, Om, w, 1, 2)
q = [e[1][j]*math.cos(phi[j]*math.pi/180.) for j in range(len(e[0]))]
p = [e[1][j]*math.sin(phi[j]*math.pi/180.) for j in range(len(e[0]))]
Pratio1 = [P[1][j]/P[0][j] for j in range(len(P[0]))]
Pratio2 = [P[2][j]/P[1][j] for j in range(len(P[0]))]

start = 0
end = len(t[0])

fig,ax = plt.subplots(4)
ax[0].plot(t[0][start:end],Pratio1[start:end], ',')
ax[0].plot(t[0][start:end],Pratio2[start:end], ',')
ax[1].plot(t[0][start:end],phi[start:end], ',')
ax[2].plot(t[0][start:end], mass[1][start:end], ',')
ax[3].plot(t[0][start:end], e[0][start:end], ',')
ax[3].plot(t[0][start:end], e[1][start:end], ',')
ax[3].plot(t[0][start:end], e[2][start:end], ',')

print("Took {0} seconds".format(time.time() - start_time))
plt.show()