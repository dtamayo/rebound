import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy as np
import math
import os
import pytools
import random
import matplotlib.pyplot as plt

import time
start_time = time.time()

#rebound.set_G(4.*math.pi**2)  

starmass = 0.55     # in solar masses
N = 3               # including central star    
a = [0.,42.,53.]    # AU
mass = 1.e-6           # solar masses
e = 1.e-8
inc = 1.e-8
#taue=1.e4
#taua=-taue*1000
outputdelta=1.

tmax = 1.e2
rebound.particle_add( rebound.Particle(m=1.) )                  # Test particle
rebound.particle_add( rebound.Particle(m=1e-3,x=1.,vy=1.) )     # Planet

# Add particles
#star = rebound.Particle(m=starmass,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
#rebound.particle_add(star)                  # Star at origin (zeros by default for pos & vel)
'''
#for j in range(1,N):
    rebound.particle_add(pytools.init_planet(star,a=a[j],e=e,omega=random.uniform(0,2*math.pi),
                                             f=random.uniform(0,2*math.pi),inc=inc,
                                             Omega=random.uniform(0,2*math.pi), mass=mass))
'''
#rebound.init_damping_forces()
#rebound.add_migration([0.,0.,0.])
#rebound.add_e_damping([0.,0.,0.])
#rebound.add_i_damping([0.,0.,0.])

rebound.move_to_center_of_momentum()
particles = rebound.particles_get()

# timestep counter
steps = 0 

last_t = -1e6

xs = np.zeros((rebound.get_N(),np.round(tmax/outputdelta)))
ys = np.zeros((rebound.get_N(),np.round(tmax/outputdelta)))
zs = np.zeros((rebound.get_N(),np.round(tmax/outputdelta)))
#ts = numpy.zeros(numpy.round(tmax/outputdelta))
N_output = 0

Om = [[],[]]
a = [[],[]]
w = [[],[]]
ml = [[],[]]
t = [[],[]]
P = [[],[]]
o = rebound.Orbit()

_xs = []
_ys = []
_zs = []
while rebound.get_t()<tmax:
    _t = rebound.get_t()
    if _t - last_t > outputdelta:
        for i in range(1,rebound.get_N()):
            o = pytools.p2orbit(particles[i],particles[0])
            t[i-1].append(_t)
            a[i-1].append(o.a)
            Om[i-1].append(o.Omega)
            w[i-1].append(o.omega)
            ml[i-1].append(o.l)
            P[i-1].append(o.P)
            xs[i,N_output] = particles[i].x
            ys[i,N_output] = particles[i].y
            zs[i,N_output] = particles[i].z
            #print(particles[i].x)
            
        N_output += 1
        last_t = _t
    _xs.append(particles[1].x)
    _ys.append(particles[1].y)
    _zs.append(particles[1].z)
    rebound.step()
'''
while rebound.get_t()<tmax:
    _t = rebound.get_t()
    if _t - last_t > outputdelta:
        for i in range(1,rebound.get_N()):
            o = pytools.p2orbit(particles[i],particles[0])
            t[i-1].append(_t)
            a[i-1].append(o.a)
            Om[i-1].append(o.Omega)
            w[i-1].append(o.omega)
            ml[i-1].append(o.l)
            P[i-1].append(o.P)
            xs[i,N_output] = particles[i].x
            ys[i,N_output] = particles[i].y
            zs[i,N_output] = particles[i].z
            #print(particles[i].x)
            
        N_output += 1
        last_t = _t
    rebound.step()'''
'''
N = 2 #number of planets
numavg = 200 # number of pts to use from the end of the simulation to average P2/P1 and resonant angles
endindex = len(t[0]) #last index in time array

pairs = []
#do this for adjacent pairs
for v in range(1,N):
    pairs.append([v,v+1])
    
for r in pairs:
    p1 = r[0]
    p2 = r[1]
    print("%d, %d"%(p1,p2))
    
    ratio = [P[p2-1][q]/P[p1-1][q] for q in range(endindex)]
    meanratio = np.mean(ratio[-numavg:])

    pytools.findres(p1,p2,t,ml,Om,w,numavg,meanratio)

figP, axP = plt.subplots(2,figsize=(12,12))
#axP[0].set_ylim([1.,2.1])

for r in pairs:
    p1 = r[0]
    p2 = r[1]
    print("%d, %d"%(p1,p2))
    
    axP[0].plot(xs[1,:],ys[1,:], ',')
    axP[1].plot(t[0],a[1], ',')
    #ratio = [P[p2-1][q]/P[p1-1][q] for q in range(endindex)]
    #axP.plot(t[0][::], ratio, ',')

print(a[0][-10:])
#axP[1].set_ylim([0.,0.1])
#axP[1].plot(t[0][::], e[0][::], ',')
#axP[1].plot(t[0][::], e[1][::], ',')
plt.show()'''
figP, axP = plt.subplots(1,figsize=(12,12))
axP.plot(xs[1,:],ys[1,:],'.')


print(_xs)
print(math.sqrt(_xs[-1]**2 + _ys[-1]**2))

plt.show()
'''#assumes no particles lost
_xs = numpy.zeros((rebound.get_N(),N_output))
_ys = numpy.zeros((rebound.get_N(),N_output))
_zs = numpy.zeros((rebound.get_N(),N_output))

_ts = ts[:N_output]
_xs = xs[::,:N_output]
_ys = ys[::,:N_output]
_zs = zs[::,:N_output]

print(numpy.sqrt(_xs[1]**2 + _ys[1]**2)[-1])'''

#pytools.plot_freq_spectrum(_ts,_xs[1],2*numpy.pi/2.e6, 2*numpy.pi/100.,log=False)
#print(2*numpy.pi/plot_freq_spectrum(_ts,_xs[1],2*numpy.pi/2.e6, 2*numpy.pi/100.,log=False))
#getPeriods(_ts,_xs[1],_ys[1],_zs[1])

#fig, ax = plt.subplots(figsize=(6,6))
#ax.scatter(_xs[1],_ys[1])
#plt.show()
print("Took {0} seconds".format(time.time() - start_time))