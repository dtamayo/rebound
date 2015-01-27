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

def jiggle(particles,delta):
    for q in range(1,rebound.get_N()):
        com = particles[0]
        for w in range(1,rebound.get_N()):
            com = pytools.get_center_of_mass(com, particles[w])
        print(com.x, com.y, com.z, com.vx, com.vy, com.vz, com.m)
        
        o = pytools.p2orbit(com, particles[q])
        o.f += np.random.uniform(-delta,delta)
        r = o.a*(1-o.e**2)/(1 + o.e*math.cos(o.f))
    
        particles[q].x = com.x + r*(math.cos(o.Omega)*math.cos(o.omega+o.f) - math.sin(o.Omega)*math.sin(o.omega+o.f)*math.cos(o.inc))
        particles[q].y = com.y + r*(math.sin(o.Omega)*math.cos(o.omega+o.f) + math.cos(o.Omega)*math.sin(o.omega+o.f)*math.cos(o.inc))
        particles[q].z = com.z + r*math.sin(o.omega+o.f)*math.sin(o.inc)
    
        n = math.sqrt(rebound.get_G()*(com.m+particles[q].m)/(o.a**3))
    
        # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
        particles[q].vx = com.vx + (n*o.a/math.sqrt(1-o.e*o.e))*((o.e+math.cos(o.f))*(-math.cos(o.inc)*math.cos(o.omega)*math.sin(o.Omega) - math.cos(o.Omega)*math.sin(o.omega)) - math.sin(o.f)*(math.cos(o.omega)*math.cos(o.Omega) - math.cos(o.inc)*math.sin(o.omega)*math.sin(o.Omega)))
        particles[q].vy = com.vy + (n*o.a/math.sqrt(1-o.e*o.e))*((o.e+math.cos(o.f))*(math.cos(o.inc)*math.cos(o.omega)*math.cos(o.Omega) - math.sin(o.omega)*math.sin(o.Omega)) - math.sin(o.f)*(math.cos(o.omega)*math.sin(o.Omega) + math.cos(o.inc)*math.cos(o.Omega)*math.sin(o.omega)))
        particles[q].vz = com.vz + (n*o.a/math.sqrt(1-o.e*o.e))*((o.e+math.cos(o.f))*math.cos(o.omega)*math.sin(o.inc) - math.sin(o.f)*math.sin(o.inc)*math.sin(o.omega))
        
        '''particles[q].x += np.random.uniform(-delta, delta)
        particles[q].y += np.random.uniform(-delta, delta)
        particles[q].z += np.random.uniform(-delta, delta)
        particles[q].vx += np.random.uniform(-delta, delta)
        particles[q].vy += np.random.uniform(-delta, delta)
        particles[q].vz += np.random.uniform(-delta, delta)    '''    
    '''ctr = 1
    while True:
        p = rebound.particle_get(ctr)
        if p is None:
            break
        print(ctr)
        ctr += 1
        p.x += np.random.uniform(-delta, delta)
        p.y += np.random.uniform(-delta, delta)
        p.z += np.random.uniform(-delta, delta)
        p.vx += np.random.uniform(-delta, delta)
        p.vy += np.random.uniform(-delta, delta)
        p.vz += np.random.uniform(-delta, delta)'''
        
rebound.set_G(4.*math.pi**2)  

starmass = 0.55     # in solar masses
N = 4              # including central star   
a0 = 42. 
afac = 1.24
a = [0.,a0,a0*afac,a0*afac**2]    # AU
mass = 1.501e-5       # solar masses
e = 1.e-8
inc = 1.e-8
taue=1.e4
taua=taue*100
taues = (0.,taue,taue,taue)
tauis = (0.,taue,taue,taue)
tauas = (0.,0.,0.,taua)

outputdelta=10.

tmax = 3.e6

# Add particles
star = rebound.Particle(m=starmass,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
rebound.particle_add(star)                  # Star at origin (zeros by default for pos & vel)

for j in range(1,N):
    rebound.particle_add(pytools.kepler_particle(m=mass,primary=star,a=a[j],anom=0.,e=e,omega=0.,
                                             inc=inc,Omega=0.))
'''    rebound.particle_add(pytools.init_planet(star,a=a[j],e=e,omega=random.uniform(0,2*math.pi),
                                             f=random.uniform(0,2*math.pi),inc=inc,
                                             Omega=random.uniform(0,2*math.pi), mass=mass))
'''
rebound.init_damping_forces()
rebound.add_migration(tauas)
rebound.add_e_damping(taues)
rebound.add_i_damping(tauis)

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

Om = []
a = []
e = []
w = []
ml = []
t = []
P = []

for q in range(N-1):
    Om.append([])
    a.append([])
    e.append([])
    w.append([])
    ml.append([])
    t.append([])
    P.append([])
    
o = rebound.Orbit()

_xs = []
_ys = []
_zs = []

flag = False

while rebound.get_t()<tmax:
    _t = rebound.get_t()
    if _t - last_t > outputdelta:
        com = particles[0]
        for i in range(1,rebound.get_N()):
            o = pytools.p2orbit(particles[i],com)#particles[0])
            t[i-1].append(_t)
            a[i-1].append(o.a)
            e[i-1].append(o.e)
            Om[i-1].append(o.Omega)
            w[i-1].append(o.omega)
            ml[i-1].append(o.l)
            P[i-1].append(o.P)
            xs[i,N_output] = particles[i].x
            ys[i,N_output] = particles[i].y
            zs[i,N_output] = particles[i].z
            pytools.get_center_of_mass(com, particles[i])    
        N_output += 1
        last_t = _t
    _xs.append(particles[1].x)
    _ys.append(particles[1].y)
    _zs.append(particles[1].z)
    rebound.step()
    if _t > 6.2e6:
        with open('starter.pickle', 'wb') as f:
            pickle.dump(([particles[q] for q in range(N)],tauas,taues,tauis),f)
        print(particles[1].x)
        break
    '''if _t > 1.3e10 and flag is False:
        print("Jiggling...")
        print(particles[1].x)
        jiggle(particles,np.pi/1000.)
        print(particles[1].x)
        flag = True'''
            

''''N = 3 #number of planets
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
axP[0].set_ylim([1.,2.1])

for r in pairs:
    p1 = r[0]
    p2 = r[1]
    print("%d, %d"%(p1,p2))
    
    ratio = [P[p2-1][q]/P[p1-1][q] for q in range(endindex)]
    axP[0].plot(t[0], ratio, ',')

#axP[1].set_ylim([0.,0.1])
#axP[1].plot(t[0][::], e[0][::], ',')
#axP[1].plot(t[0][::], e[1][::], ',')
plt.show()
'''
phi = pytools.resarg(4,-3,-1,0, ml, Om, w, 1, 2)
q = [e[1][j]*math.cos(phi[j]*math.pi/180.) for j in range(len(e[0]))]
p = [e[1][j]*math.sin(phi[j]*math.pi/180.) for j in range(len(e[0]))]
Pratio1 = [P[1][j]/P[0][j] for j in range(len(P[0]))]
Pratio2 = [P[2][j]/P[1][j] for j in range(len(P[0]))]

start = 0
end = len(t[0])

'''for ctr in range(end-start):
    mov,axmov = plt.subplots(2,figsize=(6,12))
    axmov[0].plot(t[0][start:start+ctr], phi[start:start+ctr], '.')
    axmov[1].plot(q[start:start+ctr],p[start:start+ctr],'.')
    axmov[1].set_xlim([-0.03,0.03])
    axmov[1].set_ylim([-0.03,0.03])
    plt.savefig('/Users/dtamayo/Desktop/mov/{:04d}.png'.format(ctr))

'''

fig,ax = plt.subplots(4)
ax[0].plot(t[0][start:end],Pratio1[start:end], ',')
ax[0].plot(t[0][start:end],Pratio2[start:end], ',')
ax[1].plot(t[0][start:end],phi[start:end], ',')
ax[2].plot(t[0][start:end],a[0][start:end],',')
ax[2].plot(t[0][start:end],a[1][start:end],',')
ax[2].plot(t[0][start:end],a[2][start:end],',')
ax[3].plot(t[0][start:end], e[0][start:end], ',')
ax[3].plot(t[0][start:end], e[1][start:end], ',')

print(len(P[0]))

start = 0
end = len(t[0])

fig2,ax2 = plt.subplots()
ax2.plot(q[start:end],p[start:end],',')
plt.axis('equal')
emax = max(e[1])
ax2.set_ylim([-emax,emax])
ax2.set_xlim([-emax,emax])

'''fig3,ax3 = plt.subplots(2)
ax3[0].plot(t[0][start:end], e[1][start:end], ',')
ax3[1].plot(t[0][start:end], phi[start:end], ',')
'''
print(Pratio1[-1])

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
plt.show()