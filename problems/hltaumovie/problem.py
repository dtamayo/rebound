import sys; sys.path.append('../../python_modules/')
import rebound
# Import other modules
import numpy as np
import math
import os
import pytools
import random

def check_jumps(particles,r_thresh):
    for i in range(1,rebound.get_N()):
        r = np.sqrt(particles[i].x**2 + particles[i].y**2 + particles[i].z**2)
        if r > r_thresh:
            return True

    return False

G = 4*math.pi**2
rebound.set_G(G)
starmass = 0.55

tmax = 2.e4

sun = rebound.Particle(m=starmass)
rebound.particle_add(sun)

a0s = [0.,13.6,33.3,65.1,77.3,93.0]
mass = 5.e-4
for a in a0s[1:]:
    p = pytools.kepler_particle(m=mass,primary=sun,a=a,anom=random.uniform(0,2*math.pi),e=0.,omega=0.,inc=0.,Omega=0.)
    rebound.particle_add(p)

particles = rebound.particles_get()    

N = rebound.get_N()  

taue = 1.e5
rebound.add_e_damping([taue for i in range(N)])

tprev = -100.
r_thresh = 1000.
output_dt = 10.
rebound.move_to_center_of_momentum()

while rebound.get_t()<tmax:
    _t = rebound.get_t()
    
    rebound.step()
    breakflag = check_jumps(particles, r_thresh)
        if breakflag is True:
            break
    if _t - tprev > output_dt:
        with open('m_{0:.1e}_taue_{1:.1e}.txt'.format(mass,taue), mode='a') as f:
            for i in range(1,rebound.get_N()):
            o = particles[i].get_orbit(particles[0])
            f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(i,_t,o.a,o.e,o.inc,o.Omega.o.omega,o.l,o.P,o.f))

    tprev = _t
