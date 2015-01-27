import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy as np
import math
import os
import pytools
import pickle
import time
start_time = time.time()
        
rebound.set_G(4.*math.pi**2)  

starmass = 0.55     # in solar masses
N = 6              # including central star    
a0 = 71.
afac = 1.24
a = [0.,13.6,33.3,a0,a0*afac,a0*afac**2]    # AU

e = 1.e-8
inc = 1.e-8
taue=1.e7
taua=taue*100
taues = (0.,taue,taue,taue,taue,taue)
tauis = (0.,taue,taue,taue,taue,taue)
tauas = (0.,0.,0.,0.,0.,taua)

mass = 8.773e-5 if taue == 1.e4 else 1.501e-5 # use 1.5e-5 Msun = 5 Mearth for anything
# but taue = 1e4, as these will capture.  For taue=1e4 (taua=1e6), the libration 
# timescale is too long for 5 Mearth, so have to raise to 29.2 Mearth in order to 
# shorten the libration timescale and capture.

tmax = taua/3.
outputdelta=tmax/10000. # so we get about 10k points

# Add particles
star = rebound.Particle(m=starmass,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
rebound.particle_add(star)                  # Star at origin (zeros by default for pos & vel)

for j in range(1,N):
    rebound.particle_add(pytools.kepler_particle(m=mass,primary=star,a=a[j],anom=0.,e=e,omega=0.,
                                             inc=inc,Omega=0.))

rebound.init_damping_forces()
rebound.add_migration(tauas)
rebound.add_e_damping(taues)
rebound.add_i_damping(tauis)

rebound.move_to_center_of_momentum()
particles = rebound.particles_get()

# timestep counter
steps = 0 
last_t = -1e6
N_output = 0

with open('taue_{0:.1e}.txt'.format(taue), mode='w') as f:
    pass # overwrite previous file

while rebound.get_t()<tmax:
    _t = rebound.get_t()
    if _t - last_t > outputdelta:
        with open('taue_{0:.1e}.txt'.format(taue), mode='a') as f:
            com = particles[0]
            for i in range(1,rebound.get_N()):
                o = pytools.p2orbit(particles[i],com)#particles[0])
                f.write('{0:.16e}\t{1:.16e}\t{2:.16e}\t{3:.16e}\t{4:.16e}\t{5:.16e}\t{6:.16e}\t{7:.16e}\t{8:.16e}\n'.format(_t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f))
                pytools.get_center_of_mass(com, particles[i])   
        N_output += 1
        last_t = _t
    
    rebound.step()

with open('taue_{0:.1e}_starter.bin'.format(taue), 'wb') as f:
    pickle.dump(([particles[q] for q in range(N)],tauas,taues,tauis),f)
        
print("Took {0} seconds".format(time.time() - start_time))
