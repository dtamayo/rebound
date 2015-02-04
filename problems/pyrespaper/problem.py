import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy as np
import math
import os
import pytools
import pickle
import time
import random
import matplotlib.pyplot as plt
start_time = time.time()
        
rebound.set_G(4.*math.pi**2)  

starmass = 0.55     # in solar masses
N = 6              # including central star    
a0 = 71.
afac = 1.24
a = [0.,13.6,33.3,a0,a0*afac,a0*afac**2]    # AU

e = 1.e-8
inc = 1.e-8
taue=1.e4
k = 100
taua=taue*k
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
    rebound.particle_add(pytools.kepler_particle(m=mass,primary=star,a=a[j],anom=random.uniform(0.,2*math.pi),e=e,omega=0.,
                                             inc=inc,Omega=0.))

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

phi343 = [] # res angle between planet 3 & 4, with pericenter of 3 etc.
phi344 = []
phi454 = []
phi455 = [] 
l = [0]*rebound.get_N() # list of N 0s
po = [0]*rebound.get_N()
a = [0]*rebound.get_N()
t = []
while rebound.get_t()<tmax:
    _t = rebound.get_t()
    if _t - last_t > outputdelta:
        with open('taue_{0:.1e}.txt'.format(taue), mode='a') as f:
            com = particles[0]
            for i in range(1,rebound.get_N()):
                o = pytools.p2orbit(particles[i],com)#particles[0])
                f.write('{0:.16e}\t{1:.16e}\t{2:.16e}\t{3:.16e}\t{4:.16e}\t{5:.16e}\t{6:.16e}\t{7:.16e}\t{8:.16e}\n'.format(_t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f))
                l[i] = o.l
                po[i] = o.Omega + o.omega
                a[i] = o.a
                pytools.get_center_of_mass(com, particles[i])   
            phi343.append(pytools.mod2pi(4*l[4] - 3*l[3] - po[3]))
            phi344.append(pytools.mod2pi(4*l[4] - 3*l[3] - po[4]))
            phi454.append(pytools.mod2pi(4*l[5] - 3*l[4] - po[4]))
            phi455.append(pytools.mod2pi(4*l[5] - 3*l[4] - po[5]))
            t.append(_t)
        N_output += 1
        last_t = _t
    
    rebound.step()

fig,axs = plt.subplots(4)
axs[0].set_title(r'$a_3,a_4,a_5$ = '+' {0:.1f}\t{1:.1f}\t{2:.1f}'.format(a[3],a[4],a[5]))
axs[0].plot(t,phi343)
axs[1].plot(t,phi344)
axs[2].plot(t,phi454)
axs[3].plot(t,phi455)
axs[0].set_ylabel(r'$\phi_{34}^i$')
axs[1].set_ylabel(r'$\phi_{34}^o$')
axs[2].set_ylabel(r'$\phi_{45}^i$')
axs[3].set_ylabel(r'$\phi_{45}^o$')
axs[3].set_xlabel('t (yrs)')

plt.savefig("inres/taue_{0:.1e}_k{1:}.png".format(taue,k))

with open('inres/taue_{0:.1e}_k{1:}.txt'.format(taue,k), 'w') as f:
    f.write("{0:.2f}\t{1:.2f}\t{2:.2f}\n".format(a[3],a[4],a[5]))
    f.write("{0:.2f}\t{1:.2f}\n".format(np.mean(phi343[-100:]), np.std(phi343[-100:])))
    f.write("{0:.2f}\t{1:.2f}\n".format(np.mean(phi344[-100:]), np.std(phi344[-100:])))
    f.write("{0:.2f}\t{1:.2f}\n".format(np.mean(phi454[-100:]), np.std(phi454[-100:])))
    f.write("{0:.2f}\t{1:.2f}\n".format(np.mean(phi455[-100:]), np.std(phi455[-100:])))

rebound.reset()

print(particles[1].x)
print("Took {0} seconds".format(time.time() - start_time))
