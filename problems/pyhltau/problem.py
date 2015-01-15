import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy
import math
import os
import pytools
import random
import matplotlib.pyplot as plt
import scipy.signal as signal

# function to calculate the non-keplerian forces
def additional_forces():
    com = particles[0] # calculate forces with respect to center of mass
    for i in range(1,N):
        if (tau_e[i]!=0 or tau_a[i]!=0):
            p = particles[i]
            dvx = p.vx-com.vx
            dvy = p.vy-com.vy
            dvz = p.vz-com.vz
            
            if tau_a[i]!=0.:     # Migration is on
                p.ax -=  dvx/(2.*tau_a[i])
                p.ay -=  dvy/(2.*tau_a[i])
                p.az -=  dvz/(2.*tau_a[i])
            
            # For ecc or inc damping, need h and e vectors
            if (tau_e[i]!=0 or tau_i[i]!=0): 
                mu = G*(com.m + p.m)
                dx = p.x-com.x
                dy = p.y-com.y
                dz = p.z-com.z
                
                hx = dy*dvz - dz*dvy
                hy = dz*dvx - dx*dvz
                hz = dx*dvy - dy*dvx
                h = numpy.sqrt ( hx*hx + hy*hy + hz*hz )
                v = numpy.sqrt ( dvx*dvx + dvy*dvy + dvz*dvz )
                r = numpy.sqrt ( dx*dx + dy*dy + dz*dz )
                vr = (dx*dvx + dy*dvy + dz*dvz)/r
                ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx )
                ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy )
                ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz )
                e = numpy.sqrt( ex*ex + ey*ey + ez*ez )    # eccentricity
                 
                if (tau_e[i]!=0):                           # ecc damping
                    a = -mu/( v*v - 2.*mu/r )               # semimajor axis
                    prefac1 = 1./(1.-e*e) /tau_e[i]/1.5
                    prefac2 = 1./(r*h) * numpy.sqrt(mu/a/(1.-e*e))/tau_e[i]/1.5
                    p.ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2
                    p.ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2
                    p.az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2
                
                if (tau_i[i]!=0):                           # inc damping
                    p.az += -2*dvz/tau_i[i]
                    prefac = (hx*hx + hy*hy)/h/h/tau_i[i]
                    p.ax += prefac*dvx
                    p.ay += prefac*dvy
                    p.az += prefac*dvz
        
        com = pytools.get_center_of_mass(com,particles[i])

G = 4.*math.pi**2
rebound.set_G(G)  
#rebound.set_dt(0.0001) 

starmass = 0.55     # in solar masses
N = 3               # including central star    
a = [0.,42.,53.]    # AU
M = 1.e-4           # solar masses
e = 0.01
i = 1.e-8
taue=1.e4
taua=-taue*100
outputdelta=10000.

tmax = 1.e6

# Add particles
sun = rebound.Particle(m=starmass,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
rebound.particle_add(sun)                  # Star at origin (zeros by default for pos & vel)

tau_a = numpy.zeros(N)             # include star and ignore the 0 index
tau_e = numpy.zeros(N)             # so we can use same numbers in 
tau_i = numpy.zeros(N)             # additional_forces

for j in range(1,N):
    rebound.particle_add(pytools.add_planet_3D(G,sun,M,a[j],e,i,0.,0.,0.))
    #rebound.particle_add(pytools.add_planet_3D(G,sun,M,a[j],e,i,random.uniform(0,2*math.pi),random.uniform(0,2*math.pi),random.uniform(0,2*math.pi)))
    tau_a[j] = taua if j==1 else 0.
    tau_e[j] = taue
    tau_i[j] = tau_e[j]

rebound.move_to_center_of_momentum()
particles = rebound.particles_get()
if (N != rebound.get_N()): raise ValueError("Number of objects added != # in simulation")

rebound.set_additional_forces(additional_forces)

# timestep counter
steps = 0 

last_t = -1e6

xs = numpy.zeros((rebound.get_N(),numpy.round(tmax/outputdelta)))
ys = numpy.zeros((rebound.get_N(),numpy.round(tmax/outputdelta)))
zs = numpy.zeros((rebound.get_N(),numpy.round(tmax/outputdelta)))
ts = numpy.zeros(numpy.round(tmax/outputdelta))
N_output = 0

while rebound.get_t()<tmax:
    t = rebound.get_t()
    if t - last_t > outputdelta:
        ts[N_output] = t
        for i in range(rebound.get_N()):
            xs[i,N_output] = particles[i].x
            ys[i,N_output] = particles[i].y
            zs[i,N_output] = particles[i].z
            #print(particles[i].x)
            
        N_output += 1
        last_t = t
    rebound.step()
    
#assumes no particles lost
_xs = numpy.zeros((rebound.get_N(),N_output))
_ys = numpy.zeros((rebound.get_N(),N_output))
_zs = numpy.zeros((rebound.get_N(),N_output))

_ts = ts[:N_output]
_xs = xs[::,:N_output]
_ys = ys[::,:N_output]
_zs = zs[::,:N_output]

print(numpy.sqrt(_xs[1]**2 + _ys[1]**2))

#pytools.plot_freq_spectrum(_ts,_xs[1],2*numpy.pi/2.e6, 2*numpy.pi/100.,log=False)
#print(2*numpy.pi/plot_freq_spectrum(_ts,_xs[1],2*numpy.pi/2.e6, 2*numpy.pi/100.,log=False))
#getPeriods(_ts,_xs[1],_ys[1],_zs[1])

#fig, ax = plt.subplots(figsize=(6,6))
#ax.scatter(_xs[1],_ys[1])
plt.show()