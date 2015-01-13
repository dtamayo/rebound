import sys; sys.path.append('../')
print(sys.version)
import rebound
# Import other modules
import numpy
import math
import os
import pytools
import random
import matplotlib.pyplot as plt
import scipy.signal as signal

def getPeriods(t,x,y):
    r = numpy.sqrt(x**2+y**2)
    f,Pper_spec = signal.welch(r,0.01,'flattop', scaling='spectrum')
    plt.semilogy(f,Pper_spec)
    plt.grid()
    plt.show()
    #minindices = argrelextrema(r,numpy.less)
    #xsmooth = smooth(x[minindices],window_len=10,window='hanning')
    #perieq0indices = argrelextrema(xsmooth,numpy.greater)

    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(t,r)

    plt.show()

    #precessionperiod = t[minindices][perieq0indices][1] - t[minindices][perieq0indices][0]
    #orbP = t[minindices[0][1]]-t[minindices[0][0]]
    #precPinorb = precessionperiod/orbP
# Set variables (defaults are G=1, t=0, dt=0.01)

G = 4.*math.pi**2
rebound.set_G(G)  

starmass = 0.55;    # in solar masses
Nplanets = 2    
a = [42.,53.]       # AU
M = 5.5e-5          # solar masses
e = 0.01
i = 0.01

tmax = 2.e5

# Add particles
sun = rebound.Particle(m=starmass,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
rebound.particle_add(sun)                  # Star at origin (zeros by default for pos & vel)

tau_a = numpy.zeros(Nplanets)
tau_e = numpy.zeros(Nplanets)
tau_i = numpy.zeros(Nplanets)

print("***Before adding particle to rebound***")
print(sun.x)
print(sun.y)
for j in range(Nplanets):
    p = pytools.add_planet_3D(G,sun,M,a[j],e,i,0.,0.,0.,MEAN=False)
    print(p.x)
    print(p.y)
    rebound.particle_add(p)
    #rebound.particle_add(pytools.add_planet_3D(G,sun,M,a[j],e,i,0.,0.,0.,MEAN=False))
    #rebound.particle_add(pytools.add_planet_3D(G,sun,M,a[i],e,i,random.uniform(0,2*math.pi),random.uniform(0,2*math.pi),random.uniform(0,2*math.pi)))
    tau_a[j] = -1.5e7 if j==0 else 0.
    tau_e[j] = 10000.
    tau_i[j] = tau_e[j]

particles = rebound.particles_get()
print("***After adding particle to rebound***")
for j in range(rebound.get_N()):
    print(particles[j].x)
    print(particles[j].y)
    
rebound.move_to_center_of_momentum()
# Get the particle data
# Note: this is a pointer and will automatically update as the simulation progresses
particles = rebound.particles_get()
# timestep counter
steps = 0 

last_t = 0
outputdelta=100.
xs = numpy.zeros((rebound.get_N(),numpy.round(tmax/outputdelta)))
ys = numpy.zeros((rebound.get_N(),numpy.round(tmax/outputdelta)))
ts = numpy.zeros(numpy.round(tmax/outputdelta))
N_output = 0

while rebound.get_t()<tmax:
    rebound.step()
    t = rebound.get_t()
    if t - last_t > outputdelta:
        ts[N_output] = t
        for i in range(rebound.get_N()):
            xs[i,N_output] = particles[i].x
            ys[i,N_output] = particles[i].y
            
        N_output += 1
        last_t = t

#assumes no particles lost
_xs = numpy.zeros((rebound.get_N(),N_output))
_ys = numpy.zeros((rebound.get_N(),N_output))

_ts = ts[:N_output]
_xs = xs[::,:N_output]
_ys = ys[::,:N_output]

getPeriods(_ts,_xs[1],_ys[1])

fig, ax = plt.subplots(figsize=(6,6))
ax.scatter(_xs[1],_ys[1])
plt.show()