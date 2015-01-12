import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy
import math
import os
import pytools
import random
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=numpy.ones(window_len,'d')
        else:  
                w=eval('numpy.'+window+'(window_len)')
        y=numpy.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

def getPeriods(t,x,y):
    r = numpy.sqrt(x**2+y**2)
    minindices = argrelextrema(r,numpy.less)
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
i = 1.e-8

tmax = 3.e3

# Add particles
sun = rebound.Particle(m=starmass,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
rebound.particle_add(sun)                  # Star at origin (zeros by default for pos & vel)

tau_a = numpy.zeros(Nplanets)
tau_e = numpy.zeros(Nplanets)
tau_i = numpy.zeros(Nplanets)

for i in range(Nplanets):
    rebound.particle_add(pytools.add_planet_3D(G,sun,M,a[i],0.6,0.,0.,0.,0.,MEAN=False))
    #rebound.particle_add(pytools.add_planet_3D(G,sun,M,a[i],e,i,random.uniform(0,2*math.pi),random.uniform(0,2*math.pi),random.uniform(0,2*math.pi)))
    tau_a[i] = -1.5e7 if i==0 else 0.
    tau_e[i] = 10000
    tau_i[i] = tau_e[i]
    
rebound.move_to_center_of_momentum()
# Get the particle data
# Note: this is a pointer and will automatically update as the simulation progresses
particles = rebound.particles_get()
# timestep counter
steps = 0 

last_t = 0
outputdelta=1.
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

print(N_output)
#assumes no particles lost
_xs = numpy.zeros((rebound.get_N(),N_output))
_ys = numpy.zeros((rebound.get_N(),N_output))

_ts = ts[:N_output]
_xs = xs[::,:N_output]
_ys = ys[::,:N_output]

getPeriods(_ts[:100],_xs[1][:100],_ys[1][:100])

fig, ax = plt.subplots(figsize=(6,6))
ax.scatter(_xs[1],_ys[1])
plt.show()