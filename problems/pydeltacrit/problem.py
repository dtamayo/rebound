import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy as np
import math
import os
import pytools
import random
import pickle
from interruptible_pool import InterruptiblePool
import time
import getopt

start_time = time.time()

def check_jumps(a, particles, dr_thresh):
    com = particles[0]
    for i in range(1,rebound.get_N()):
        o = pytools.p2orbit(particles[i],com)

        if o.a*o.e > dr_thresh:
            print("Planet {0} had a*e > {1} AU\n".format(i,dr_thresh))
            return True

        if a[i] == 0.:
            a[i] = o.a
            
        if math.fabs(o.a - a[i]) > dr_thresh:
            print("Planet {0} had semimajor axis jump by more than {1} AU".format(i, dr_thresh))
            return True
        
        a[i] = o.a
        com = pytools.get_center_of_mass(com, particles[i]) 
    return

def integrate(args):    
    Delta,nRH,j,it,dr_thresh,folder = args
     
    rebound.reset()
    G = 4*math.pi**2
    rebound.set_G(G)
    
    starmass = 0.55
    a1 = 71.2
    tmax = 1.e6

    daOvera = (j/(j+1)*(1-Delta))**(-2./3.)-1.
    M = 1.5*(daOvera/nRH)**3*starmass
    a2 = a1*(1.+ daOvera)
    
    sun = rebound.Particle(m=starmass)
    rebound.particle_add(sun)

    a0s = [0.,a1,a2]

    for a in a0s[1:]:
        p = pytools.kepler_particle(m=M,primary=sun,a=a,anom=random.uniform(0,2*math.pi),e=0.,omega=0.,inc=0.,Omega=0.)
        rebound.particle_add(p)

    particles = rebound.particles_get()    

    N = rebound.get_N()  
    
    #rebound.add_e_damping([taue for i in range(N)])

    tprev = -2.
    dtthresh = 1.
    while rebound.get_t()<tmax:
        _t = rebound.get_t()
        rebound.step()
        breakflag = check_jumps(a0s,particles, dr_thresh)
        if breakflag is True:
            break
        if _t - tprev < dtthresh:
            break
    print(_t)
    with open(folder + '/N_{0:3f}.txt'.format(nRH), mode='a') as f:
        f.write("{0:.4e}\t{1:.4e}\n".format(Delta,_t))

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "D:n:j:", ["Delta=", "nRH=", "j="])
    except getopt.GetoptError:
        print('Have to pass a Delta, nRH and j')
        sys.exit(2)
        
    for opt, arg in opts:
        if opt in ("-D", "--Delta"):
            Delta = float(arg)   
        elif opt in ("-n", "--nRH"):
            nRH = float(arg)
        elif opt in ("-j", "--j"):
            j = float(arg)
    
    folder = "{0}to{1}".format(int(j+1), int(j))
    try:
        os.mkdir(folder)
    except OSError:
        pass
    
    dr_thresh = 5.   
    args = []
    for it in range(24):
        args.append((Delta,nRH,j,it,dr_thresh,folder))
     
    pool = InterruptiblePool()
    pool.map(integrate, args)

if __name__ == "__main__":
    main(sys.argv[1:])   
      


