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
            print("Planet {} had a*e > {} AU\n".format(i,dr_thresh))
            return True

        if a[i] == 0.:
            a[i] = o.a
            
        if math.fabs(o.a - a[i]) > dr_thresh:
            print("Planet {} had semimajor axis jump by more than {} AU".format(i, dr_thresh))
            return True
        
        a[i] = o.a
        com = pytools.get_center_of_mass(com, particles[i]) 
    return

def integrate(args):    
    mass, taue, ctr, dr_thresh,folder = args
    
    G = 4*math.pi**2
    rebound.set_G(G)
    starmass = 0.55
    
    tmax = 1.e7

    sun = rebound.Particle(m=starmass)
    rebound.particle_add(sun)

    a0s = [0.,8.,20.,48.,69.,100.]#[0.,13.6,33.3,71.2,93.0]

    for a in a0s[1:]:
        p = pytools.kepler_particle(m=mass,primary=sun,a=a,anom=0.,e=0.,omega=0.,inc=0.,Omega=0.)
        rebound.particle_add(p)

    particles = rebound.particles_get()    

    N = rebound.get_N()  
    
    rebound.add_e_damping([taue for i in range(N)])

    while rebound.get_t()<tmax:
        _t = rebound.get_t()
        rebound.step()
        breakflag = check_jumps(a0s,particles, dr_thresh)
        if breakflag is True:
            break
    
    with open(folder+'/eos/m_{0:.1e}_taue_{1:.1e}.txt'.format(mass,taue), mode='a') as f:
        f.write("{0}\t{1:.3e}\n".format(ctr,_t))

def main(argv,folder):
    try:
        opts, args = getopt.getopt(argv, "m:t:", ["mass=", "taue="])
    except getopt.GetoptError:
        print('Have to pass a mass and taue to this function')
        sys.exit(2)
        
    for opt, arg in opts:
        if opt in ("-m", "--mass"):
            mass = float(arg)   
        elif opt in ("-t", "--taue"):
            taue = float(arg)
    
    dr_thresh = 5.   
    args = []
    for it in range(24):
        args.append((mass,taue,it,dr_thresh,folder))
      
    pool = InterruptiblePool()
    pool.map(integrate, args)

if __name__ == "__main__":
    folder = "queenssystem"
    try:
        os.mkdir(folder)
    except OSError:
        pass
    try:
        os.mkdir(folder+"/eos")
    except OSError:
        pass
    main(sys.argv[1:],folder)   
      


