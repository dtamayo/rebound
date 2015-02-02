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
    mthresh, taue, ctr, dr_thresh = args
    try:
        with open('/Users/dtamayo/Desktop/starters/m_{0:.2e}_taue_{1:.2e}_{2:1d}.pickle'.format(mthresh,taue,ctr), 'r') as f:
            data = pickle.load(f)
    except (IOError, OSError) as e:
        with open('eos/m_{0:.1e}_taue_{1:.1e}.txt'.format(mthresh,taue), mode='a') as f:
            f.write("{0}\t{1:.3e}\n".format(ctr,0.))
        return
        
    atrack = [0.,0.,0.,0.,0.,0.]
    rebound.reset()
    
    outputdelta=10.
    rebound.set_G(4.*math.pi**2)  
    tmax = 1.e6
    
    prevparticles = data[0]
    tauas = data[1]
    taues = data[2]
    tauis = data[3]
    
    for p in prevparticles: # data[0] is a list of the particles
        rebound.particle_add(p)

    N = rebound.get_N()  
    rebound.init_damping_forces()
    rebound.add_migration(tauas)
    rebound.add_e_damping(taues)
    rebound.add_i_damping(tauis)
    
    particles = rebound.particles_get()

    e = []
    t = []
    a = []
    P = []
    
    for q in range(N-1):
        e.append([])
        t.append([])   
        a.append([])
        P.append([])     
        
    last_t = -1e6
    N_output = 0

    while rebound.get_t()<tmax:
        _t = rebound.get_t()
        if _t - last_t > outputdelta:
            com = particles[0]
            for i in range(1,N):
                o = pytools.p2orbit(particles[i],com)
                t[i-1].append(_t)
                e[i-1].append(o.e)
                a[i-1].append(o.a)
                P[i-1].append(o.P)
                pytools.get_center_of_mass(com, particles[i])    
            N_output += 1
            last_t = _t
        rebound.step()
        breakflag = check_jumps(atrack,particles, dr_thresh)
        if breakflag is True:
            break
    
    with open('eos/m_{0:.1e}_taue_{1:.1e}.txt'.format(mthresh,taue), mode='a') as f:
        f.write("{0}\t{1:.3e}\n".format(ctr,_t))

def main(argv):
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
        args.append((mass,taue,it,dr_thresh))
      
    pool = InterruptiblePool()
    pool.map(integrate, args)

if __name__ == "__main__":
    try:
        os.mkdir("eos")
    except OSError:
        pass
    main(sys.argv[1:])   
      


