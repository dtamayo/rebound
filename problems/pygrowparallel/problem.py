import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy as np
import math
import os
import pytools
import random
import pickle
import getopt
import time

start_time = time.time()

def make_restarts():
    #make the first one the original simulation
    ps = []
    for q in range(N):
        ps.append(particles[q])
        if q == 1:
            print(particles[q].x)
    with open('starters/m_{0:.2e}_taue_{1:.2e}_{2:2d}.bin'.format(mthresh,taues[1],0), 'wb') as f:
        pickle.dump((ps,tauas,taues,tauis),f)
    #now perturb n_restarts-1 times
    for i in range(n_restarts-1):
        com = particles[0] # just use star as com.  We're then jiggling the heliocentric f (diff between heliocentric and jacobi shouldn't matter)
        ps = []
        ps.append(particles[0])
        for q in range(1,N):
            o = pytools.p2orbit(particles[q],com)
            o.f += np.random.uniform(-delta,delta)
            p = pytools.kepler_particle(particles[q].m,com,o.a,o.f,o.e,o.omega,o.inc,o.Omega)
            ps.append(p)
            if q == 1:
                print(p.x)
        with open('starters/m_{0:.2e}_taue_{1:.2e}_{2:2d}.bin'.format(mthresh,taues[1],i+1), 'wb') as f:
            pickle.dump((ps,tauas,taues,tauis),f)

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "m:", ["mass="])
    except getopt.GetoptError:
        print('Have to pass a mass to this function (and only a mass)')
        sys.exit(2)
        
    for opt, arg in opts:
        if opt in ("-m", "--mass"):
            mass = float(arg)
    
    print(mass)

    its = 3
    
    for taue in range(4,8):
        for i in range(its):
            success = grow(mass,10**taue)
            if success is True:
                break
  
def grow(mthresh, taue):
    n_restarts = 3
    delta = 2e-2

    outputdelta=20.
    rebound.set_G(4.*math.pi**2)  

    tmax = 2.e7 # something big, will break out

    with open('taue_{0:.1e}_starter.bin'.format(taue), 'rb') as f:
        data = pickle.load(f)

        prevparticles = data[0]
        tauas = data[1]
        taues = data[2]
        tauis = data[3]

    for p in prevparticles: # data[0] is a list of the particles
        rebound.particle_add(p)

        rebound.init_damping_forces()
        rebound.add_migration(tauas)
        rebound.add_e_damping(taues)
        rebound.add_i_damping(tauis)

    particles = rebound.particles_get()

    N = rebound.get_N()

    last_t = -1e6
    tprev = -1 # for keeping track of timestep getting too small separate from output
    
    while rebound.get_t()<tmax:
        _t = rebound.get_t()
        if _t - last_t > outputdelta:
            o = pytools.p2orbit(particles[i],particles[0])
            tlib = 0.078*(particles[i].m/0.55)**(-2./3.)*o.P
            deltaM = particles[i].m * outputdelta / 10. / tlib  
        last_t = _t
        
        for i in range(1,N):
            particles[i].m += deltaM
        if particles[3].m > mthresh:
            break
        if _t - tprev < 0.1: # timestep got too small (close encounter)
            return False
        tprev = _t
        rebound.step()

    if rebound.get_N() < 6:
        return False
    if particles[3].m < mthresh:
        return False
        
    make_restarts()
    with open('final_a_taue_{0:.1e}.txt'.format(taue), mode='a', encoding='utf-8') as f:
        f.write('{0:.2e}'.format(mthresh))
        for i in range(1,rebound.get_N()):
            o = pytools.p2orbit(particles[i].m,particles[0])
            f.write('\t{0:.2f}'.format(o.a))
        f.write('\n')
    return True

print("Took {0} seconds".format(time.time() - start_time))

if __name__ == "__main__":
    try:
        os.mkdir("starters")
    except OSError:
        pass
    main(sys.argv[1:])