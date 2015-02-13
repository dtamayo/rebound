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
from interruptible_pool import InterruptiblePool
import getopt
import matplotlib.pyplot as plt

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

def initialize(args):    
    mass,taue,k,a0,atarget,a_error,iterctr,last_mds,it,folder,dr_thresh = args
    
    print("Iteration {0}:".format(iterctr))
    # if we're below the following masses, then we're stable even out of resonance (and the below code
    # starts with a larger mass value to grow from), so automatically set to 1e6
    if taue == 1.e4:
        if mass < 2e-4:
            with open(folder+'/eos/m_{0:.1e}_taue_{1:.1e}.txt'.format(mass,taue), mode='a') as f:
                f.write("{0}\t{1:.3e}\t{2:.1f}\t{3:.1f}\t{4:.1f}\n".format(it,1.e6,0.,0.,0.,0))
            return
    else:
        if mass < 6.e-5:
            with open(folder+'/eos/m_{0:.1e}_taue_{1:.1e}.txt'.format(mass,taue), mode='a') as f:
                f.write("{0}\t{1:.3e}\t{2:.1f}\t{3:.1f}\t{4:.1f}\n".format(it,1.e6,0.,0.,0.,0))
            return
        
    rebound.reset()
    rebound.set_G(4.*math.pi**2)  

    starmass = 0.55     # in solar masses
    N = 5              # including central star    

    e = 1.e-8
    inc = 1.e-8
    taua=taue*k
    taues = [taue]*N
    tauas = [0.]*(N-1)
    tauas.append(taua)

    m0 = 2.e-4 if taue == 1.e4 else 6.e-5 # use 3.e-5 Msun = 10 Mearth for anything
# but taue = 1e4, as these will capture.  For taue=1e4 (taua=1e6), the libration 
# timescale is too long for 10 Mearth, so have to raise to ~39 Mearth in order to 
# shorten the libration timescale and capture.  This is stable in non-resonant case,
# so no problem starting this high, since masses below this would only be more stable

    tmax = taua/2.
    outputdelta=tmax/1000. # so we get about 10k points

    # Add particles
    star = rebound.Particle(m=starmass,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
    rebound.particle_add(star)                  # Star at origin (zeros by default for pos & vel)

    for j in range(1,N):
        rebound.particle_add(pytools.kepler_particle(m=m0,primary=star,a=a0[j],anom=random.uniform(0.,2*math.pi),e=e,omega=0.,
                                             inc=inc,Omega=0.))

    rebound.add_migration(tauas)
    rebound.add_e_damping(taues)


    rebound.move_to_center_of_momentum()
    particles = rebound.particles_get()

    # timestep counter
    steps = 0 
    last_t = -1e6
    N_output = 0

    phi343 = [] # res angle between planet 3 & 4, with pericenter of 3 etc.
    phi344 = []
    
    l = [0]*rebound.get_N() # list of N 0s
    po = [0]*rebound.get_N()
    a = [[] for i in range(rebound.get_N())]
    t = []

    while rebound.get_t()<tmax:
        _t = rebound.get_t()
        if _t - last_t > outputdelta:
            com = particles[0]
            for i in range(1,rebound.get_N()):
                o = pytools.p2orbit(particles[i],com)#particles[0])
                l[i] = o.l
                po[i] = o.Omega + o.omega
                a[i].append(o.a)
                pytools.get_center_of_mass(com, particles[i])   
            phi343.append(pytools.mod2pi(3*l[4] - 2*l[3] - po[3]))
            phi344.append(pytools.mod2pi(3*l[4] - 2*l[3] - po[4]))
            
            t.append(_t)
            N_output += 1
            last_t = _t
    
        rebound.step()
    
    o = pytools.p2orbit(particles[rebound.get_N()-1],particles[0])
    tlib = 0.078*(particles[N-1].m/0.55)**(-2./3.)*o.P
    tlib_ind = int(tlib // outputdelta)

    lib_thresh = 0.5

    if np.std(phi343[-tlib_ind:]) > lib_thresh or np.std(phi344[-tlib_ind:]) > lib_thresh:
        print("Iteration {0}: Didn't capture into resonance".format(iterctr))
        fig,axs = plt.subplots(2)
        axs[0].set_title(r'After Capture: $a_3,a_4$ = '+' {0:.1f}\t{1:.1f}'.format(np.mean(a[3][-tlib_ind:]),np.mean(a[4][-tlib_ind:])))
        axs[0].plot(t,phi343)
        axs[0].plot(t,phi344)
        axs[0].set_ylabel(r'$\phi s$')

        for i in range(1,N):
            axs[1].plot(t,a[i])
        axs[1].set_ylabel('a')
        axs[1].set_xlabel('t (yrs)')

        plt.savefig(folder+"/capinresprobs/m_{0:.1e}_taue_{1:.1e}_k{2:.1e}.png".format(mass,taue,k))
        plt.show()
        with open(folder+'/capinresprobs.txt', 'a') as f:
            f.write("{0:.1e}\t{1:1e}\t{2:3f}\t{3:.2f}\t{4:.2f}\t{5:.2f}\t{6:.2f}\t{7:.2f}\t{8:.2f}\n"
            .format(mass,taue,k,np.mean(a[3][-tlib_ind:]),np.mean(a[4][-tlib_ind:]),np.mean(phi343[-tlib_ind:]), np.std(phi343[-tlib_ind:]),
            np.mean(phi344[-tlib_ind:]), np.std(phi344[-tlib_ind:])))
        iterctr += 1
        initialize((mass,taue,k,a0,atarget,a_error,iterctr,last_mds,it,folder,dr_thresh))
      
    tmax += 1e7
    last_t = -1
    tprev = -1
    outputdelta=20.

    phi343 = [] # res angle between planet 3 & 4, with pericenter of 3 etc.
    phi344 = []
    
    l = [0]*rebound.get_N() # list of N 0s
    po = [0]*rebound.get_N()
    a = [[] for i in range(rebound.get_N())]
    t = []
    
    breakFlag = False
    #integrate until we reach the mass we want to start with (mass)
    while rebound.get_t()<tmax:
        _t = rebound.get_t()
        if _t - last_t > outputdelta:
            o = pytools.p2orbit(particles[N-1],particles[0])
            tlib = 0.078*(particles[N-1].m/0.55)**(-2./3.)*o.P
            deltaM = particles[N-1].m * outputdelta / 10. / tlib  
            for i in range(1,N):
                particles[i].m += deltaM
            if particles[3].m > mass:
                break
            
            com = particles[0]
            for i in range(1,rebound.get_N()):
                o = pytools.p2orbit(particles[i],com)#particles[0])
                l[i] = o.l
                po[i] = o.Omega + o.omega
                a[i].append(o.a)
                pytools.get_center_of_mass(com, particles[i])   
            phi343.append(pytools.mod2pi(3*l[4] - 2*l[3] - po[3]))
            phi344.append(pytools.mod2pi(3*l[4] - 2*l[3] - po[4]))
            t.append(_t)   
            last_t = _t
        if _t > 1.e3 and _t - tprev < 0.1: # timestep got too small (close encounter)
            breakFlag = True
            break
        tprev = _t
        rebound.step()
        
    if rebound.get_N() < 5:
        breakFlag = True
    if particles[3].m < mass:
        breakFlag = True
    
    if breakFlag is True:
        print("Went unstable while growing")
        fig,axs = plt.subplots(2)
        axs[0].set_title(r'After Growth: $a_3,a_4$ = '+' {0:.1f}\t{1:.1f}'.format(np.mean(a[3][-tlib_ind:]),np.mean(a[4][-tlib_ind:])))
        axs[0].plot(t,phi343)
        axs[0].plot(t,phi344)
        axs[0].set_ylabel(r'$\phi s$')

        for i in range(1,N):
            axs[1].plot(t,a[i])
        axs[1].set_ylabel('a')
        axs[1].set_xlabel('t (yrs)')

        plt.savefig(folder+"/growprobs/m_{0:.1e}_taue_{1:.1e}_k{2:.1e}.png".format(mass,taue,k))
        plt.show()
        with open(folder+'/growprobs.txt', 'a') as f:
            f.write("{0:.1e}\t{1:1e}\t{2:3f}\t{3:.2f}\t{4:.2f}\t{5:.2f}\t{6:.2f}\t{7:.2f}\t{8:.2f}\n"
            .format(mass,taue,k,np.mean(a[3][-tlib_ind:]),np.mean(a[4][-tlib_ind:]),np.mean(phi343[-tlib_ind:]), np.std(phi343[-tlib_ind:]),
            np.mean(phi344[-tlib_ind:]), np.std(phi344[-tlib_ind:])))

        # If we can't get it into resonance, set eos t=0
        with open(folder+'/eos/m_{0:.1e}_taue_{1:.1e}.txt'.format(mass,taue), mode='a') as f:
            f.write("{0}\t{1:.3e}\t{2:.1f}\t{3:.1f}\t{4:.1f}\t{5}\n".format(it,0.,-1,-1,-1,0))
            return
      
    af = [0]
    for i in range(1,rebound.get_N()):
        af.append(np.mean(a[i][-tlib_ind:]))
    diffs = [atarget[i] - af[i] for i in range(5)]
    iterfac = []
    
    print("Diffs = {0}".format(diffs[3:5]))
    
    iterflag = False
    
    for i in range(3,5):
        if math.fabs(diffs[i]) > a_error:
            iterfac.append(atarget[i]/af[i])
            iterflag = True
    
    if iterflag is True:
        meanfac = np.mean(iterfac)
        mds = np.mean([i**2 for i in diffs[3:5]]) # mean diffs squared
        for i in range(3,5):
            a0[i] *= meanfac
        print("iterfac = {0}, mds = {1}".format(meanfac,mds))
        if mds/last_mds > 1.: #if the squared differences for outer planets don't decrease by at least 20%,quit
            print("Trying to jiggle a0")
            last_mds = 1.6
            for i in range(3,5):
                a0[i] += 1.
            iterctr += 1
            return initialize((mass,taue,k,a0,atarget,a_error,iterctr,last_mds,it,folder,dr_thresh)) # if you don't return, returns within function called recursively won't work!

        iterctr += 1
        last_mds = mds
        if iterctr > 4:
            print("Ran past 4 iterations.  Exiting...")
            with open(folder+'/iterationprobs.txt', 'a') as f:
                f.write("{0:.1e}\t{1:1e}\t{2:1e}\t{3}\n".format(mass,taue,k,diffs))
            # If we can't get it into resonance, set eos t=0
            with open(folder+'/eos/m_{0:.1e}_taue_{1:.1e}.txt'.format(mass,taue), mode='a') as f:
                f.write("{0}\t{1:.3e}\t{2:.1f}\t{3:.1f}\t{4:.1f}\t{5}\n".format(it,0.,-3,-3,-3,0))
            return
        
        return initialize((mass,taue,k,a0,atarget,a_error,iterctr,last_mds,it,folder,dr_thresh)) # if you don't return, returns within function called recursively won't work!
        #google python return not working
    
    a0s = [0.]*rebound.get_N()
    rebound.set_t(0.)
    tmax = 1.e6

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
    with open(folder+'/eos/m_{0:.1e}_taue_{1:.1e}.txt'.format(mass,taue), mode='a') as f:
        f.write("{0}\t{1:.3e}\t{2:.1f}\t{3:.1f}\t{4:.1f}\t{5}\n".format(it,_t,diffs[-3],diffs[-2],diffs[-1],1))
    
    return
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
    
    a0 = 95.
    afac = 1.45
    a0 = [0.,13.6,33.3,a0,a0*afac]    # AU
    atarget = [0.,13.6,33.3,71.2,93.0]
    a_error = 3
    iterctr = 0
    last_mds = 1.e6
    k=100
    dr_thresh = 5.  
     
    args = []
    for it in range(24):
        args.append((mass,taue,k,a0,atarget,a_error,iterctr,last_mds,it,folder,dr_thresh))
    
    initialize(args[0])
    #pool = InterruptiblePool()
    #pool.map(initialize, args)

if __name__ == "__main__":
    folder = "4pres"
    try:
        os.mkdir(folder)
    except OSError:
        pass
    try:
        os.mkdir(folder+"/eos")
    except OSError:
        pass
    try:
        os.mkdir(folder+"/growprobs")
    except OSError:
        pass
    try:
        os.mkdir(folder+"/capinresprobs")
    except OSError:
        pass

    main(sys.argv[1:],folder)   