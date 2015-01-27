#!/usr/bin/python

# Import the rebound module
import sys; sys.path.append('../')
import rebound
# Import other modules
import math

def get_center_of_mass(p1, p2):
    com = rebound.Particle()
    com.x   = p1.x*p1.m + p2.x*p2.m        
    com.y   = p1.y*p1.m + p2.y*p2.m
    com.z   = p1.z*p1.m + p2.z*p2.m
    com.vx  = p1.vx*p1.m + p2.vx*p2.m
    com.vy  = p1.vy*p1.m + p2.vy*p2.m
    com.vz  = p1.vz*p1.m + p2.vz*p2.m
    com.m   = p1.m + p2.m
    if (com.m>0.):
        com.x  /= com.m
        com.y  /= com.m
        com.z  /= com.m
        com.vx /= com.m
        com.vy /= com.m
        com.vz /= com.m
    else:
        raise ValueError("Both particles had zero mass")
    
    return com

# mod2pi helper function.
TWOPI = 2.*math.pi
def mod2pi(f):
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

def get_E(e,M):
    E = M if e<0.8 else math.pi   
        
    F = E - e*math.sin(M) - M
    for i in range(100):
        E = E - F/(1.0-e*math.cos(E))
        F = E - e*math.sin(E) - M
        if math.fabs(F)<1e-16:
            break
    E = mod2pi(E)
    return E

def kepler_particle(m, 
                primary,    # central body (rebound.Particle object)
                a,          # semimajor axis
                anom=0.,  # anomaly
                e=0.,     # eccentricity
                omega=0., # argument of pericenter    
                inc=0.,   # inclination
                Omega=0., # longitude of ascending node
                MEAN=False):    # mean anomaly
    """Returns a particle structure initialized with the passed set of 
    orbital elements. Mass (m), primary and 'a' are required (see Parameters
    below, and any orbital mechanics text, e.g., Murray & Dermott
    Solar System Dynamics for definitions). Other values default to zero.
    All angles should be passed in radians. Units are set by the 
    gravitational constant G (default = 1.). If MEAN is set to True, anom is
    taken as the mean anomaly, rather than the true anomaly.  
    
    Usage
    _____
    primary = rebound.Particle(m=1.) # particle with unit mass at origin & v=0 
    
    # test particle (m=0) with specified elements using mean anomaly
    p = kepler_particle(m=0.,primary=primary,a=2.5, anom=math.pi/2,e=0.3,
                        omega=math.pi/6,inc=math.pi/3,Omega=0.,MEAN=True)
       
    # m=0.1 particle on circular orbit math.pi/4 from x axis in xy plane
    p = kepler_particle(0.1,primary,2.5,math.pi/4)   
                   
    Parameters
    __________
    m       : (float)            Mass of the particle
    primary : (rebound.Particle) Particle structure for the central body
    a       : (float)            Semimajor axis
    anom    : (float)            True anomaly (default).
                                 Mean anomaly if MEAN is set to True
    e       : (float)            Eccentricity
    omega   : (float)            Argument of pericenter
    inc     : (float)            Inclination (to xy plane)
    Omega   : (float)            Longitude of the ascending node
    MEAN    : (boolean)          If False (default), anom = true anomaly
                                 If True, anom = mean anomaly
    
    Returns
    _______
    A rebound.Particle structure initialized with the given orbital parameters
    """
    
    if not(0.<=e<1.): raise ValueError('e must be in range [0,1)') 
    # not sure if these equations work for parabolic/hyperbolic obits
    if not(0.<=inc<=math.pi): raise ValueError('inc must be in range [0,pi]')
    
    if MEAN is True: # need to calculate f
        E = get_E(e,anom)        
        f = mod2pi(2.*math.atan(math.sqrt((1.+ e)/(1. - e))*math.tan(0.5*E)))
    else:
        f = anom
    
    cO = math.cos(Omega)
    sO = math.sin(Omega)
    co = math.cos(omega)
    so = math.sin(omega)
    cf = math.cos(f)
    sf = math.sin(f)
    ci = math.cos(inc)
    si = math.sin(inc)
        
    r = a*(1.-e**2)/(1.+e*cf)
    
    # Murray & Dermott Eq. 2.122
    x  = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    y  = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    z  = primary.z + r*(so*cf+co*sf)*si
    
    n = math.sqrt(rebound.get_G()*(primary.m+m)/(a**3))
    v0 = n*a/math.sqrt(1.-e**2)
    
    # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)

    return rebound.Particle(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)

TINY=1.e-308
MIN_REL_ERROR = 1.e-12
# from tools.c.  Converts cartesian elements to orbital elements.

def p2orbit(p, primary,verbose=False):
    """ Returns a rebound.Orbit object with the keplerian orbital elements
    corresponding to p (rebound.Particle) around the central body primary
    (rebound.Particle). Edge cases will return values set to None. If
    verbose is set to True (default=False), error messages are printed
    when a breakout condition is met.
    
    Usage
    _____
    orbit = p2orbit(p,primary)
    print(orbit.e) # gives the eccentricity
    
    orbit = p2orbit(p,primary,verbose=True) # will print out error msgs
    
    Parameters
    __________
    p        : (rebound.Particle) particle for which orbital elements are sought
    primary  : (rebound.Particle) central body
    verbose  : (boolean)          If set to True, will print out error msgs
    
    Returns
    _______
    A rebound.Orbit object (with member variables for the orbital elements)
    """
    o = rebound.Orbit()
    if primary.m <= TINY: 
        if verbose is True:
            print("Star has no mass.")
        return o                            # all values set to None 
      
    dx = p.x - primary.x
    dy = p.y - primary.y
    dz = p.z - primary.z
    o.r = math.sqrt ( dx*dx + dy*dy + dz*dz )
    if o.r <= TINY: 
        if verbose is True:
            print('Particle and primary positions are the same.')
        return o
   
    dvx = p.vx - primary.vx
    dvy = p.vy - primary.vy
    dvz = p.vz - primary.vz
    v = math.sqrt ( dvx*dvx + dvy*dvy + dvz*dvz )
    
    mu = rebound.get_G()*(p.m+primary.m)
    o.a = -mu/( v*v - 2.*mu/o.r )               # semi major axis
    
    h0 = (dy*dvz - dz*dvy)                      # angular momentum vector
    h1 = (dz*dvx - dx*dvz)
    h2 = (dx*dvy - dy*dvx)
    o.h = math.sqrt ( h0*h0 + h1*h1 + h2*h2 )   # abs value of angular momentum 
    if o.h/o.r/v <= MIN_REL_ERROR: 
        if verbose is True:
            print('Particle orbit is radial.')
        return o
    
    vr = (dx*dvx + dy*dvy + dz*dvz)/o.r
    e0 = 1./mu*( (v*v-mu/o.r)*dx - o.r*vr*dvx )
    e1 = 1./mu*( (v*v-mu/o.r)*dy - o.r*vr*dvy )
    e2 = 1./mu*( (v*v-mu/o.r)*dz - o.r*vr*dvz )
    o.e = math.sqrt( e0*e0 + e1*e1 + e2*e2 )   # eccentricity
                  
    o.P = 2.*math.pi*math.sqrt( o.a*o.a*o.a/mu )  # period
    o.inc = math.acos( h2/o.h )               # inclination (wrt xy-plane)
                                                # if pi/2<i<pi it's retrograde
    n0 = -h1                                    # node vector 
    n1 =  h0                                    # in xy plane => no z component     
    n = math.sqrt( n0*n0 + n1*n1 )
    er = dx*e0 + dy*e1 + dz*e2               
    if (n/o.r/v<=MIN_REL_ERROR or o.inc<=MIN_REL_ERROR):# we are in the xy plane
        o.Omega=0.
        if (o.e <= MIN_REL_ERROR):              # omega not defined for circular orbit
            o.omega = 0.
        else:
            if (e1>=0.):
                o.omega=math.acos(e0/o.e)
            else:
                o.omega = 2.*math.pi-math.acos(e0/o.e)
    else:                                        
        if (o.e <= MIN_REL_ERROR):
            o.omega = 0.
        else:
            if (e2>=0.):                        # omega=0 if perictr at asc node
                o.omega=math.acos(( n0*e0 + n1*e1 )/(n*o.e))  
            else:
                o.omega=2.*math.pi-math.acos(( n0*e0 + n1*e1 )/(n*o.e)) 
        if (n1>=0.):
            o.Omega = math.acos(n0/n)
        else:
            o.Omega=2.*math.pi-math.acos(n0/n)  # Omega=longitude of asc node
                                                # taken in xy plane from x axis
    
    if (o.e<=MIN_REL_ERROR):                    # circular orbit
        o.f=0.                                  # f has no meaning
        o.l=0.
    else:
        cosf = er/(o.e*o.r)
        cosea = (1.-o.r/o.a)/o.e
        
        if -1.<=cosf<=1.:                       # failsafe
            o.f = math.acos(cosf)
        else:
            o.f = math.pi/2.*(1.-cosf)
        
        if -1.<=cosea<=1.:  
            ea  = math.acos(cosea)  
        else:
            ea = math.pi/2.*(1.-cosea)     
        
        if (vr<0.):
            o.f=2.*math.pi-o.f    
            ea =2.*math.pi-ea
    
        o.l = ea -o.e*math.sin(ea) + o.omega+ o.Omega  # mean longitude
    
    return o
