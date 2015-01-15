#!/usr/bin/python

# Import the rebound module
import sys; sys.path.append('../')
import rebound
# Import other modules
import math
import numpy
import scipy.signal as signal
import matplotlib.pyplot as plt

def plot_freq_spectrum(t,vals,wmin,wmax,nfreq=1000,log=True):
    if log == True:
        lnminw = numpy.log(wmin)
        lnmaxw = numpy.log(wmax)
        logw = numpy.linspace(lnminw,lnmaxw,nfreq)
        w = numpy.asarray([numpy.e**i for i in logw])
    else:
        w = numpy.linspace(wmin,wmax,nfreq)
    
    pgram = signal.lombscargle(t,vals,w)
    
    fig,ax = plt.subplots()
    if log == True:
        ax.semilogx(w,numpy.sqrt(4*pgram/vals.shape[0]))
    else:
        ax.plot(w,numpy.sqrt(4*pgram/vals.shape[0]))
    plt.show()

def find_freq_peak(t,vals,wmin,wmax,nfreq=1000,showplot=False):   
    w = numpy.linspace(wmin,wmax,nfreq)
    pgram = signal.lombscargle(t,vals,w)
    
    if showplot==True:
        fig,ax = plt.subplots()
        ax.plot(w,numpy.sqrt(4*pgram/vals.shape[0]))
        plt.show()
        
    maxindex = numpy.argmax(pgram)
    return w[maxindex] 

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
    E = M if e<0.8 else numpy.pi   
        
    F = E - e*math.sin(M) - M
    for i in range(100):
        E = E - F/(1.0-e*math.cos(E))
        F = E - e*math.sin(E) - M
        if math.fabs(F)<1e-16:
            break
    E = mod2pi(E)
    return E

# Converts orbital elements to cartesian and returns a particle
# if Mean is set to True, the anomaly is taken as the mean anomaly
# if Mean is set to False, the anomaly is taken as the true anomaly
def add_planet_3D(G,com,mass,a,e,inc,Omega,omega,anom,MEAN=True):
    if MEAN==False:
        f = anom
    else:
        E = get_E(e,anom)        
        f = mod2pi(2.*math.atan(math.sqrt((1. + e)/(1. - e))*math.tan(0.5*E)))
    
    r = a*(1-e**2)/(1 + e*math.cos(f))
    
    # Murray & Dermott Eq. 2.122
    _x  = com.x + r*(math.cos(Omega)*math.cos(omega+f) - math.sin(Omega)*math.sin(omega+f)*math.cos(inc))
    _y  = com.y + r*(math.sin(Omega)*math.cos(omega+f) + math.cos(Omega)*math.sin(omega+f)*math.cos(inc))
    _z  = com.z + r*math.sin(omega+f)*math.sin(inc)
    
    n = math.sqrt(G*(com.m+mass)/(a**3))
    
    # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    _vx = com.vx + (n*a/math.sqrt(1-e*e))*((e+math.cos(f))*(-math.cos(inc)*math.cos(omega)*math.sin(Omega) - math.cos(Omega)*math.sin(omega)) - math.sin(f)*(math.cos(omega)*math.cos(Omega) - math.cos(inc)*math.sin(omega)*math.sin(Omega)))
    _vy = com.vy + (n*a/math.sqrt(1-e*e))*((e+math.cos(f))*(math.cos(inc)*math.cos(omega)*math.cos(Omega) - math.sin(omega)*math.sin(Omega)) - math.sin(f)*(math.cos(omega)*math.sin(Omega) + math.cos(inc)*math.cos(Omega)*math.sin(omega)))
    _vz = com.vz + (n*a/math.sqrt(1-e*e))*((e+math.cos(f))*math.cos(omega)*math.sin(inc) - math.sin(f)*math.sin(inc)*math.sin(omega))

    return rebound.Particle(
        m  = mass, 
        x = _x, y = _y, z = _z,
        vx = _vx, vy = _vy, vz = _vz)

#define TINY 1.0e-12
# from tools.c (changed star to ref to make the naming more general)
# converts cartesian elements to orbital elements.  ref should be a Particle
# structure with mass, r and v (e.g., the central body or the center of mass
# in Jacobi coordinates)

def p2orbit(G, p, ref):
    o = rebound.Orbit()
    mu = G*(p.m+ref.m)
    p.x -= ref.x
    p.y -= ref.y
    p.z -= ref.z
    p.vx -= ref.vx
    p.vy -= ref.vy
    p.vz -= ref.vz
    h0 = (p.y*p.vz - p.z*p.vy)                  # angular momentum vector
    h1 = (p.z*p.vx - p.x*p.vz)
    h2 = (p.x*p.vy - p.y*p.vx)
    o.h = numpy.sqrt ( h0*h0 + h1*h1 + h2*h2 )  # abs value of angular moment 
    v = numpy.sqrt ( p.vx*p.vx + p.vy*p.vy + p.vz*p.vz )
    o.r = numpy.sqrt ( p.x*p.x + p.y*p.y + p.z*p.z )
    vr = (p.x*p.vx + p.y*p.vy + p.z*p.vz)/o.r
    e0 = 1./mu*( (v*v-mu/o.r)*p.x - o.r*vr*p.vx )
    e1 = 1./mu*( (v*v-mu/o.r)*p.y - o.r*vr*p.vy )
    e2 = 1./mu*( (v*v-mu/o.r)*p.z - o.r*vr*p.vz )
    o.e = numpy.sqrt( e0*e0 + e1*e1 + e2*e2 )   # eccentricity
    o.a = -mu/( v*v - 2.*mu/o.r )               # semi major axis
    o.P = 2.*numpy.pi*numpy.sqrt( o.a*o.a*o.a/mu )  # period
    o.inc = numpy.arccos( h2/o.h )              # inclination (wrt xy-plane)
                                                # if pi/2<i<pi it's retrograde
    n0 = -h1                                    # node vector 
    n1 =  h0                                    # in xy plane => no z component     
    n = numpy.sqrt( n0*n0 + n1*n1 )
    er = p.x*e0 + p.y*e1 + p.z*e2               
    if (n<=1.e-30 or o.inc<=1.e-30):            # we are in the xy plane
        o.Omega=0.
        if (e1>=0.):
            o.omega=numpy.arccos(e0/o.e)
        else:
            o.omega = 2.*numpy.pi-numpy.arccos(e0/o.e)
    else:                                        
        if (e2>=0.):                            # omega=0 if perictr at asc node
            o.omega=numpy.arccos(( n0*e0 + n1*e1 )/(n*o.e))  
        else:
            o.omega=2.*numpy.pi-numpy.arccos(( n0*e0 + n1*e1 )/(n*o.e)) 
        if (n1>=0.):
            o.Omega = numpy.arccos(n0/n)
        else:
            o.Omega=2.*numpy.pi-numpy.arccos(n0/n)  # Omega=longitude of asc node
                                                # taken in xy plane from x axis

    o.f = er/(o.e*o.r)
    ea = (1.-o.r/o.a)/o.e
    if (o.f>1. or o.f<-1.):                     # failsafe
        o.f = numpy.pi - numpy.pi * o.f
        ea  = numpy.pi - numpy.pi * ea
    else:
        o.f = numpy.arccos(o.f)                 # true anom=0 if obj at perictr
        ea  = numpy.arccos(ea)                  # eccentric anomaly
    
    if (vr<0.):
        o.f=2.*numpy.pi-o.f    
        ea =2.*numpy.pi-ea
    
    o.l = ea -o.e*numpy.sin(ea) + o.omega+ o.Omega  # mean longitude
    if (o.e<=1.e-10):                           # circular orbit
        o.omega=0.
        o.f=0.                                  # f has no meaning
        o.l=0.

    return o

if __name__ == '__main__':
    sun = rebound.Particle( m=1.,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
    G = 1.
    params = (G,sun,0.,1.759, 0.851, 0.882, 5.852, 6.139, 5.445,True)
    p = add_planet_3D(*params)
    o = p2orbit(G,p,sun)