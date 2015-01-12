#!/usr/bin/python

# Import the rebound module
import sys; sys.path.append('../')
import rebound
# Import other modules
import math

# mod2pi helper function.
TWOPI = 2.*math.pi
def mod2pi(f):
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

# Converts orbital elements to cartesian and returns a particle
# if Mean is set to True, the anomaly is taken as the mean anomaly
# if Mean is set to False, the anomaly is taken as the true anomaly
def add_planet_3D(G,com,mass,a,e,i,Omega,omega,anom,MEAN=True):
    if MEAN==False:
        f = anom
    else:
        M = anom
        E = M
        if e>0.8:
            E = np.pi
        F = E - e*math.sin(M) - M
        for i in xrange(100):
            E = E - F/(1.0-e*math.cos(E))
            F = E - e*math.sin(E) - M
            if math.fabs(F)<1e-16:
                break
        E = mod2pi(E)
        f = mod2pi(2.*math.atan(math.sqrt((1. + e)/(1. - e))*math.tan(0.5*E)))
    
    r = a*(1-e**2)/(1 + e*math.cos(f));
    
    # Murray & Dermott Eq. 2.122
    _x  = com.x + r*(math.cos(Omega)*math.cos(omega+f) - math.sin(Omega)*math.sin(omega+f)*math.cos(i));
    _y  = com.y + r*(math.sin(Omega)*math.cos(omega+f) + math.cos(Omega)*math.sin(omega+f)*math.cos(i));
    _z  = com.z + r*math.sin(omega+f)*math.sin(i);
    
    n = math.sqrt(G*(com.m+mass)/(a**3));
    
    # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    _vx = com.vx + (n*a/math.sqrt(1-e*e))*((e+math.cos(f))*(-math.cos(i)*math.cos(omega)*math.sin(Omega) - math.cos(Omega)*math.sin(omega)) - math.sin(f)*(math.cos(omega)*math.cos(Omega) - math.cos(i)*math.sin(omega)*math.sin(Omega)));
    _vy = com.vy + (n*a/math.sqrt(1-e*e))*((e+math.cos(f))*(math.cos(i)*math.cos(omega)*math.cos(Omega) - math.sin(omega)*math.sin(Omega)) - math.sin(f)*(math.cos(omega)*math.sin(Omega) + math.cos(i)*math.cos(Omega)*math.sin(omega)));
    _vz = com.vz + (n*a/math.sqrt(1-e*e))*((e+math.cos(f))*math.cos(omega)*math.sin(i) - math.sin(f)*math.sin(i)*math.sin(omega));

    return rebound.Particle(
        m  = mass, 
        x = _x, y = _y, z = _z,
        vx = _vx, vy = _vy, vz = _vz)

if __name__ == '__main__':
    sun = rebound.Particle( m=1.,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
    p = add_planet_3D(1.,sun,0.,1.,0.01,0.,0.,0.,0.,MEAN=False)
    print(p.x, p.y, p.z, p.vx, p.vy, p.vz)
