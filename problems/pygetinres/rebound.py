from ctypes import *
# Try to load libias15 from the obvioud places it could be in.

try:
    xrange          # If there's a name error we're using python 3.x
except NameError:
    xrange = range  # xrange = range in python 3.x

try:
    libias15 = CDLL('../../shared/libias15.so', RTLD_GLOBAL)
except:
    try:
        libias15 = CDLL('../shared/libias15.so', RTLD_GLOBAL)
    except:
        try:
            libias15 = CDLL('shared/libias15.so', RTLD_GLOBAL)
        except:
            print("Cannot find library 'libias15.so'. Check path set in 'rebound.py'.")
            raise

# Defines the same data structure as in particle.h
class Particle(Structure):
    _fields_ = [("x", c_double),
                ("y", c_double),
                ("z", c_double),
                ("vx", c_double),
                ("vy", c_double),
                ("vz", c_double),
                ("ax", c_double),
                ("ay", c_double),
                ("az", c_double),
                ("m", c_double) ]
    
# Defines the same data structure as in tools.h
class Orbit():
    def __init__(self):
        self.a      =   None    # semimajor axis
        self.r      =   None    # radial distance from reference
        self.h      =   None    # angular momentum 
        self.P      =   None    # orbital period
        self.l      =   None    # mean longitude = Omega + omega + M
        self.e      =   None    # eccentricity
        self.inc    =   None    # inclination
        self.Omega  =   None    # longitude of ascending node
        self.omega  =   None    # argument of perihelion
        self.f      =   None    # true anomaly

# Set function pointer for additional forces

AFF = CFUNCTYPE(None)
fp = None
def set_additional_forces(func):
    global fp  # keep references
    fp = AFF(func)
    libias15.set_additional_forces(func)

def init_damping_forces():
    libias15.init_damping_forces()
    
# Setter/getter of parameters and constants
def set_G(G):
    c_double.in_dll(libias15, "G").value = G

def get_G():
    return c_double.in_dll(libias15, "G").value

def set_dt(dt):
    c_double.in_dll(libias15, "dt").value = dt

def get_dt():
    return c_double.in_dll(libias15, "dt").value

def set_t(t):
    c_double.in_dll(libias15, "t").value = t

def set_min_dt(t):
    c_double.in_dll(libias15, "integrator_min_dt").value = t

def get_t():
    return c_double.in_dll(libias15, "t").value

def megno_init(delta):
    libias15.integrator_megno_init(c_double(delta))

def get_megno():
    libias15.integrator_megno.restype = c_double
    return libias15.integrator_megno()

def get_N():
    return c_int.in_dll(libias15,"N").value 

def add_migration(tau_a):
    arr = (c_double * get_N())(*tau_a)
    libias15.add_migration(byref(arr))
    
def add_e_damping(tau_e):
    arr = (c_double * get_N())(*tau_e)
    libias15.add_e_damping(byref(arr))
    
def add_i_damping(tau_i):
    arr = (c_double * get_N())(*tau_i)
    libias15.add_i_damping(byref(arr))
    
# Setter/getter of particle data
def set_particles(particles):
    c_int.in_dll(libias15,"N").value = len(particles)
    arr = (Particle * len(particles))(*particles)
    libias15.setp(byref(arr))

def particles_add(particles):
    particle_add(particles)

def particle_add(particles):
    if isinstance(particles,list):
        for particle in particles:
            libias15.particles_add(particle)
    else:
       libias15.particles_add(particles)

def particle_get(i):
    N = get_N() 
    if i>=N:
        return None
    getp = libias15.particle_get
    getp.restype = Particle
    _p = getp(c_int(i))
    return _p

def particles_get_array():
    N = get_N() 
    particles = []
    for i in xrange(0,N):
        particles.append(particle_get(i))
    return particles


def particles_get():
    N = c_int.in_dll(libias15,"N").value 
    getp = libias15.particles_get
    getp.restype = POINTER(Particle)
    return getp()


# Tools
def move_to_center_of_momentum():
    libias15.tools_move_to_center_of_momentum()

def reset():
    libias15.reset()

# Integration
def step():
    libias15.step()

def integrate(tmax):
    libias15.integrate(c_double(tmax))


