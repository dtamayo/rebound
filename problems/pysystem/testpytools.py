import pytools
import unittest
import rebound
import random
import numpy

def almost_equal_wrap_2pi(val1,val2, places):
    diff = val2-val1
    diff2 = diff + 2*numpy.pi if diff < 0 else diff - 2*numpy.pi
    return True if min(numpy.fabs(diff), numpy.fabs(diff2)) < 10**(-places) else False

class cartesian_to_orbital(unittest.TestCase):
    sun = rebound.Particle( m=1.,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
    G = 1.
    N_random_tests = 10 # num of random orbits to test in test_rand_r_to_orb_(f or M)
    
    def test_aew2pi(self):
        places=15
        cases = ((0.,10**(-16), True),
                 (0.,10**(-16), True),
                 (0.,10**(-14), False),
                 (0.,10**(-14), False),
                 (0.1,0.2, False))
        for case in cases:
            self.assertIs(almost_equal_wrap_2pi(case[0], 2*numpy.pi - case[1], places),case[2], '{}'.format(case))
            self.assertIs(almost_equal_wrap_2pi(2*numpy.pi - case[1], case[0], places),case[2], '{}'.format(case))
            
    def test_keplers_eq(self):
        '''test Kepler's equation'''
        zero_to_pi = numpy.linspace(0,numpy.pi,10,endpoint=True)
        pi_to_2pi = numpy.linspace(numpy.pi,2*numpy.pi,10,endpoint=False)
        e = numpy.asarray([0.,0.7,0.9,0.999])
        for ecc in e:
            # E & M match at 0,pi, so if M is in range 0<M<pi (or pi<M<2pi), so should E
            # Also check that Kepler's equation is satisfied to abs precision of 1e-15
            for M in zero_to_pi:
                E = pytools.get_E(ecc,M)
                err = E-ecc*numpy.sin(E)-M
                self.assertTrue(0.<= E <= numpy.pi, 
                                "E & M not in same half-plane: E=%.3f, M=%.3f"%(E,M))
                self.assertAlmostEqual(numpy.fabs(err),0.,places=14,msg=
                "Kepler's Eq not satisfied: e=%.3f, M=%.3f, abs. error=%.3e"%(ecc,M,err))
            for M in pi_to_2pi:
                E = pytools.get_E(ecc,M)
                err = E-ecc*numpy.sin(E)-M
                self.assertTrue(numpy.pi<= E <= 2*numpy.pi,
                                "E & M not in same half-plane: E=%.3f, M=%.3f"%(E,M))       
                self.assertAlmostEqual(numpy.fabs(err),0.,places=14,msg=
                "Kepler's Eq not satisfied: e=%.3f, M=%.3f, abs. error=%.3e"%(ecc,M,err))

    def test_r_to_orb_defaults(self):
        '''test conversion from orbital elements to cartesian and back 
        when not all orbital elements are passed to init_planet'''
        specified_cases = ((self.sun,1.,),
                       (self.sun,1.,0.01)
                       )
        results = ((1.,0.,0.,0.,0.,0.),
                   (1.,0.01,0.,0.,0.,0.))
        places=12
        for j in range(len(specified_cases)):
            p = pytools.init_planet(*specified_cases[j])
            o = pytools.p2orbit(p,self.sun)
            self.assertAlmostEqual(results[j][0],o.a,places=places)
            self.assertAlmostEqual(results[j][1],o.e,places=places)
            self.assertAlmostEqual(results[j][2],o.inc,places=places)
            self.assertIs(almost_equal_wrap_2pi(results[j][3],pytools.mod2pi(o.Omega),places), True, '{}'.format(specified_cases[j]))
            self.assertIs(almost_equal_wrap_2pi(results[j][4],pytools.mod2pi(o.omega),places), True, '{}'.format(specified_cases[j]))
            self.assertIs(almost_equal_wrap_2pi(results[j][5],o.f,places), True, '{}'.format(specified_cases[j]))
            
    def test_specified_r_to_orb(self):
        '''test conversion from orbital elements to cartesian and back with specified cases'''
        specified_cases = ((self.sun,1.,0.01,0.,0.,0.,0.,0.),
                       (self.sun,1.,0.999,2.,3.,3.,0.,0.),
                       (self.sun,1.759,0.851,0.882,1.287,1.728,5.445,0.),
                       (self.sun,42.,1.e-8,0.,0.,1.e-8,0.,0.)
                       )
        places=12
        for params in specified_cases:
            #print("%.3f, %.3f, %.3f, %.3f, %.3f, %.3f"%(params[3:9]))
            p = pytools.init_planet(*params)
            o = pytools.p2orbit(p,self.sun)
            self.assertAlmostEqual(params[1],o.a,places=places, msg='{}'.format(params))
            self.assertAlmostEqual(params[2],o.e,places=places, msg='{}'.format(params))
            self.assertAlmostEqual(params[5],o.inc,places=places, msg='{}'.format(params))
            self.assertIs(almost_equal_wrap_2pi(params[6],pytools.mod2pi(o.Omega),places), True, '{}'.format(params))
            self.assertIs(almost_equal_wrap_2pi(params[3],pytools.mod2pi(o.omega),places), True, '{}'.format(params))
            self.assertIs(almost_equal_wrap_2pi(params[4],o.f,places), True, '{}'.format(params))
            
    def test_rand_r_to_orb_f(self):    
        '''test conversion from orb. elements to cart. and back with random cases w/ true anom'''
        places=12
        for q in range(self.N_random_tests):
            _a=random.uniform(1.,2.)
            _e=random.uniform(0.,1.)
            _inc=random.uniform(0,numpy.pi)
            _Omega=random.uniform(0,2*numpy.pi)
            _omega=random.uniform(0,2*numpy.pi)
            _f=random.uniform(0,2*numpy.pi)
            #print("%.3f, %.3f, %.3f, %.3f, %.3f, %.3f"%(_a,_e,_inc,_Omega,_omega,_f))
            p = pytools.init_planet(self.sun,a=_a, e=_e, omega=_omega, f=_f, inc=_inc, 
                                    Omega=_Omega, mass=0.)
            o = pytools.p2orbit(p,self.sun)
            self.assertAlmostEqual(_a,o.a,places=12)
            self.assertAlmostEqual(_e,o.e,places=12)
            self.assertAlmostEqual(_inc,o.inc,places=12)
            self.assertIs(almost_equal_wrap_2pi(_Omega,pytools.mod2pi(o.Omega),places), True)
            self.assertIs(almost_equal_wrap_2pi(_omega,pytools.mod2pi(o.omega),places), True)
            self.assertIs(almost_equal_wrap_2pi(_f,pytools.mod2pi(o.f),places), True)
            
    def test_rand_r_to_orb_f(self):    
        '''test conversion from orb. elements to cart. and back with random cases w/ mean anom'''
        places=12
        for q in range(self.N_random_tests):
            _a=random.uniform(1.,2.)
            _e=random.uniform(0.,1.)
            _inc=random.uniform(0,numpy.pi)
            _Omega=random.uniform(0,2*numpy.pi)
            _omega=random.uniform(0,2*numpy.pi)
            _M=random.uniform(0,2*numpy.pi)
            #print("%.3f, %.3f, %.3f, %.3f, %.3f, %.3f"%(_a,_e,_inc,_Omega,_omega,_f))
            p = pytools.init_planet(self.sun,a=_a, e=_e, omega=_omega, inc=_inc, 
                                    Omega=_Omega, mass=0., M=_M)
            o = pytools.p2orbit(p,self.sun)
            self.assertAlmostEqual(_a,o.a,places=12)
            self.assertAlmostEqual(_e,o.e,places=12)
            self.assertAlmostEqual(_inc,o.inc,places=12)
            self.assertIs(almost_equal_wrap_2pi(_Omega,pytools.mod2pi(o.Omega),places), True)
            self.assertIs(almost_equal_wrap_2pi(_omega,pytools.mod2pi(o.omega),places), True)
            self.assertIs(almost_equal_wrap_2pi(_M,pytools.mod2pi(o.l-o.Omega-o.omega),places), True)

if __name__ == "__main__":
    unittest.main()