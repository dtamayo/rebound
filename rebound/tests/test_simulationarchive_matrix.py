import rebound
import unittest
import warnings

class TestSimulationArchiveMatrix(unittest.TestCase):
    pass

def runSimulation(tmax=40., restart=False, interval=None, safe_mode=True, integrator="ias15",G=1., testparticle=0):
    if restart:
        sim = rebound.Simulation.from_archive("test.bin")
    else:
        sim = rebound.Simulation()
        sim.G = G
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        if safe_mode==False:
            sim.ri_whfast.safe_mode = 1
            sim.ri_mercurius.safe_mode = 1
        if testparticle>0:
            if testparticle==1:
                sim.testparticle_type=0
            if testparticle==2:
                sim.testparticle_type=1
            sim.add(m=1e-4,a=1.2,e=0.04,omega=0.21,M=1.41,inc=0.21,Omega=1.1)
            sim.N_active = sim.N-1 # one test particle
        if interval:
            sim.automateSimulationArchive("test.bin", interval, deletefile=True)
    sim.integrate(tmax,exact_finish_time=0)
    return sim

def compareSim(test,sim1,sim2):
    test.assertEqual(sim1.N,sim2.N)
    test.assertEqual(sim1.N_active,sim2.N_active)
    test.assertEqual(sim1.N_var,sim2.N_var)
    test.assertEqual(sim1.t,sim2.t)
    test.assertEqual(sim1.G,sim2.G)
    for i in range(sim1.N):
        test.assertEqual(sim1.particles[i].r,sim2.particles[i].r)
        test.assertEqual(sim1.particles[i].m,sim2.particles[i].m)
        test.assertEqual(sim1.particles[i].x,sim2.particles[i].x)
        test.assertEqual(sim1.particles[i].y,sim2.particles[i].y)
        test.assertEqual(sim1.particles[i].z,sim2.particles[i].z)
        test.assertEqual(sim1.particles[i].vx,sim2.particles[i].vx)
        test.assertEqual(sim1.particles[i].vy,sim2.particles[i].vy)
        test.assertEqual(sim1.particles[i].vz,sim2.particles[i].vz)

def create_test(params):
    def doTest(self): 
        runSimulation(40., restart=False, interval=10., **params)
        sim1 = runSimulation(80., restart=True, **params)
        sim2 = runSimulation(80., restart=False, **params)
        compareSim(self,sim1,sim2)
    return doTest


for integrator in ["ias15","whfast","leapfrog","janus","mercurius"]:
    for safe_mode in [True,False]:
        for G in [1.,0.9]:
            for testparticle in [0,1,2]: # no test particle, passive, semi-active
                params = {'safe_mode':safe_mode,
                        'integrator':integrator,
                        'G':G, 
                        'testparticle':testparticle}
                test_method = create_test(params)
                name = "test"
                for key in params:
                    name += "_"+key+":"+str(params[key])
                test_method.__name__ = name
                setattr(TestSimulationArchiveMatrix, name,test_method)

if __name__ == "__main__":
    unittest.main()