# distutils: language = c++
# distutils: extra_compile_args = -O3

####################################################

# Numpy and math
import numpy as np

# C++ Standard Library vector class
from libcpp.vector cimport vector
from libcpp cimport bool

####################################################

# Import the parts of the C++ KCTF class we need
cdef extern from "kctf.hpp":
    cdef cppclass _Tidal_Body "Tidal_Body":
        double m,P,theta0,phi0,R,rg,k2,dt
        _Tidal_Body() except +
    cdef cppclass _KCTF "KCTF":
        double mu1,mu2,beta1,beta2
        double x, J20_const,J21_const, m0,m1,m2, mindist,mine
        double R0,R1, C0,C1, k20,k21, K0c,K1c, gammac,degrc
        double J20c,J21c, alpha10c,alpha11c, alpha20c,alpha21c
        vector[double] e2,G0, y,dy
        _KCTF() except +
        _KCTF(_Tidal_Body,_Tidal_Body,double,double,double,double,double,double,double) except +
        void print_line()
        bool is_sane(bool)
        void evolve(double)

####################################################

# Python KCTF class
cdef class KCTF:
    
    # Contain a copy of the C++ class
    cdef _KCTF *thisptr
    def __cinit__(self):
        self.thisptr = new _KCTF()
    def __del__(self):
        del self.thisptr

    # Setup the system
    def setup(self, b0,b1,double m20,double a1,double e1m,double I,double w,double a2,double e2m):
        cdef _Tidal_Body b0c,b1c
        b0c.m,b0c.P,b0c.theta0,b0c.phi0,b0c.R,b0c.rg,b0c.k2,b0c.dt = (
                b0.m,b0.P,b0.theta0,b0.phi0,b0.R,b0.rg,b0.k2,b0.dt)
        b1c.m,b1c.P,b1c.theta0,b1c.phi0,b1c.R,b1c.rg,b1c.k2,b1c.dt = (
                b1.m,b1.P,b1.theta0,b1.phi0,b1.R,b1.rg,b1.k2,b1.dt)
        self.thisptr = new _KCTF(b0c,b1c,m20,a1,e1m,I,w,a2,e2m)

    # Make a pickle
    def __reduce__(self):
        s,d = self.thisptr, {}
        d['doubles'] = [ s.x, s.J20_const,s.J21_const, s.m0,s.m1,s.m2, s.mindist,s.mine,
                s.R0,s.R1, s.C0,s.C1, s.k20,s.k21, s.K0c,s.K1c, s.gammac,s.degrc,
                s.J20c,s.J21c, s.alpha10c,s.alpha11c, s.alpha20c,s.alpha21c,
                s.mu1,s.mu2,s.beta1,s.beta2 ]
        d['e2'] = np.array(s.e2)
        d['G0'] = np.array(s.G0)
        d['y'] = np.array(s.y)
        d['dy'] = np.array(s.dy)
        return (KCTF, (), d)
    
    # Get out of a pickle
    def __setstate__(self, d):
        s = self.thisptr
        s.x, s.J20_const,s.J21_const, s.m0,s.m1,s.m2, s.mindist,s.mine = d['doubles'][:8]
        s.R0,s.R1, s.C0,s.C1, s.k20,s.k21, s.K0c,s.K1c, s.gammac,s.degrc = d['doubles'][8:18]
        s.J20c,s.J21c, s.alpha10c,s.alpha11c, s.alpha20c,s.alpha21c = d['doubles'][18:24]
        s.mu1,s.mu2,s.beta1,s.beta2 = d['doubles'][24:]
        s.e2,s.G0,s.y,s.dy = (d['e2'],d['G0'],d['y'],d['dy'])

    # Print Line
    def print_line(self):
        self.thisptr.print_line()

    # Check if anything is wrong
    def is_sane(self, verbose=False):
        cdef bool v2 = verbose
        return self.thisptr.is_sane(v2)

    # Integrate the system forward
    def evolve(self, double tgoal):
        self.thisptr.evolve(tgoal)

    # The time variable
    property x:
        def __get__(self): return self.thisptr.x
        def __set__(self,double v): self.thisptr.x = v

    # The minimum eccentricty
    property mine:
        def __get__(self): return self.thisptr.mine
        def __set__(self,double v): self.thisptr.mine = v

    # Allow the J2 of body 0 to be constant
    property J20_const:
        def __get__(self): return self.thisptr.J20_const
        def __set__(self,double v): self.thisptr.J20_const = v

    # Allow the J2 of body 1 to be constant
    property J21_const:
        def __get__(self): return self.thisptr.J21_const
        def __set__(self,double v): self.thisptr.J21_const = v

    '''# Evolve system
    def evolve(self, double tgoal, precision=None):
        if tgoal==self.t:
            return 0
        elif precision==None: 
            return self.thisptr.evolve_rkn(tgoal,1e-12)
        else:
            return self.thisptr.evolve_rkn(tgoal,float(precision))'''
