import numpy as np
import scipy.optimize
import sys
sys.path.append("../build/lib")
import libflow
import grid_utils, hickshenne
import yaml

class CG_Design:
    def __init__(self, inputs, airfoil0, cp_desired):
        assert(airfoil0.shape[0] == cp_desired.size+1)
        self.inputs     = inputs
        self.cp_desired = cp_desired
        self.airfoil0   = airfoil0
        self.dvars      = np.zeros((4,3))
        self.dvars[0,:] = [0.250,  0.5,   0.750] # lower loc
        self.dvars[1,:] = [0.250,  0.5,   0.750] # upper loc
        self.euler      = None
        self.adjoint    = None
        self.case       = None
        self.grid       = None
        self.xy         = None

    def _cost(self):
        if(self.euler is None):
            print "Trying to get cost of missing Euler object"
            raise
        cp1   = self.euler.pressure()
        cp1   = cp1[:,1] # only pressure, not x
        M     = self.inputs['mach']
        p_inf = 1.0 / 1.4
        dynp  = 0.5 * 1.0 * M * M
        pd    = self.cp_desired*dynp + p_inf
        p     = cp1*dynp + p_inf
        return np.float64(np.sum(0.5*(p-pd)**2))

    def _solve_euler(self):
        airfoil = hickshenne.perturb(self.airfoil0, self.dvars)
        mg  = libflow.MeshGen(airfoil, self.inputs['ktot'], 0.001)
        print "doing possion0 with dvars: ", self.dvars
        mg.poisson(600)
        xy  = mg.get_mesh()
        self.xy = xy
        del mg
        grid  = libflow.Grid(xy, 1)
        euler = libflow.Euler(grid, yaml.dump(self.inputs))
        euler.go()
        self.grid  = grid
        self.euler = euler

    def _solve_adjoint(self):
        if(self.euler is None):
            print "Trying to get adjoint of missing Euler object"
            raise
        adjoint  = libflow.Adjoint(self.euler)
        # adjoint  = libflow.ADadj(self.euler)
        self.adjoint = adjoint
        print "adjoint saved?"
        adjoint.init(self.cp_desired)
        adjoint.take_steps(10000)

    def _hash(self, dvars):
        return hash(str(dvars))
        
    def do_euler(self,dvars):
        key = self._hash(dvars)
        if self.case != key:
            # we haven't yet run this case. do it
            self.case  = key
            self.dvars = dvars
            self._solve_euler()
        cst = self._cost()
        print "returning cost: ", cst
        return cst

    def do_adjoint(self, dvars, eps=1e-8):
        key = self._hash(dvars)
        #
        # Check that we have an euler solution
        # 
        if self.case != key:
            print "uh oh, need to do euler first!"
            # this euler case hasn't been run, do this first
            self.do_euler(dvars)
        #
        # Solve the adjoint equations
        #
        self._solve_adjoint()
        #
        # Get sensitivity to each variable
        #
        sens  = np.zeros(self.dvars.size)
        # for i in range(sens.shape[0]/2, sens.shape[0]):
        for i in range(sens.shape[0]/2, sens.shape[0]/2+1):
            dvars = self.dvars.flatten() # flat copy of array
            dvars[i] += eps
            # the new airfoil:
            tmp_airfoil = hickshenne.perturb(self.airfoil0, dvars.reshape(self.dvars.shape))
            # the new mesh:
            mg     = libflow.MeshGen(tmp_airfoil, self.inputs['ktot'], 0.001)
            print "doing possion1 with dvars: ",  dvars.reshape(self.dvars.shape)
            mg.poisson(600)
            tmp_xy = mg.get_mesh()
            # the delta mesh
            d_xy = tmp_xy - self.xy
            sens[i] = self.adjoint.sens_xd(d_xy)/eps
        print "returning sensitivities: ", sens.reshape(self.dvars.shape)
        return sens.reshape(self.dvars.shape)




           
