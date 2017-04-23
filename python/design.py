import numpy as np
import sys
sys.path.append('../build/lib')
import libflow
import grid_utils, euler_utils, hickshenne
import yaml

def solve_euler(airfoil,inputs):
    mg  = libflow.MeshGen(airfoil, inputs['ktot'], 0.001)
    mg.poisson(600)
    xy  = mg.get_mesh()
    del mg
    #
    grid  = libflow.Grid(xy, 1)
    euler = libflow.Euler(grid, yaml.dump(inputs))
    euler.go()
    return euler

def solve_p(airfoil,inputs):
    euler = solve_euler(airfoil,inputs)
    p = euler.pressure()
    del euler
    return p

def find_cost(p1, p2):
    return np.float64(np.sum(0.5*(p1-p2)**2))

class Design:
    def __init__(self, obj, inputs):
        assert("desired_pressure" in obj)
        assert("base_airfoil" in obj)
        self.desired      = obj['desired_pressure']
        self.design_vars  = obj['design_vars']
        self.inputs       = inputs
        self.base_airfoil = obj['base_airfoil']
        self.airfoils     = []
        self.pressures    = []
        self.dynp         = 0.5*inputs['mach']*inputs['mach']
        self.pinf         = 1.0/1.4
        # add our airfoil to the list and the pressure too
        airfoil1 = hickshenne.perturb(self.base_airfoil, self.design_vars)
        self.airfoils.append(airfoil1)
        pressure = solve_p(airfoil1, inputs)
        pressure = pressure*self.dynp + self.pinf
        self.desired = self.desired*self.dynp + self.pinf
        self.pressures.append(pressure)
        print "start cost: ", find_cost(pressure, self.desired)

    def brute_sensitivities(self, eps):
        s        = np.zeros(self.design_vars.size)
        pressure = self.pressures[-1]
        cost     = find_cost(pressure, self.desired)
        
        # for i in range(s.shape[0]/2, s.shape[0]):
        for i in range(s.shape[0]/2, s.shape[0]/2+1):
            dvars = self.design_vars.flatten() # return flat copy of array
            dvars[i] += eps
            airfoil1 = hickshenne.perturb(self.base_airfoil, dvars.reshape(self.design_vars.shape))
            new_pressure = solve_p(airfoil1, self.inputs)
            new_pressure = new_pressure*self.dynp + self.pinf
            new_cost     = find_cost(new_pressure, self.desired)
            s[i] = (new_cost - cost) / eps
            print "%d: (%e - %e) / %e = %e"%(i,new_cost, cost, eps, s[i])
        
        return s.reshape(self.design_vars.shape)

    def second_deriv(self, eps, i):
        s        = np.zeros(self.design_vars.size)
        pressure = self.pressures[-1]
        cost     = find_cost(pressure, self.desired)
        
        dvars = self.design_vars.flatten() # return flat copy of array
        dvars[i] += eps
        airfoil1 = hickshenne.perturb(self.base_airfoil, dvars.reshape(self.design_vars.shape))
        new_pressure = solve_p(airfoil1, self.inputs)
        cost1        = find_cost(new_pressure, self.desired)

        dvars = self.design_vars.flatten() # return flat copy of array
        dvars[i] -= eps
        airfoil1 = hickshenne.perturb(self.base_airfoil, dvars.reshape(self.design_vars.shape))
        new_pressure = solve_p(airfoil1, self.inputs)
        cost2        = find_cost(new_pressure, self.desired)

        d2 = (cost2 - 2*cost + cost1)/(eps*eps)
        print d2

        return d2

