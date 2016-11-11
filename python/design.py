import numpy as np
import sys
sys.path.append('../build/lib')
import libflow
import grid_utils, euler_utils
import yaml

def solve_p(airfoil,inputs):
    mg  = libflow.MeshGen(airfoil, inputs['ktot'], 5.0)
    mg.poisson(500)
    xy  = mg.get_mesh()
    del mg

    grid  = libflow.Grid(xy, 1)
    euler = libflow.Euler(grid, yaml.dump(inputs))
    euler.go()
    p = euler.pressure()
    return p

class Design:
    def __init__(self, obj, inputs):
        assert(obj.hasKey("desired_pressure"))
        assert(obj.hasKey("base_airfoil"))
        self.desired     = obj['desired_pressure']
        self.airfoils    = [obj['base_airfoil']]
        self.design_vars = obj['design_vars']
        
        mg  = libflow.MeshGen(airfoil, ktot, 5.0)
        mg.poisson(500)
        self.base_xy  = mg.get_mesh()
        del mg

        grid  = libflow.Grid(self.base_xy, 1)
        euler = libflow.Euler(grid, yaml.dump(inputs))
        euler.take_steps(1000)
        self.base_p = euler.pressure()
        del euler
        del grid

        

