import numpy as np
import sys
sys.path.append("/home/dylan/work/archive/adjoint-euler/build/lib/")
sys.path.append("/home/dylan/work/archive/adjoint-euler/python")
import libflow
import grid_utils, euler_utils, hickshenne, cg_design, design
import yaml
from matplotlib import pyplot as plt
import naca
print "ok"

# -----------------------------------------------------------------
# Inputs
# 
ktot = 60
jtot = 181

# -----------------------------------------------------------------
# Starting Airfoil Surface
#
if(jtot%2==0):
    quit('jtot must be odd')
half         = (jtot-1)/2
x, y         = naca.naca4('2312', half, True, True)
x, y         = x[::-1], y[::-1]
airfoil0      = np.zeros((jtot,2))
airfoil0[:,0] = x
airfoil0[:,1] = y

# -----------------------------------------------------------------
# Get the CFD inputs
#
inputs = euler_utils.read_inputs("input.yaml")
inputs['ktot'] = ktot

# -----------------------------------------------------------------
# Perturb to get an airfoil
#
dvars2 = np.zeros((4,3))
dvars2[0,:] = [ 0.250,  0.500,   0.750] # lower loc
dvars2[1,:] = [ 0.250,  0.500,   0.750] # upper loc
dvars2[2,:] = [ 0.000,  0.000,   0.000] # lower mag
dvars2[3,:] = [ 0.000,  0.000,   0.000] # upper mag
airfoil2   = hickshenne.perturb(airfoil0, dvars2)
pressure2  = design.solve_p(airfoil2, inputs)
np.savetxt("desired.dat", pressure2)

