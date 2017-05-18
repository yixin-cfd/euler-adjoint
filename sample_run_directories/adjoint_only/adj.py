import numpy as np
import sys
sys.path.append("/home/dylan/work/archive/adjoint-euler/build/lib/")
sys.path.append("/home/dylan/work/archive/adjoint-euler/python")
import libflow
import grid_utils, euler_utils, hickshenne, cg_design
import yaml
from matplotlib import pyplot as plt
import naca
import scipy.optimize
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
x, y         = naca.naca4('0012', half, True, True)
x, y         = x[::-1], y[::-1]
airfoil0      = np.zeros((jtot,2))
airfoil0[:,0] = x
airfoil0[:,1] = y

# -----------------------------------------------------------------
# Get the CFD inputs
#
inputs = euler_utils.read_inputs('input.yaml')
inputs['ktot'] = ktot

# -----------------------------------------------------------------
# Read Desired distribution
#
desired = np.loadtxt('desired.dat', usecols=(1,)) # only pressure

# ----------------------------------------------------------------
# Design Class Prep
#
DC = cg_design.CG_Design(inputs, airfoil0, desired)

dvars = np.loadtxt("start.dat")

cost       = DC.do_euler(dvars)
sens       = DC.do_adjoint(dvars)
sens       = sens/cg_design.CFACTOR

print "Sensitivity: ", sens


# guess = np.array([0.001])

# opts = {"disp"       : True,
#         "return_all" : True,
#         "maxiter"    : 2,
# }
# result = scipy.optimize.minimize(f, guess, args=(DC,), method="cg", jac=df, options=opts)
# print result
