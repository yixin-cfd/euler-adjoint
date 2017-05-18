import numpy as np
import sys, time
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

def f(dvars, *args):
    DC    = args[0]
    cost  = DC.do_euler(dvars.reshape((4,-1)))
    return cost

def df(dvars, *args):
    DC    = args[0]
    sens  = DC.do_adjoint(dvars.reshape((4,-1)))
    return np.ravel(sens)

guess = np.zeros((4,3))
guess[0,:] = [0.250,  0.5,   0.750] # lower loc
guess[1,:] = [0.250,  0.5,   0.750] # upper loc
guess[2,:] = [0.000,  0.0,   0.000] # lower mag
guess[3,:] = [0.000,  0.0,   0.000] # upper mag

opts = {"disp"       : True,
        "maxiter"    : 4,
}


t1 = time.time()

result = scipy.optimize.minimize(f, np.ravel(guess), args=(DC,), method="SLSQP", jac=df, options=opts)
print result

t2 = time.time()
diff = t2-t1

print "time is %f seconds"%(diff)

