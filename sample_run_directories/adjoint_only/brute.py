import numpy as np
import sys
sys.path.append("/home/dylan/work/archive/2dflow/build/lib/")
sys.path.append("/home/dylan/work/archive/2dflow/python")
import libflow
import grid_utils, euler_utils, hickshenne, design
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

# starting variables
dvars = np.loadtxt('start.dat')
print dvars.shape
print dvars

# -----------------------------------------------------------------
# Create the Design class
#
obj = {
    "desired_pressure" : desired,
    "base_airfoil"     : airfoil0,
    "design_vars"      : dvars
    }

D         = design.Design(obj, inputs)

delta = 1.0e-8
s     = D.brute_sensitivities(delta)

print "Sensitivity: ", s
