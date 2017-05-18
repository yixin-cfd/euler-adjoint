import numpy as np
import sys
sys.path.append("/home/dylan/work/archive/2dflow/build/lib/")
sys.path.append("/home/dylan/work/archive/2dflow/python")
import libflow
import grid_utils, euler_utils
import yaml
from matplotlib import pyplot as plt
import naca
print "ok"


# ----------------------------------------------------------------
# Airfoil Surface
#
ktot = 60
half = 93
x, y = naca.naca4('0012', half, True, True)
x = x[::-1]
y = y[::-1]
jtot = half*2 + 1

airfoil      = np.zeros((jtot,2))
airfoil[:,0] = x
airfoil[:,1] = y

# ----------------------------------------------------------------
# Mesh Generation
#
mg  = libflow.MeshGen(airfoil, ktot, 0.001)
mg.poisson(600)
xy  = mg.get_mesh()

# ----------------------------------------------------------------
# Convert to CFD Grid
#
nghost = 2
grid   = libflow.Grid(xy, nghost)

# Start CFD
inputs = euler_utils.read_inputs("input.yaml")
euler  = libflow.Euler(grid, yaml.dump(inputs))

euler.go()

euler.write_solution('solution.dat')

pdata = euler.pressure()

np.savetxt("pressure.dat", pdata)

del euler
del grid
