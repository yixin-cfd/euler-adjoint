import numpy as np


def load_grid(filename):
    with open(filename, "r") as f:
        jtot, ktot = [int(x) for x in f.readline().split()]
    data = np.loadtxt(filename, skiprows=1)
    data = data.reshape((2,ktot,jtot))
    reordered = np.zeros((ktot, jtot, 2))
    reordered[:, :, 0] = data[0, :, :]
    reordered[:, :, 1] = data[1, :, :]
    return reordered
