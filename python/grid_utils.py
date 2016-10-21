from matplotlib import pyplot as plt
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


def plot2d(filename):
    with open(filename, "r") as f:
        jtot, ktot = [int(x) for x in f.readline().split()]
    
    data = np.loadtxt(filename, skiprows=1)
    data = data.reshape((2, ktot, jtot))
    
    for k in range(ktot):
        plt.plot(data[0,k,:], data[1,k,:], '-b')
    for j in range(jtot):
        plt.plot(data[0,:,j], data[1,:,j], '-b')
    plt.show()
