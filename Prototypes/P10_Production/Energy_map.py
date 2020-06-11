import sys
import VIS_col as VC
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from md_tools import *


def get_bin2d(mins, deltas, value, size = 10):
    bins = np.floor(size*(value-mins)/deltas)
    bin_no = size * bins[0] + bins[1]
    return bin_no.astype(int)


saves = load("e_map.pickle")
if "a" not in saves:
    saves["a"] = load()
    save(saves, "e_map.pickle")
if "values" not in saves:
    values = np.array([]).reshape(0, len(saves["a"].startVIS.get_CVs()))
    for string in saves["a"].strings:
        for bead in string.SO:
            ori = bead.origin.get_CVs()
            for traj in bead.trajs:
                val = traj.get_CVs()
                rows = np.array([ori, val])
                values = np.concatenate((values,rows), axis = 0)
    saves["values"] = values
    save(saves, "e_map.pickle")
values = saves["values"]
mins = np.amin(values, axis = 0)
maxs = np.amax(values, axis = 0)
print(mins,maxs)
print(values.shape)
deltas = (maxs - mins)*1.0000001
size = 3


t_matrix = np.zeros((size**2,size**2))
for row in range(0,values.shape[0],2):
    i = get_bin2d(mins, deltas, values[row, :],size)
    j = get_bin2d(mins, deltas, values[row + 1, :],size)
    t_matrix[i, j] += 1

#print(t_matrix)
#input("pause")
v, vec = np.linalg.eig(t_matrix)
#print(vec)
#print(v)
#input("pause")
p = vec[0,:]
T = 300
kB = sp.constants.Boltzmann
E = -kB * T * np.log(p)
print(E.shape)

EM = np.reshape(E, (size,size))
print(E.shape, EM.shape)
x = np.linspace(mins[0], mins[0] + deltas[0], size)
y = np.linspace(mins[1], mins[1] + deltas[1], size)

print(x.shape, y.shape)
X,Y = np.meshgrid(x,y)
print(X.shape, Y.shape)
#np.nanmin(a[a != -np.inf])
levels = np.linspace(np.nanmin(E[E != -np.inf]),
                     np.nanmax(E[E != np.inf]),
                     size)
plt.contour(X,
            Y,
            EM,
            levels = levels,
            colors = ["red"])
plt.show()
print(levels)
print(EM)
