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

def get_bins1d(mins, deltas, value, size):
    bins = np.floor(size*(value-mins)/deltas)
    return bins.astype(int)

def normalise(matrix):
    b = np.sum(matrix, axis = 1)+1e-9
    return (matrix.T/b).T


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
np.save("swarm_coords.npy", values)
mins = np.amin(values, axis = 0)
maxs = np.amax(values, axis = 0)
print(mins,maxs)
print(values.shape)
deltas = (maxs - mins)*(1+1e-11)
size = 3
size2 = 11

t_matrix = np.zeros((size**2,size**2))
phi_t_matrix = np.zeros((size2,size2))
psi_t_matrix = np.zeros((size2,size2))
for row in range(0,values.shape[0],2):
    i = get_bin2d(mins, deltas, values[row, :],size)
    j = get_bin2d(mins, deltas, values[row + 1, :],size)
    t_matrix[i, j] += 1
    bins_i = get_bins1d(mins, deltas, values[row, :],size2)
    bins_j = get_bins1d(mins, deltas, values[row + 1, :],size2)
    phi_t_matrix[bins_i[0], bins_j[0]] += 1
    psi_t_matrix[bins_i[1], bins_j[1]] += 1

log("t_matrix", "emap.log")
log(str(t_matrix), "emap.log")
log("phi_t_matrix", "emap.log")
log(str(phi_t_matrix), "emap.log")
log("psi_t_matrix", "emap.log")
log(str(phi_t_matrix), "emap.log")
#input("pause")
t_matrix = normalise(t_matrix)
phi_t_matrix = normalise(phi_t_matrix)
psi_t_matrix = normalise(psi_t_matrix)

v, vec = np.linalg.eig(t_matrix)
phi_v, phi_vec = np.linalg.eig(phi_t_matrix)
psi_v, psi_vec = np.linalg.eig(psi_t_matrix)
#print(vec)
#print(v)
#input("pause")
p = vec[0,:]
phi_p = phi_vec[0, :]
psi_p = psi_vec[0, :]
T = 300
kB = sp.constants.Boltzmann
E = -kB * T * np.log(p)
E_phi = -kB * T * np.log(phi_p)
E_psi = -kB * T * np.log(psi_p)

print(E.shape)

EM = np.reshape(E, (size,size))
print(E.shape, EM.shape)
x = np.linspace(mins[0], mins[0] + deltas[0], size)
y = np.linspace(mins[1], mins[1] + deltas[1], size)

x2 = np.linspace(mins[0], mins[0] + deltas[0], size2)
y2 = np.linspace(mins[1], mins[1] + deltas[1], size2)

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
plt.figure()
plt.plot(x2, E_phi)
plt.xlabel("phi")

plt.figure()
plt.plot(y2, E_psi)
plt.xlabel("psi")

print(levels)
print(EM)
print(phi_v, phi_vec, sep="\n\n")
plt.show()
