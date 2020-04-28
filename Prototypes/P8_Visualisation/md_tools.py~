import gmxapi as gmx
import numpy as np
import scipy.interpolate as interp
import scipy.integrate as integr
import VIS
import os
import pickle
import shutil

CV_RES = {"define": "-DPOSRES_CV"}
ITP_HEADERS = {}
ITP_WRITE = {}


def append_dict(main, extra):
    for key in extra:
        if key in ["include"] and key in main:
            main[key].union(extra[key])
        main[key] = extra[key]



def append_mdp(ordered_keys, all_parameters, name, all_keys, override = False):
    keys, parameters = read_mdp(name)
    for key in keys:
        if key in ["include"]:
            if key not in all_keys:
                all_keys.add(key)
                ordered_keys.append(key)
                all_parameters[key] = parameters[key]
            else:
                all_parameters[key].union(parameters[key])
        elif key not in all_keys or all_paremters[key] == "" or override:
            if override and key in all_keys and  all_paremters[key]!="":
                ordered_keys.remove(key)
            all_keys.add(key)
            ordered_keys.append(key)
            all_parameters[key] = parameters[key]
        else:
            raise KeyError("Conflicting .mdp file parameters")

def calc_splines(CVs):
    polys = []
    t = np.array(list(range(CVs.shape[0])))
    for i in range(CVs.shape[1]):
        polys.append(interp.CubicSpline(t,CVs[:,i],bc_type = 'natural'))
    return polys

def create_string(start, end, intermediaries, delta, saves):
    starter = VIS.VIS.fresh_copy(start)
    n = len(intermediaries)
    if "path" in saves:
        cv_span = saves["cv_span"]
        path = saves["path"]
    else:
        cv_span={}
        for i in range(len(delta)):
            cv_span["pull_coord" + str(i+1) +"_rate"] = delta[i]/500
        starter.steered(new_parameters = cv_span)
        path = starter.split_traj()
        saves["cv_span"] = cv_span
        saves["path"] = path
        save(saves)
    if "cv_traj" in saves:
        cv_traj = saves["cv_traj"]
    else:
        cv_traj = get_angles(path = path)
        saves["cv_traj"] = cv_traj
        save(saves)
    n_traj, n_targ = normalise_lin(cv_traj, intermediaries, starter, delta)
    indexes =  find_matches(n_traj, n_targ)
    #print(indexes)
    VIS_collection = []
    for i in indexes:
        VIS_collection.append(VIS.VIS(path+"conf"+str(i)+".gro"))
    return VIS_collection

def create_VIS (file, CV_values=None):
    vis = VIS.VIS(file)
    print("Start angle create_VIS: ", file, vis.get_CVs(), sep="\n" )
    vis.EM(restrain_cv = True)
    print("Start angle create_VIS: ", file, vis.get_CVs(), sep="\n" )
    # prepare VIS
    """ #
    vis.box()
    vis.solvate()
    vis.ions()
    vis.EM()
    vis.nvt()
    vis.npt()
    vis.run()
    """
    return vis

def denormalise(CVs, mins, deltas):
    denormed = np.zeros(CVs.shape)
    for i in range(len(deltas)):
        denormed[:, i] = CVs[:, i] * deltas[i] + mins[i]
    return denormed

def find_matches(trajectories, targets):
    distances = np.ones(targets.shape[0]*len(targets))
    indexes = [None]*len(targets)
    j=0
    for i in range(len(trajectories)):
        d1 = np.sum((trajectories[i, :] - targets[j, :])**2)
        d2 = np.sum((trajectories[i, :] - targets[j+1, :])**2)
        d3 = np.sum((trajectories[i, :] - targets[j+2, :])**2)
        if distances[j] > d1:
            distances[j] = d1
            indexes[j] = i
        if distances[j+1] > d2:
            distances[j+1] = d2
            indexes[j+1] = i
        if distances[j+2] > d3:
            distances[j+2] = d3
            indexes[j+2] = i
        if d3 < d2 and j < len(targets) -3:
            j += 1
    return indexes

def get_angle(topology_file, index_file):
    angler = gmx.commandline_operation(executable = "gmx",
                                       arguments = ["gangle",
                                                    "-g1", "dihedral",
                                                    "-group1", "dihedrals"],
                                       input_files = {"-s": topology_file,
                                                      "-n": index_file},
                                       output_files = {"-oall": "temp.xvg"})
    angler.run()
    err = angler.output.erroroutput.result()
    if err != "":
        print("get_angle:\n", err)
    result = open("temp.xvg").read().split("\n")[-2]
    angles = result.split()[1:]
    #print(angles)
    angles = [float(angle) for angle in angles]
    os.remove("temp.xvg")
    return angles

# default local and 501 steps
def get_angles(path = "", steps = 501):
    angles = []
    for i in range(steps):
        angles.append(get_angle(path+"conf" + str(i) + ".gro", "test.ndx"))
    return(angles)

def get_cartesian(topology_file, atom_nos = None):
    # 9 and 15
    if atom_nos is None:
        atom_nos = [9, 15]
    file = open(topology_file)
    lines = file.read().split("\n")
    atoms = [lines[atom_no+1] for atom_no in atom_nos]
    coords = [[float(coord) for coord in atom.split()[-3:]] for atom in atoms]
    return coords

def get_extremes(CVs):
    mins = np.min(CVs, axis = 0)
    maxs = np.max(CVs, axis = 0)
    deltas = maxs - mins
    return mins, deltas

def itp_write_dihedral(indexes, value, file):
    row = ["{:10}".format(i) for i in indexes]
    row = "".join(row)
    row += "{:10}{:10}   0   1".format(1,value)+"\n"
    file.write(row)

def itp_write_distance(indexes, vlaue, file):
    pass

def linear_interpolation(start, end, parts = 10):
    states = [None]*(parts-1)
    delta = []
    for i in range(len(start)):
        delta.append(end[i]-start[i])
        if delta[i] < -180:
            delta[i] =  360 + delta[i]
        elif delta[i] > 180:
            delta[i] = delta[i] - 360
    for i in range(1, parts):
        state = []
        for j in range(len(start)):
            state.append(start[j]+delta[j]*i/parts)
        states[i-1]=state
    return states, delta

def load(file_name = "debug.pickle"):
    try:
        return pickle.load(open(file_name,"rb"))
    except FileNotFoundError:
        return {}

def mdp_create(file_name, new_parameters = None, old_file = ""):      # for the .mdp-file
    if new_parameters is None:
        new_parameters = {}
    keys, parameters = [],{}
    if old_file != "":
        keys, parameters = read_mdp(old_file)
    if new_parameters != {}:
        keys.append("\n;Custom_parameters\n")
        parameters["\n;Custom_parameters\n"] = ""
    for key in new_parameters.keys():
        if key in ["include"]:
            if key not in keys:
                keys.insert(0, key)
                parameters[key] = new_parameters[key]
            else:
                parameters[key].union(new_parameters[key])
        elif key not in parameters:
            keys.append(key)
            parameters[key] = new_parameters[key]
        else:
            parameters[key] = new_parameters[key]
    file = open(file_name,"w")
    for key in keys:
        if key not in ["include"]:
            value = parameters[key]
            if value == "":
                file.write(key + "\n")
            else:
                file.write("{:<25} = {}\n".format(key,value))
        else:
            for value in parameters[key]:
                file.write("{:<25} = {}\n".format(key,value))

def normalise(CVs, mins, deltas):
    normed = np.zeros(CVs.shape)
    for i in range(len(deltas)):
        normed[:, i] = (CVs[:, i] - mins[i]) / deltas[i]
    return normed

def normalise_lin(cv_traj, intermediaries, start, delta):
    # Initial normalisation before first string is created,
    # specialised for normalising those values
    traj = np.array(cv_traj)
    targets = np.array(intermediaries)
    start = np.array(start.get_CVs())
    delta = np.array(delta)
    print ("traj", traj,"delta", delta,"start", start, sep="\n")
    for i in range(len(delta)):
        print(i)
        traj[:, i] = (traj[:, i] - start[i]) / delta[i]
        targets[:, i] = (targets[:, i] - start[i]) / delta[i]
    return traj, targets

def read_mdp(file_name):
    file =  open(file_name)
    keys=[]
    parameters = {}
    for row in file.read().split("\n"):
        if "=" not in row.split(";")[0]:
            key = row
            value = ""
        else:
            row = row.split("=")
            key = row[0].strip()
            value = row[1].strip()
        if key not in ["include", "define"]:
            keys.append(key)
            parameters[key]=value
        else:
            if key not in keys:
                keys.append(key)
                parameters[key]=set()
            parameters[key].add(value)
    file.close()
    return keys, parameters

def read_multi_mdp(file_names, override = False):
    all_keys = set()
    ordered_keys = []
    all_parameters = {}
    for name in file_names:
        append_mdp(ordered_keys, all_parameters, name, all_keys, override)
    return ordered_keys, all_parameters

def reparam_norm_cvs(splines, n):
    dists = []
    for i in range(n-1):
        dists.append(integr.quad(lambda x: np.sqrt(np.sum(np.array([spline(x, 1)**2 for spline in splines]))), i, i+1)[0])
    dists = np.array(dists)
    full_dist = np.sum(dists)
    reparam = np.zeros((n,len(splines)))
    so_far = 0
    index = 0
    reparam[0,:] = np.array([spline(0) for spline in splines])
    for i in range(1, n-1):
        target = i/(n-1)*full_dist
        while so_far + dists[index] < target:
            so_far += dists[index]
            index += 1
        fraction = (target - so_far) / dists[index]
        reparam[i,:] = np.array([spline(index+fraction) for spline in splines])
    reparam[n-1,:] = np.array([spline(n-1) for spline in splines])
    return reparam

def reparameterise(CVs):
    mins, deltas = get_extremes(CVs)
    norm_CVs = normalise(CVs, mins, deltas)
    splines = calc_splines(norm_CVs)
    new_norm_CVs = reparam_norm_cvs(splines, CVs.shape[0])
    new_CVs = denormalise(new_norm_CVs, mins, deltas)
    return new_CVs

def save(data, file_name = "debug.pickle"):
    pickle.dump(data,open(file_name, "wb"))

def update_topol_file(file_name = "topol.top"):
    name_parts = file_name.split(".")
    new_name = ".".join(name_parts[:-1])+"_orig.top"
    try:
        open(new_name)
    except FileNotFoundError:
        shutil.copyfile(file_name, new_name)
        new_file = open(file_name, "a")
        new_file.write('''
#ifdef POSRES_CV
#include "cv_restraints.itp"
#endif
''')
        new_file.close()


ITP_WRITE["dihedral"] = itp_write_dihedral
ITP_WRITE["distance"] = itp_write_distance
ITP_HEADERS["dihedral"] = "; ai   aj   ak   al  type  phi  dphi  kfac"
if __name__ == "__main__":
    mdp_create("simple.mdp")
