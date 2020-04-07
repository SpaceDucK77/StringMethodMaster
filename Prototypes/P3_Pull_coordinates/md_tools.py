import gmxapi as gmx
import numpy as np
import os


def mdp_create(file_name, new_parameters = None, old_file = ""):      # for the .mdp-file
    if new_parameters is None:
        new_parameters = {}
    keys, parameters = [],{}
    if old_file != "":
        keys, parameters = read_mdp(old_file)
    orig_no = len(keys)
    for key in new_parameters.keys():
        if key not in parameters:
            keys.append(key)
        parameters[key] = new_parameters[key]
    file = open(file_name,"w")
    for par_count in range(len(keys)):
        if par_count == orig_no:
            file.write("\n;Custom_parameters\n")
        key = keys[par_count]
        value = parameters[key]
        if value == "":
            file.write(key + "\n")
        else:
            file.write("{:25} = {}\n".format(key,value))


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
        keys.append(key)
        parameters[key]=value
    file.close()
    return keys, parameters

def read_multi_mdp(file_names, override = False):
    all_keys = set()
    ordered_keys = []
    all_parameters = {}
    for name in file_names:
        append_mdp(ordered_keys, all_parameters, name, override)
    return ordered_keys, all_parameters

def append_dict(main, extra):
    for key in extra:
        main[key] = extra[key]



def append_mdp(ordered_keys, all_parameters, name, override = False):
    keys, parameters = read_mdp(name)
    for key in keys:
        if key not in all_keys or all_paremters[key] == "" or override:
            if override and key in all_keys and  all_paremters[key]!="":
                ordered_keys.remove(key)
            all_keys.add(key)
            ordered_keys.append(key)
            all_parameters[key] = parameters[key]
        else:
            raise KeyError("Conflicting .mdp file parameters")

def normalise(cv_traj, intermediaries, start, delta):
    # normed = normalise(cv_traj, intermediaries, start, end, delta)
    traj = np.array(cv_traj)
    targets = np.array(intermediaries)
    print(traj)
    print(targets)
    start = np.array(start.collection)
    delta = np.array(delta)
    for i in range(len(delta)):
        traj[i, :] = (traj[i, :] - start[i]) / delta[i]
        targets[i, :] = (targets[i, :] - start[i]) / delta[i]
    return traj, targets


def get_angle(topology_file, index_file):
    angler = gmx.commandline_operation(executable = "gmx",
                                       arguments = ["gangle",
                                                    "-g1", "dihedral",
                                                    "-group1", "dihedrals"],
                                       input_files = {"-s": topology_file,
                                                      "-n": index_file},
                                       output_files = {"-oall": "temp.xvg"})
    angler.run()
    result = open("temp.xvg").read().split("\n")[-2]
    angles = result.split()[1:]
    #print(angles)
    angles = [float(angle) for angle in angles]
    os.remove("temp.xvg")
    return angles


def get_cartesian(topology_file, atom_nos = None):
    # 9 and 15
    if atom_nos is None:
        atom_nos = [9, 15]
    file = open(topology_file)
    lines = file.read().split("\n")
    atoms = [lines[atom_no+1] for atom_no in atom_nos]
    coords = [[float(coord) for coord in atom.split()[-3:]] for atom in atoms]
    return coords

# default local and 501 steps
def get_angles(path = "", steps = 501):
    angles = []
    for i in range(steps):
        angles.append(get_angle(path+"conf" + str(i) + ".gro", "test.ndx"))
    return(angles)

def find_matches(trajectories, targets):
    print(trajectories,targets, sep="\n")
    distancses = np.ones_like(targets)*len(targets)
    indexes = [None]*len(targets)
    j=0
    for i in range(len(trajectories)):
        d1 = (trajectories[:, i] - targets[:, j])**2
        d2 = (trajectories[:, i] - targets[:, j+1])**2
        d3 = (trajectories[:, i] - targets[:, j+2])**2
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

if __name__ == "__main__":
    mdp_create("simple.mdp")
