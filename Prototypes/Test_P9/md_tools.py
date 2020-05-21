import gmxapi as gmx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as interp
import scipy.integrate as integr
import plotly.express as px
import VIS
import os
import pickle
import shutil

CV_RES = {"define": "-DPOSRES_CV"}
ITP_HEADERS = {}
ITP_WRITE = {}


def append_dict(main, extra, mkeys = None, ekeys = None):
    for key in extra:
        if key in ["include"] and key in main:
            main[key].union(extra[key])
        main[key] = extra[key]
    if mkeys != None:
        for key in ekeys:
            if key not in mkeys:
                mkeys.append(key)



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

def backup_file(file_name, copy):
    done = False
    backup_no = 1
    f1, f2 = file_name.rsplit(".",1)
    try:
        open(file_name)
    except FileNotFoundError:
        return
    while not done:
        bu_name = f1 + "_backup_" + str(backup_no) + "_." + f2
        try:
            open(bu_name)
            backup_no += 1
        except FileNotFoundError:
            if copy:
                shutil.copyfile(file_name, bu_name)
            else:
                shutil.move(file_name, bu_name)
            done = True

def calc_splines(CVs):
    polys = []
    t = np.array(list(range(CVs.shape[0])))
    for i in range(CVs.shape[1]):
        polys.append(interp.CubicSpline(t,CVs[:,i],bc_type = 'natural'))
    return polys

def create_panda(start,
                 end,
                 init_states,
                 iterations,
                 CV_names,
                 CV_index = ()):
    if CV_index == ():
        CV_index = tuple(range(len(CV_names) + 2))
    else:
        CV_index = tuple([0, 1] + [index + 2 for index in CV_index])
    set_length = (len(init_states) + 2)
    base = np.zeros(((len(iterations) + 1) * set_length, len(CV_names) + 2))
    base[:set_length, 0] = np.ones((set_length))
    base[:set_length, 1] = np.array(list(range(set_length))).T
    base[:set_length, 2:] = np.append([start],
                                      np.append(init_states,
                                                [end],
                                                axis = 0),
                                      axis = 0)
    for i in range(len(iterations)):
        j = i+1
        base[j * set_length : (j + 1) * set_length, 0] = (j + 1) * np.ones((set_length))
        base[j * set_length : (j + 1) * set_length, 1] = np.array(list(range(set_length))).T
        base[j * set_length : (j + 1) * set_length, 2:] = iterations[i].get_CVs()
    dict_base = {}
    names = ["Iter", "VISno"] + list(CV_names)
    for index in CV_index:
        dict_base[names[index]] = base[:,index]
    frame = pd.DataFrame(dict_base)
    frame = frame[names]
    return frame

def create_string(start, end, intermediaries, delta, saves, opposites):
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
        print("---- Steered No crash ----")
        path = starter.split_traj()
        saves["cv_span"] = cv_span
        saves["path"] = path
        #save(saves)
    if "cv_traj" in saves:
        cv_traj = saves["cv_traj"]
    else:
        cv_traj = get_CVs(path = path,
                          CVs = starter.CV_def,
                          index_file = starter.index_file)
        saves["cv_traj"] = cv_traj
        #save(saves)
    startCV = start.get_CV_coll()
    cv_traj2 =  np.concatenate((np.array(cv_traj),
                                np.array([startCV["dihedral"] + startCV["distance"]])))
    mins, deltas = get_extremes(cv_traj2, opposites)
    #print(mins, deltas, sep = "\n")
    n_dihedrals = len(startCV["dihedral"])
    n_traj = normalise(np.array(cv_traj), mins, deltas, n_dihedrals)
    n_targ = normalise(np.array(intermediaries), mins, deltas, n_dihedrals)
    indexes =  find_matches(n_traj, n_targ)
    VIS_collection = []
    for i in indexes: # FIXED? CALL TO CONSTRUCTOR
        bead = VIS.VIS.fresh_copy(starter)
        bead.configuration_file = path+"conf"+str(i)+".gro"
        VIS_collection.append(bead)
    return VIS_collection

def create_VIS (c_file,
                CV_definitions = None,
                topology_file = "topol.top",
                index_file = "index.ndx",
                pull_groups = None,
                solvated = 0):
    vis = VIS.VIS(c_file,
                  CV_definitions = CV_definitions,
                  topology_file = topology_file,
                  index_file = index_file,
                  pull_groups = pull_groups,
                  solvated = solvated)
    #vis.EM(restrain_cv = True)
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

def delta_angle(angle, origin):
    angle -= origin
    if angle < -180:
        angle += 360
    elif angle > 180:
        angle -= 360
    return angle

def delta_angles(angle_array, origin):
    angle_array = angle_array - origin
    for i in range(len(angle_array)):
        angle = angle_array[i]
        if angle < -180:
            angle += 360
        elif angle > 180:
            angle -= 360
        angle_array[i] = angle
    return angle_array

def denormalise(newCVs, oldCVs, mins, deltas, dih_no):
    denormed = np.zeros(newCVs.shape)
    drift = np.zeros(newCVs.shape)
    for i in range(len(deltas)):
        drift[:, i] = newCVs[:, i] * deltas[i]
        if i < dih_no:
            denormed[:, i] = delta_angles(drift[:, i], -mins[i])
        else:
            denormed[:, i] = drifts[:, i] + mins[i]
    return denormed, drift

def dictionarise(text):
    rows = text.split("\n")
    result = {}
    for row in rows:
        no_comment = row.split("#")[0].strip()
        if len(no_comment) > 0:
            parts = no_comment.split("=")
            result[parts[0].strip()] = parts[1].strip()
    return result

def find_greatest_d_angle(angles):
    lo = 0
    hi = 0
    max_delta = 0
    angles2 = np.sort(angles)
    begin = 0
    for i in range(len(angles2)):
        d_angle = (360 + angles2[i] - angles2[i - 1]) % 360
        #print (360 - d_angle, d_angle, angles2[i] - angles2[i-1], angles2[i], angles2[i - 1], sep = "\t")
        if d_angle > 180:
            #input("pause")
            return angles2[i], 360 - d_angle
    raise ValueError("path not clear")



def find_matches(trajectories, targets):
    distances = np.ones(targets.shape[0]*(len(targets)**4))
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

def get_angle(conf_file, index_file):
    angler = gmx.commandline_operation(executable = "gmx",
                                       arguments = ["gangle",
                                                    "-g1", "dihedral",
                                                    "-group1", "dihedrals"],
                                       input_files = {"-s": conf_file,
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

def get_CVs(path = "",
            CVs = None,
            steps = 501,
            index_file = "index.ndx"):
    if CVs == None:
        raise Exception("CV dictionary needed to extract CVs from .gro files")
    CV_matrix = []
    for i in range(steps):
        conf = path+"conf" + str(i) + ".gro"
        CV_list = get_angle(conf, index_file)
        for dist in CVs["distances"]:
            CV_list += get_distance (conf, index_file, dist)
        CV_matrix.append(CV_list)
    return CV_matrix


def get_distance(conf_file, index_file, cv):
    # gmx distance -s em.gro -n indexB.ndx -select 'com of group "first_A" plus com of group "first_B"' -oall dist_test.xvg
    gr_sel = "com of group \"" + cv.pull_groups[0].get_name() +\
             "\" plus com of group \"" +cv.pull_groups[1].get_name() + "\""
    distancer = gmx.commandline_operation(executable = "gmx",
                                          arguments = ["distance",
                                                       "-select",
                                                       gr_sel],
                                          input_files = {"-s": conf_file,
                                                         "-n": index_file},
                                          output_files = {"-oall": "temp.xvg"})
    distancer.run()
    err = distancer.output.erroroutput.result()
    if err != "":
        print("get_distance:\n", err)
    result = open("temp.xvg").read().split("\n")[-2]
    dist = result.split()[-1]
    os.remove("temp.xvg")
    return [float(dist)]

'''
def get_cartesian(topology_file, atom_nos = None):
    # 9 and 15
    if atom_nos is None:
        atom_nos = [9, 15]
    file = open(topology_file)
    lines = file.read().split("\n")
    atoms = [lines[atom_no+1] for atom_no in atom_nos]
    coords = [[float(coord) for coord in atom.split()[-3:]] for atom in atoms]
    return coords
'''
def get_extremes(CVs, opposites):
    print(CVs.shape)
    mins = np.min(CVs, axis = 0)
    maxs = np.max(CVs, axis = 0)
    deltas = maxs - mins
    for i in range(len(opposites)):
        if opposites[i]:
            base, delta = find_greatest_d_angle(CVs[:, i])
            mins[i] = base
            deltas[i] = delta
    return mins, deltas

def itp_write_dihedral(indexes, value, file):
    row = ["{:10}".format(i) for i in indexes]
    row = "".join(row)
    row += "{:10}{:10}   0   1".format(1,value)+"\n"
    file.write(row)

def itp_write_distance(indexes, value, file):
    pass

def linear_interpolation(start, end, parts = 10, no_dih = -1):
    states = [None]*(parts-1)
    delta = []
    opposites = []
    for i in range(len(start)):
        delta.append(end[i]-start[i])
        opposites.append(False)
        if -1 == no_dih or i < no_dih:
            if delta[i] < -180:
                delta[i] =  360 + delta[i]
                opposites[i] = True
            elif delta[i] > 180:
                delta[i] = delta[i] - 360
                opposites[i] = True
    for i in range(1, parts):
        state = []
        for j in range(len(start)):
            deltaj = start[j]+delta[j]*i/parts
            if -1 == no_dih or j < no_dih:
                if deltaj < -180:
                    deltaj =  360 + deltaj
                elif deltaj > 180:
                    deltaj = deltaj - 360
            state.append(deltaj)
        states[i-1]=state
    return states, delta, opposites

def load(file_name = "debug.pickle"):
    try:
        return pickle.load(open(file_name,"rb"))
    except FileNotFoundError:
        return {}

def log(msg, file_name = "log.txt"):
    f = open(file_name, "a")
    f.write(msg + "\n")
    f.close()

def make_index(c_file, o_file):
    backup_file(o_file, copy = False)
    maker =  gmx.commandline_operation(executable = "gmx",
                                       arguments = ["make_ndx"],
                                       input_files = {"-f": c_file},
                                       output_files = {"-o": o_file},
                                       stdin = "q\n")
    maker.run()
    print("make_index:\n", maker.output.erroroutput.result())


def mdp_create(file_name,
               new_parameters = None,
               new_keys = None,
               old_file = ""):      # for the .mdp-file
    if new_parameters is None:
        new_parameters = {}
    keys, parameters = [],{}
    if old_file != "":
        keys, parameters = read_mdp(old_file)
    if new_parameters != {}:
        keys.append("\n;Custom_parameters\n")
        parameters["\n;Custom_parameters\n"] = ""
        if new_keys == None:
            new_keys = new_parameters.keys()
    for key in new_keys:
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

def normalise(CVs, mins, deltas, CV_angles):
    normed = np.zeros(CVs.shape)
    for i in range(len(deltas)):
        if i < CV_angles:
            normed[:, i] = delta_angles(CVs[:, i], mins[i]) / deltas[i]
        else:
            normed[:, i] = (CVs[:, i] - mins[i]) / deltas[i]
    return normed

def normalise_lin(cv_traj, intermediaries, start, delta, opposites):
    # Initial normalisation before first string is created,
    # specialised for normalising those values
    # Shall be removed?
    traj = np.array(cv_traj)
    targets = np.array(intermediaries)
    start = np.array(start.get_CVs())

    delta = np.array(delta)
    print(start)
    for target in targets:
        print(target)
    input("pause")
    print("\n\ndelta", delta,"start", start, sep="\n")
    for i in range(len(delta)):
        if not opposites[i]:
            traj[:, i] = (traj[:, i] - start[i]) / delta[i]
            targets[:, i] = (targets[:, i] - start[i]) / delta[i]
        else:
            traj[:, i] = delta_angles(traj[:, i], start[i]) / delta[i]
            targets[:, i] = delta_angles(targets[:, i], start[i]) / delta[i]
    print ("\n\ntraj_norm")
    for one_traj in traj:
        print(one_traj)
    print ("\ntarget_norm")
    for one_target in targets:
        print(one_target)
    print("delta", delta,"start", start, sep="\n")
    return traj, targets

def sp_dim(number):
    cols = (4/3*number) ** 0.5
    rows = int(3/4*cols + 0.5)
    cols = int(cols + 0.5)
    return (rows, cols)

def plot_iterations_2D(phie,
                       psie,
                       init_states,
                       iterations,
                       CV_index = (0, 1),
                       spdim = (2, 3),
                       select = 1,
                       savedir = None,
                       CV_names = ["phi", "psi"]):
    iters = len(iterations)
    plt.figure()
    plt.plot(phie,psie,"ro-", label = "Initial")
    #print(CV_index)
    #print("init_states: ", init_states)
    init_x = init_states[:, CV_index[0]]
    init_y = init_states[:, CV_index[1]]
    plt.plot(init_x, init_y, "x")
    for i,string in enumerate(iterations[::select]):
        string.plot_CVs_2D(plt,
                           CV_index = CV_index,
                           label = str(i*select+1))
    plt.xlabel(CV_names[0])
    plt.ylabel(CV_names[1])
    plt.title("Every " + str(select) + "th iteration, linear interpolation")
    plt.legend()
    if savedir != None:
        plt.savefig(savedir+"/" + \
                    "".join([name.capitalize() for name in CV_names])\
                    + "every" + str(select) + "iterations.png",
                    dpi = 300)
        plt.close()


    subplots = spdim[0] * spdim[1]
    ppf = iters // subplots # ppf: plots per figure

    plt.figure()
    for i in range(subplots):
        plt.subplot(spdim[0], spdim[1], 1+i)
        plt.plot(phie,psie,"ro-", label = "Initial")
        plt.plot(init_x, init_y, "x")
        for j in range(ppf*i,ppf*(i+1)):
            string = iterations[j]
            string.plot_CVs_2D(plt,
                               CV_index = CV_index,
                               label = str(j+1))
        plt.title("Iterations " + str(ppf * i + 1) + " to " + str(ppf * (i + 1)))
        plt.xlabel(CV_names[0])
        plt.ylabel(CV_names[1])
        plt.legend()
    if savedir != None:
        plt.savefig(savedir+"/" + \
                    "".join([name.capitalize() for name in CV_names])\
                    + "iterations" + str(1) + "To" + str(subplots * ppf) +\
                    ".png",
                    dpi = 300)
        plt.close()
    if savedir == None:
        plt.show()

def plot_iter_splines_2D(phie,
                         psie,
                         init_states,
                         iterations,
                         CV_index = (0, 1),
                         select = 0,
                         ppf = 1,
                         savedir = None,
                         CV_names = ["phi", "psi"]):
    iters = len(iterations)
    if select is 0:
        select = range(iters)
    for i in select:
        if i % ppf == 0:
            plt.figure()
            plt.plot(phie,psie,"ro-", label = "Initial", color = 'k')
            plt.plot(init_states[:, CV_index[0]], init_states[:, CV_index[1]], "x")
        string = iterations[i]
        string.plot_spline_curve(plt,
                                 CV_index = CV_index,
                                 label = str(i+1))
        if (i - 1) % ppf == 0:
            plt.xlabel(CV_names[0])
            plt.ylabel(CV_names[0])
            plt.title("Iterations: " + str(i)+ "to " + str(i + ppf - 1))
            plt.legend()

            if savedir != None:
                plt.savefig(savedir+"/" + \
                            "".join([name.capitalize() for name in CV_names])\
                            + "Spline2DIter" + str(i) + "To" + \
                            str(i + ppf - 1) + ".png",
                            dpi = 300)
                plt.close()
    if savedir == None:
        plt.show()

def plot_sim_par_coords(panda_frame, savedir = None):
    fig = px.parallel_coordinates(panda_frame, color = "VISno")
    # fig = px.parallel_coordinates
    if savedir == None:
        fig.show()
    else:
        fig.write_html(savedir + "/parrallel_coordinates_plot.html")


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
        if key not in ["include"]:
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

def reparameterise(CVs, opposites, dih_no):
    mins, deltas = get_extremes(CVs, opposites)
    norm_CVs = normalise(CVs, mins, deltas, dih_no)
    splines = calc_splines(norm_CVs)
    new_norm_CVs = reparam_norm_cvs(splines, CVs.shape[0])
    new_CVs, drifts = denormalise(new_norm_CVs, CVs, mins, deltas, dih_no)
    return new_CVs, drifts, splines, mins, deltas

def save(data, file_name = "debug.pickle"):
    pickle.dump(data,open(file_name, "wb"))

def update_index_file(file_name = "index.ndx",
                      CVs = {},
                      pull_groups = {}):

    backup_file(file_name, copy = True)
    new_file = open(file_name, "a")
    new_file.write("\n\n[ dihedrals ]\n")
    for dihedral in CVs["dihedrals"]:
        new_file.write(dihedral.index_add() + "\n")
    for pull_group in pull_groups:
        new_file.write("\n" + pull_group.index_add() + "\n")
    new_file.close()

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
ITP_HEADERS["distance"] = "; ai   aj   type   index   type'      low     up1     up2     fac"
if __name__ == "__main__":
    mdp_create("simple.mdp")