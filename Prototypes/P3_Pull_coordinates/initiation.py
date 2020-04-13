# Prototype 3
# Preparing a string of Virtual Initial States(VIS) for MD-runs
# By: Marko Petrovic
# Date: 2020-03-30

# Command example:
# gmx gangle -s end.gro -n test.ndx -g1 dihedral -group1 0 -oall test.xvg

# Remeber to normalise CVs before working pathing and path coordinates.
import gmxapi as gmx
import subprocess as sp
import numpy as np
import os
import pickle
from VIS import *
from md_tools import *



def linear_interpolation(start, end, parts = 10):
    states = [None]*(parts-1)
    delta = []
    for i in range(len(start)):
        delta.append(end[i]-start[i])
    for i in range(1, parts):
        state = []
        for j in range(len(start)):
            state.append(start[j]+delta[j]*i/parts)
        states[i-1]=state
    return states, delta


#def non_linear_interpolation():

def create_VIS(file, CV_values=None):
    vis = VIS(file, CV_values)
    vis.EM()
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

def create_string(start, end, intermediaries, delta, saves):
    n = len(intermediaries)
    if "path" in saves:
        cv_span = saves["cv_span"]
        path = saves["path"]
    else:
        cv_span={}
        for i in range(len(delta)):
            cv_span["pull_coord" + str(i+1) +"_rate"] = delta[i]/500
        start.steered(new_parameters = cv_span)
        path = start.split_traj()
        saves["cv_span"] = cv_span
        saves["path"] = path
        save(saves)
    if "cv_traj" in saves:
        cv_traj = saves["cv_traj"]
    else:
        cv_traj = get_angles(path = path)
        saves["cv_traj"] = cv_traj
        save(saves)
    n_traj, n_targ = normalise(cv_traj, intermediaries, start, delta)
    indexes =  find_matches(n_traj, n_targ)
    print(indexes)
    VIS_collection = []
    for i in indexes:
        VIS_collection.append(VIS(path+"conf"+str(i)+".gro"))
    return VIS_collection


def load():
    try:
        return pickle.load(open("debug.pickle","rb"))
    except FileNotFoundError:
        return {}

def save(data):
    pickle.dump(data,open("debug.pickle","wb"))

# Specific run with few variables and many constants.
def p3_run():
    # 1: Get Initial coordinates and CV angles
    debug = input("debug (y/n)")=="y"
    #debug = True
    if debug:
        saves = load()
    else:
        saves={}

    if "start_angles" in saves:
        start_angles = saves["start_angles"]
    else:
        start_angles = get_angle("start.gro", "test.ndx")
        saves["start_angles"] = start_angles
        save(saves)
    #stationary = get_cartesian("start.gro")

    # 2: Get Final angles
    if "end_angles" in saves:
        end_angles = saves["end_angles"]
    else:
        end_angles = get_angle("end.gro", "test.ndx")
        saves["end_angles"] = end_angles
        save(saves)

    # 3: Create start and end VISs
    if "startVIS" in saves and "endVIS" in saves:
        startVIS = saves["startVIS"]
        endVIS = saves["endVIS"]
    else:
        startVIS = create_VIS("start.gro", start_angles)#, stationary)
        endVIS = create_VIS("end.gro", end_angles)#, stationary)
        saves["startVIS"] = startVIS
        saves["endVIS"] = endVIS
        save(saves)

    # 4: Calculate CV VIS angles
    if "CVstates" in saves:
        CVstates, delta = saves["CVstates"], saves["delta"]
    else:
        CVstates, delta = linear_interpolation(start_angles, end_angles)
        saves["CVstates"], saves["delta"] = CVstates, delta
        save(saves)
    #input("old functioning?")

    #print(delta)
    # 5: Create string of VISs through coordinate pulling
    create_string(startVIS, endVIS, CVstates, delta, saves)


    return

if __name__ == "__main__":
    p3_run()
