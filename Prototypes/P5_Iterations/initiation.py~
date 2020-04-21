# Prototype 4
# Preparing a string of Virtual Initial States(VIS) for MD-runs
# By: Marko Petrovic
# Date: 2020-04-12


import gmxapi as gmx
import subprocess as sp
import numpy as np
import os
import VIS
import VIS_col as vc
from md_tools import *

class Debug(Exception):
    pass

def all_trajs(saves):
    for i in range(9):
        print("\nSwarm {:<3}".format(i))
        for j in range(20):
            print(saves["first_string"].SO[i].trajs[j].get_CVs())

def debug_func(saves):
    comm = input("debug_func: >>>")
    while comm != "":
        try:
            print(eval(comm))
        except KeyError:
            print ("Variable missing")
        except:
            print("invalid command?")
        comm = input("debug_func: >>>")
        
# Specific run with few variables and many constants.
def p3_run(saves):
    # 0: Update topology file for cv_restraints
    update_topol_file()
    # 1: Get Initial coordinates and CV angles

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


    # 5: Create string of VISs through coordinate pulling
    if "intermed" not in saves:
        saves["intermed"] = create_string(startVIS, endVIS, CVstates, delta, saves)
        save(saves)
    if "intermedEM" not in saves:
        for vis in saves["intermed"]:
            vis.EM(restrain_cv = True)
        saves["intermedEM"] = True
        save(saves)
        
    return saves


def p4a_run():
    #reuse = input("use previous run (y/n)")=="y"
    reuse = False
    #reuse = True
    #debug = True
    debug = False
    if reuse:
        saves = load()
        p3_run(saves)
    else:
        saves = p3_run({})
    if "one" not in saves:
        saves["one"] = saves["intermed"][0]
        save(saves)
    if "first_swarm" not in saves or debug:
        saves["first_swarm"] = vc.VIS_swarm(saves["one"])
        save(saves)
    saves["first_swarm"].run_swarm()
    save(saves)
    aver = saves["first_swarm"].averages()
    aver = np.mean(aver, axis = 0)

def p4b_run():
    #reuse = input("use previous run (y/n)")=="y"
    reuse = True
    #reuse = True
    debug = True
    #debug = False
    if reuse:
        saves = load()
        p3_run(saves)
    else:
        saves = p3_run({})
    if "first_string" not in saves:
        saves["first_string"] = vc.VIS_string(start = saves["startVIS"],
                                              end = saves["endVIS"],
                                              string = saves["intermed"],
                                              no_CVs = 2)
        save(saves)
    saves["first_string"].create_SO()
    if "b_act" not in saves:
        saves["b_act"] = str(saves["first_string"])
    save(saves)
    end_loc = saves["first_string"].run_swarms()
    save(saves)
    print("before, intended:")
    print(np.array(saves["CVstates"]))
    print("before, actual:")
    print(saves["b_act"])
    print("\nafter:")
    print(end_loc)

    if debug:
        all_trajs(saves)
        #debug_func(saves)

if __name__ == "__main__":
    p4b_run()
