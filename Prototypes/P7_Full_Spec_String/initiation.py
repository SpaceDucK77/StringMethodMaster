# Prototype 4
# Preparing a string of Virtual Initial States(VIS) for MD-runs
# By: Marko Petrovic
# Date: 2020-04-12


import gmxapi as gmx
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
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
        print("Calculatinging start end CV values")
        start_angles = get_angle("start.gro", "test.ndx")
        saves["start_angles"] = start_angles
        save(saves)
    #stationary = get_cartesian("start.gro")
    # 2: Get Final angles
    if "end_angles" in saves:
        end_angles = saves["end_angles"]
    else:
        print("Calculating end end CV values")
        end_angles = get_angle("end.gro", "test.ndx")
        saves["end_angles"] = end_angles
        save(saves)

    # 3: Create start and end VISs
    if "startVIS" in saves and "endVIS" in saves:
        startVIS = saves["startVIS"]
        endVIS = saves["endVIS"]
    else:
        print("Creating end point VISs")
        startVIS = create_VIS("start.gro", start_angles)#, stationary)
        endVIS = create_VIS("end.gro", end_angles)#, stationary)
        saves["startVIS"] = startVIS
        saves["endVIS"] = endVIS
        save(saves)

    # 4: Calculate CV VIS angles
    if "CVstates" in saves:
        CVstates, delta = saves["CVstates"], saves["delta"]
    else:
        print("Calculating starting VIS CV values")
        CVstates, delta = linear_interpolation(start_angles, end_angles)
        saves["CVstates"], saves["delta"] = CVstates, delta
        save(saves)
    #input("old functioning?")


    # 5: Create string of VISs through coordinate pulling
    if "intermed" not in saves:
        print("Creating initial string of VISs")
        saves["intermed"] = create_string(startVIS, endVIS, CVstates, delta, saves)
        save(saves)
        
    return saves


def p4b_run(saves):
    saves = p3_run(saves)
    if "first_string" not in saves:
        print("Creating VIS string object")
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
    """
    print("before, intended:")
    print(np.array(saves["CVstates"]))
    print("before, actual:")
    print(saves["b_act"])
    print("\nafter:")
    print(end_loc)
    """
    return saves


def p5_run(saves):
    saves = p4b_run(saves)
    # print("Start angles:", saves["start_angles"], saves["startVIS"].get_CVs(), "\n", sep="\n")
    # print("End angles:", saves["end_angles"], saves["endVIS"].get_CVs(), "\n", sep="\n")
    saves["first_string"].prep_new_CVs()
    save(saves)
    if "second_string" not in saves:
        print("Creating a second iteration string")
        #x = saves["first_string"].prep_new_string()
        #print(x)
        #input()
        dist, second_string = saves["first_string"].prep_new_string()
        saves["second_string"] = second_string
        saves["dist"] = dist
    save(saves)
    return(saves)

def p7_run(saves = None):
    #reuse = input("use previous run (y/n)")=="y"
    reuse = True
    #reuse = True
    debug = True
    #debug = False
    if reuse:
        saves = load()
    else:
        saves = {}
    saves = p5_run(saves)
    if "iterations" not in saves:
        print("calculating new strings")
        iterations = [saves["second_string"]]
        distances = [saves["dist"]]
        dist = distances[-1]
        i = 0
        while dist > 5 and i < 30:
            iterations[-1].create_SO()
            iterations[-1].run_swarms()
            iterations[-1].prep_new_CVs()
            dist, new_string = iterations[-1].prep_new_string()
            iterations.append(new_string)
            distances.append(dist)
            i = i + 1
        saves["iterations"] = iterations
        saves["distances"] = distances
        save(saves)
    sa = saves["start_angles"]
    ea = saves["end_angles"]
    phie = np.array([sa[0], ea[0]])
    psie = np.array([sa[1], ea[1]])
    #plt.plot(phi,psi)
    
    # istates = np.array(saves["CVstates"])
    save(saves)
    #input("end")
    
if __name__ == "__main__":
    p7_run()
