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





    
# Specific run with few variables and many constants.
def p3_run(saves):
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
    return saves

def run_swarm(VIS):
    pass

def p4_run():
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

if __name__ == "__main__":
    p4_run()
