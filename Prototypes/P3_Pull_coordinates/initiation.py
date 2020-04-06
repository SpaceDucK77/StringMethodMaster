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
    return states
        

#def non_linear_interpolation():

def create_VIS(file, CVs):
    vis = VIS(file)
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

def create_string(start, end, intermediaries):
    n = len(intermediaries)
    
    

# Specific run with few variables and many constants.
def p3_run():
    # 1: Get Initial coordinates and CV angles
    start_angles = get_angle("start.gro", "test.ndx")
    stationary = get_cartesian("start.gro")
    
    # 2: Get Final angles
    end_angles = get_angle("end.gro", "test.ndx")
    
    # 3: Create start and end VISs
    startVIS = create_VIS("start.gro", start_angles)#, stationary)
    endVIS = create_VIS("end.gro", end_angles)#, stationary)
    
    # 4: Calculate CV VIS angles
    CVstates = linear_interpolation(start_angles, end_angles)
    
    
    # 5: Create string of VISs through coordinate pulling
    create_string(startVIS, endVIS, CVstates)


    return

if __name__ == "__main__":
    p3_run()
