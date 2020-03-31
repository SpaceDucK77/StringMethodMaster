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


def linear_interpolation(start, end, parts = 10):
    pass

#def non_linear_interpolation():

def create_VIS():
    vis = VIS()

def create_string():
    pass

# Specific run with few variables and many constants.
def p3_run():
    # 1: Get Initial coordinates and CV angles
    start_angles = get_angle("start.gro", "test.ndx")
    stationary = get_cartesian("start.gro")
    
    # 2: Get Final angles
    end_angles = get_angle("end.gro", "test.ndx")
    # 3: Create start and end VISs
    # 4: Calculate CV VIS angles
    # 5: Create string of VISs through coordinate pulling



    return

if __name__ == "__main__":
    p3_run()
