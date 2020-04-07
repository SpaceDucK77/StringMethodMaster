import gmxapi as gmx
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
        file.write("{:25} = {}\n".format(key,value))


def read_mdp(file_name):
    file =  open(file_name)
    keys=[]
    parameters = {}
    for row in file.read().split("\n"):
        if "=" not in row.split(";")[0]:
            continue
        row = row.split("=")
        key = row[0].strip()
        value = row[1].strip()
        keys.append(key)
        parameters[key]=value
    file.close()
    return keys, parameters


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


if __name__ == "__main__":
    mdp_create("simple.mdp")
