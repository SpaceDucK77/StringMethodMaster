# Virtual Initial State class. Should probably have all
# data for an MD run.
# By: Marko Petrovic
# Date: 2020-03-30

import gmxapi as gmx
import subprocess as sp
import numpy as np
import os
import shutil
from md_tools import *

class DeprecatedError(Exception):
    pass

class VIS:
    # Virtual Initial State for Gromacs MD

    def __init__(self,
                 conf,
                 collection = None,
                 CV_definitions = None,
                 topology_file = "topol.top",
                 index_file = "index.ndx",
                 pull_groups = None,
                 solvated = 0,
                 restrain_file = "cv_restraints.itp",
                 generate_CV_mdp_dict = True):
        self.mdp_settings = {} # Specific .mdp settings for this system.
        self.configuration_file = conf # .gro file_name with configuration
        # self.start_conf = conf # for reference

        """ The following two aren't meant to be edited. """
        self.CV_keys = [] # names of CVs
        self.CV_info = {} # type: dihedral/distance,...

        self.CV_def = CV_definitions
        self.pull_groups = pull_groups

        self.is_ready = False # if the VIS is ready for a standard run
        self.collection = collection # the swarm, string, or full set this VIS is associated with.
        self.restrain_file = restrain_file
        self.history = [] #directories for all runs
        self.is_solvated = solvated # 1: solvated, 0: vacuum
        # self.cv_index_file = index_file # CV index file name, to be replaced
        self.index_file = index_file
        self.topology_file = topology_file
        self.latest_path = "" # Directory path of the latest run.
        
        self.CV_mdp_list = None
        self.CV_mdp_dict = None
        if generate_CV_mdp_dict == True:
            self.generate_CV_mdp_dict()

    def box(self):
        # Creates MD-simulation box
        file_name = self.configuration_file
        box = gmx.commandline_operation(executable ="gmx",
                                        arguments = [ "editconf","-c", "-d", "1.0", "-bt", "cubic"],
                                        input_files = {"-f": file_name},
                                        output_files= {"-o": "box" + file_name})
        box.run()
        shutil.move("box" + file_name, file_name)
        # os.remove("#" + file_name  + "#")


    def create_restrain_file(self):
        raise DeprecatedError("create_restrain_file no longer used")
    """
        cv_coll = self.get_CV_coll(self.index_file)
        t_file = open(self.restrain_file,"w")
        #in_file = open(self.index_file)
        cv_index = {}
        for collection in cv_coll:
            t_file.write("[ " + collection + "_restraints ]\n")
            t_file.write(ITP_HEADERS[cv_type] + "\n")
            for cv in cv_coll[collection]:
                ITP_WRITE[collection](cv.)
        for row in in_file.read().split("\n"):
            row = row.strip()
            if row == "":
                continue
            if row[0] == "[":
                cv_type = row.strip("[]").strip()[:-1]
                t_file.write("[ " + cv_type + "_restraints ]\n")
                t_file.write(ITP_HEADERS[cv_type] + "\n")
                cv_index[cv_type] = 0
            else:
                ITP_WRITE[cv_type]([int(s) for s in row.strip().split()],
                                   cv_coll[cv_type][cv_index[cv_type]],
                                   t_file)
                cv_index[cv_type] += 1
        in_file.close()
        t_file.close()"""

    def delete_runs(self):
        for path in self.history:
            log("clearing: " + path)
            for file_name in os.listdir(path):
                if ".xtc" in file_name or ".gro" in file_name:
                    continue
                file_name = path + file_name
                if os.path.isfile(file_name):
                    os.remove(file_name)
            

    def EM(self, new_parameters = None, new_keys = None, restrain_cv = False):
        # Sets up for an EM MD run preparation
        if new_parameters == None:
            new_parameters = {}
            new_keys = []
        if new_keys == None:
            new_keys = list(new_parameters.keys())
        sim_name = "em"
        additions = {}
        if restrain_cv:
            sim_name += "_r"
            newer_p = self.CV_mdp_dict.copy()
            append_dict(newer_p, new_parameters,
                        new_keys, self.CV_mdp_list)
            new_parameters = newer_p
            extra_in = {"-r": self.configuration_file, "-n": "index.ndx"}
            additions = {"in": extra_in}
 
        self.mdp_settings[sim_name] = new_parameters
        self.mdp_settings[sim_name + "_keys"] = new_keys
        prep = self.single_prep(name = sim_name,
                                template = "mdp_templates/minim.mdp.templ",
                                restrained = restrain_cv,
                                additions = additions)
        self.single_run(prep, sim_name)

    def fresh_copy(original, collection = None):
        copy = VIS(conf = original.configuration_file,
                   collection =  original.collection,
                   CV_definitions = original.CV_def,
                   topology_file = original.topology_file,
                   index_file = original.index_file,
                   pull_groups = original.pull_groups,
                   solvated =  original.is_solvated,
                   generate_CV_mdp_dict = False)
        #copy.mdp_settings = original.mdp_settings
        copy.CV_keys = original.CV_keys
        copy.CV_info = original.CV_info
        copy.is_ready = original.is_ready
        copy.latest_path = original.latest_path
        copy.CV_mdp_list = original.CV_mdp_list
        copy.CV_mdp_dict = original.CV_mdp_dict
        return copy

    def generate_CV_mdp_dict(self):
        self.CV_mdp_list = []
        slist = self.CV_mdp_list
        self.CV_mdp_dict = {}
        sdict = self.CV_mdp_dict
        slist.append("\n; Pull code")
        sdict["\n; Pull code"] = ""
        slist.append("pull")
        sdict["pull"] = "yes"
        slist.append("pull_ncoords")
        sdict["pull_ncoords"] = str(sum([len(self.CV_def[i]) for i in self.CV_def]))
        slist.append("pull_ngroups")
        sdict["pull_ngroups"] = str(len(self.pull_groups))
        for i, group in enumerate(self.pull_groups):
            key = "pull_group" + str(i+1) + "_name"
            slist.append(key)
            sdict[key] = group.get_name()
        i = 0
        pc = "pull_coord"
        keys = ["_type", "_geometry", "_dim",
                "_groups", "_start", "_rate", "_k"]
        for cv_type in self.CV_def:
            for cv in self.CV_def[cv_type]:
                st_dict = {"_type": "umbrella",
                           "_geometry": None,
                           "_dim": "Y Y Y",
                           "_groups": None,
                           "_start": "yes",
                           "_rate": "0",
                           "_k": "1000"}.copy()
                i += 1
                st_dict["_geometry"] = cv.geometry()
                st_dict["_groups"] = cv.mdp_groups()
                for stat in keys:
                    slist.append(pc + str(i) + stat)
                    sdict[pc + str(i) + stat] = st_dict[stat]

    def get_CV_coll(self, index_file = "index.ndx"):
        cv_coll = {}
        cv_coll["dihedral"] = get_angle(self.configuration_file, index_file)
        cv_coll["distance"] = []
        for cv in self.CV_def["distances"]:
            cv_coll["distance"].append(get_distance(self.configuration_file,
                                                    index_file, cv))
        return cv_coll

    # Replace with get_CV_coll later
    def get_CVs(self, index_file = "index.ndx"):
        # Extracts CVs as a list
        cv_coll = self.get_CV_coll()
        return cv_coll["dihedral"] + cv_coll["distance"]

    def ions(self, new_parameters = None, new_keys = None, ions = None):
        # Adding ions to box
        if new_parameters == None:
            new_parameters = {}
            new_keys = []
        if new_keys == None:
            new_keys = list(new_parameters.keys())
        #ion_top =  ".".join(self.configuration_file.split(".")[:-1]) + ".tpr"
        if ions is None or "pname" not in ions:
            pname = "NA" #get_pname()
        else:
            pname = ions["pname"]
        if ions is None or "nname" not in ions:
            nname = "CL" #get_nname()
        else:
            nname = ions["nname"]
        sim_name = "ions"
        self.mdp_settings[sim_name] = new_parameters
        self.mdp_settings[sim_name + "_keys"] = new_keys
        prep = self.single_prep(name = sim_name,
                                template = "mdp_templates/ions.mdp.templ")
        genion = gmx.commandline_operation(executable = "gmx",
                                           arguments =  ["genion",
                                                         "-pname", pname,
                                                         "-nname", nname,
                                                         "-neutral"],
                                           input_files = {"-s": prep.output.file["-o"],
                                                          "-p": "topol.top"},
                                           output_files = {"-o": "ion_"+self.configuration_file},
                                           stdin = "SOL")
        genion.run()
        self.configuration_file = "ion_"+self.configuration_file
        self.make_index()
        print("ionate:\n", genion.output.erroroutput.result())

    def make_index(self):
        make_index(self.configuration_file, self.index_file)
        update_index_file(self.index_file, self.CV_def, self. pull_groups)
                                                        
        
    def npt(self, new_parameters = None, new_keys = None):
        # Sets up for a constant preassure equilibration run preparation
        if self.is_solvated == 0:
            return
        if new_parameters == None:
            new_parameters = {}
            new_keys = []
        if new_keys == None:
            new_keys = list(new_parameters.keys())
        sim_name = "npt"
        additions = {}
        newer_p = self.CV_mdp_dict.copy()
        append_dict(newer_p, new_parameters,
                    new_keys, self.CV_mdp_list)
        new_parameters = newer_p
        self.mdp_settings[sim_name] = new_parameters
        self.mdp_settings[sim_name + "_keys"] = new_keys
        additions["in"] = {"-t": self.latest_path + "state.cpt",
                           "-r": self.configuration_file,
                           "-n": self.index_file}
        prep = self.single_prep(name = sim_name,
                                template = "mdp_templates/npt.mdp.templ",
                                additions = additions)
        self.single_run(prep, sim_name)

    def nvt(self, new_parameters = None, new_keys = None):
        # Sets up for a constant volume equilibration run preparation
        if self.is_solvated == 0:
            return        
        if new_parameters == None:
            new_parameters = {}
            new_keys = []
        if new_keys == None:
            new_keys = list(new_parameters.keys())
        sim_name = "nvt"
        additions = {}
        newer_p = self.CV_mdp_dict.copy()
        append_dict(newer_p, new_parameters,
                    new_keys, self.CV_mdp_list)
        new_parameters = newer_p
        
        self.mdp_settings[sim_name] = new_parameters
        self.mdp_settings[sim_name + "_keys"] = new_keys
        additions["in"] = {"-r": self.configuration_file,
                           "-n": self.index_file}
        prep = self.single_prep(name = sim_name,
                                template = "mdp_templates/nvt.mdp.templ",
                                additions = additions)
        self.single_run(prep, sim_name)
        
    def run(self):
        # Sets up for a standard MD run preparation
        pass

    def single_prep(self,
                    name,
                    #log = False,
                    template = "",
                    additions = None,
                    restrained = False):
        # Preperes a run and runs it.
        mdp_create(file_name = name + ".mdp",
                   new_parameters = self.mdp_settings[name],
                   new_keys = self.mdp_settings[name + "_keys"],
                   old_file = template)
        executable = "gmx"
        arguments = ["grompp"]
        input_files = {"-f": name + ".mdp",
                       "-c": self.configuration_file,
                       "-p": "topol.top"}
        output_files = {"-o" : name + ".tpr"}
        print("input_files:\n", input_files, "\noutput_files:\n", output_files)
        if additions:
            if "in" in additions:
                append_dict(input_files, additions["in"])
            if "out" in additions:
                append_dict(output_files, additions["out"])

        prep =  gmx.commandline_operation(executable = executable,
                                          arguments = arguments,
                                          input_files = input_files,
                                          output_files = output_files)
        prep.run()
        print("prep "+ name + ":\n", prep.output.erroroutput.result())
        return prep

    def single_run(self, prep, name, log = True):
        mdrun = gmx.read_tpr(prep.output.file["-o"])
        md = gmx.mdrun(mdrun)
        md.run()
        path = md.output.trajectory.result()
        path = path[:path.rfind("/") + 1]
        self.configuration_file = path + "confout.gro"
        if log:
            self.history.append(path)
        self.latest_path = path
        shutil.move(name + ".mdp", path + name + ".mdp")
        shutil.move(name + ".tpr", path + name + ".tpr")
        #os.remove(name+".mdp")

    def solvate(self, new_parameters = None):
        # Solvates box
        if 0 == self.is_solvated:
            file_name = self.configuration_file
            shutil.copyfile(self.topology_file, "unsol" + self.topology_file)
            solvate = gmx.commandline_operation(executable = "gmx",
                                                arguments = ["solvate"],
                                                input_files = {"-cp": file_name,
                                                               "-p": self.topology_file},
                                                output_files = {"-o": "solv_"+file_name})

            solvate.run()
            print("solvate:\n", solvate.output.erroroutput.result())
            self.configuration_file = "solv_"+file_name
            self.is_solvated = 1

    def split_traj(self):
        traj = gmx.commandline_operation(executable = "gmx",
                                  arguments = ["trjconv", "-sep"],
                                  input_files = {"-s": self.latest_path + "topol.tpr",
                                                 "-f": self.latest_path + "traj_comp.xtc"},
                                  output_files = {"-o": self.latest_path + "conf.gro"},
                                  stdin = "0")
        traj.run()
        print("split_traj:\n", traj.output.erroroutput.result())
        return self.latest_path

        # gmx trjconv -s pull.tpr -f pull.xtc -o conf.gro -sep

    def steered(self, new_parameters = None, new_keys = None):
        # Sets up for a steered run preparation
        # Command line: gmx grompp -f md_pull_test.mdp -c confout.gro -p topol.top -r confout.gro -n index.ndx -o pull.tpr -maxwarn 1
        if new_parameters == None:
            new_parameters = {}
            new_keys = []
        if new_keys == None:
            new_keys = list(new_parameters.keys())
        sim_name = "steered"
        extra_in = {"-r": self.configuration_file,
                    "-n": self.index_file}
        additions = {"in": extra_in}
        newer_p = self.CV_mdp_dict.copy()
        append_dict(newer_p, new_parameters,
                    new_keys, self.CV_mdp_list)
        new_parameters = newer_p
        self.mdp_settings[sim_name] = new_parameters
        self.mdp_settings[sim_name + "_keys"] = new_keys
        file_add = ""
        if self.is_solvated == 0:
            file_add = "_v"
        template = "mdp_templates/md_steered" + file_add + ".mdp.templ"
        prep = self.single_prep(name = "steered",
                                template = template,
                                additions = additions)
        self.single_run(prep, sim_name)

    def swarm(self, new_parameters = None, new_keys = None):
        # Shorter MD run to be part of swarm of trajectories
        if new_parameters == None:
            new_parameters = {}
            new_keys = []
        if new_keys == None:
            new_keys = list(new_parameters.keys())
        sim_name = "swarm"
        additions = {}
        self.mdp_settings[sim_name] = new_parameters
        self.mdp_settings[sim_name + "_keys"] = new_keys
        file_add = ""
        if self.is_solvated == 0:
            file_add = "_v"
        template = "mdp_templates/md" + file_add + ".mdp.templ"
        prep = self.single_prep(name = "swarm",
                                template = template,
                                additions =  additions)
        self.single_run(prep, sim_name)
