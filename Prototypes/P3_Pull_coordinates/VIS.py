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

class VIS:
    # Virtual Initial State for Gromacs MD

    def __init__(self, conf, collection = None):
        self.mdp_settings = {}
        self.configuration_file = conf
        self.CV_keys = []
        self.CV_info = {}
        self.isReady = False
        self.collection = collection
        self.history = []
        self.isSolvated = False
        self.cv_index_file = "test.ndx"
        self.latest_path = ""

    def box(self):
        # To define MD Box
        pass


    def get_CVs(self, index_file = "test.ndx"):
        # Extracts Alanine dipeptide phi and chi angles
        return get_angle(self.configuration_file, index_file)

    def solvate(self, new_parameters={}):
        # Solvates box
        pass
        #self.isSolvated = True

    def ions(self, new_parameters={}):
        # Adding ions to box
        pass

    def EM(self, new_parameters={}):
        # Sets up for an EM MD run preparation
        self.mdp_settings["em"] = new_parameters
        self.single_run(name = "em",
                        template = "mdp_templates/minim.mdp.templ")

    def single_run(self, name, log = False, template = "", additions = None):
        # Preperes a run and runs it.
        mdp_create(file_name = name + ".mdp",
                   new_parameters = self.mdp_settings[name],
                   old_file = template)
        executable = "gmx"
        arguments = ["grompp"]
        input_files = {"-f": name + ".mdp",
                       "-c": self.configuration_file,
                       "-p": "topol.top"}
        output_files = {"-o" : name + ".tpr"}
        #print(name, "", executable, arguments, input_files, output_files,"","",sep="\n")
        if additions:
            if "in" in additions:
                append_dict(input_files, additions["in"])
            if "out" in additions:
                append_dict(output_files, additions["out"])
        #print(name, "", executable, arguments, input_files, output_files,"","",sep="\n")
        prep =  gmx.commandline_operation(executable = executable,
                                          arguments = arguments,
                                          input_files = input_files,
                                          output_files = output_files)
        prep.run()
        print("prep "+ name + ":\n", prep.output.erroroutput.result())

        mdrun = gmx.read_tpr(prep.output.file["-o"])
        md = gmx.mdrun(mdrun)
        md.run()
        path = md.output.trajectory.result()
        path = path[:path.rfind("/") + 1]
        self.configuration_file = path + "confout.gro"
        if log:
            self.history.append(self.configuration_file)
        self.latest_path = path
        shutil.move(name + ".mdp", path + name + ".mdp")
        #os.remove(name+".mdp")

    def nvt(self, new_parameters={}):
        # Sets up for a constant volume equilibration run preparation
        pass

    def npt(self, new_parameters={}):
        # Sets up for a constant preassure equilibration run preparation
        pass

    def steered(self, new_parameters={}):
        # Sets up for a steered run preparation
        # Command line: gmx grompp -f md_pull_test.mdp -c confout.gro -p topol.top -r confout.gro -n index.ndx -o pull.tpr -maxwarn 1
        extra_in = {"-r": self.configuration_file, "-n": "index.ndx"}
        additions = {"in": extra_in}
        self.mdp_settings["steered"] = new_parameters
        self.single_run(name = "steered",
                        template = "mdp_templates/md_steered_specific.mdp.templ",
                        additions = additions)

    def run(self):
        # Sets up for an MD run preparation
        pass

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
