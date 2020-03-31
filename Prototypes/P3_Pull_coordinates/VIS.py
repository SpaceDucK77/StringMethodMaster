# Virtual Initial State class. Should probably have all
# data for an MD run.
# By: Marko Petrovic
# Date: 2020-03-30

import gmxapi as gmx
import subprocess as sp
import numpy as np
import os
from md_tools import *

class VIS:

    def __init__(self, conf, collection = None):
        self.mdp_settings = {}
        self.configuration_file = conf
        self.CV_keys = []
        self.CV_info = {}
        self.isReady = False
        self.collection = collection
        self.history = []

    def box(self):
        pass


    def get_CVs(self):
        return get_angle(self.configuration_file, "test.ndx")
    
    def solvate(self):
        pass

    def ions(self):
        pass

    def EM(self):
        keys, self.mdp_settings["em"] = read_mdp("mdp_templates/minim.mdp.templ")
        self.single_run("em")

    def single_run(self, name, log = False):
        mdp_create(name + ".mdp", self.mdp_settings[name])
        prep = gmx.commandline_operation(executable = "gmx", 
                                         arguments = ["grompp"],
                                         input_files = {"-f": name + ".mdp",
                                                        "-c": self.configuration_file,
                                                        "-p": "topol.top"},
                                         output_files = {"-o" : name + ".tpr"})
        prep.run()
        print("prep", prep.output.erroroutput.result())

        mdrun = gmx.read_tpr(prep.output.file["-o"])
        md = gmx.mdrun(mdrun)
        md.run()
        path = md.output.trajectory.result()
        self.structure_file = path[:path.rfind("/") + 1] + "confout.gro"
        if log:
            self.history.append(self.configuration_file)
        os.remove(name+".mdp")

    def nvt(self):
        pass

    def npt(self):
        pass

    def run(self):
        pass
