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
    # Virtual Initial State for Gromacs MD

    def __init__(self, conf, collection = None):
        self.mdp_settings = {}
        self.configuration_file = conf
        self.CV_keys = []
        self.CV_info = {}
        self.isReady = False
        self.collection = collection
        self.history = []

    def box(self):
        # To define MD Box
        pass


    def get_CVs(self):
        # Extracts Alanine dipeptide phi and chi angles
        return get_angle(self.configuration_file, "test.ndx")
    
    def solvate(self):
        # Solvates box
        pass

    def ions(self):
        # Adding ions to box
        pass

    def EM(self):
        # Sets up for an EM MD run preparation 
        self.single_run(name = "em",
                        template = "mdp_templates/minim.mdp.templ")

    def single_run(self, name, log = False, template = ""):
        # Preperes a run and runs it.
        mdp_create(fil_name = name + ".mdp",
                   new_parameters = self.mdp_settings[name],
                   old_file = template)
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
        # Sets up for a constant volume equilibration run preparation
        pass

    def npt(self):
        # Sets up for a constant preassure equilibration run preparation
        pass

    def steered(self):
        # Sets up for a steered run preparation
        pass

    def run(self):
        # Sets up for an MD run preparation
        pass
