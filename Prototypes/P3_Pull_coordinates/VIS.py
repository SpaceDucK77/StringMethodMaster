# Virtual Initial State class. Should probably have all
# data for an MD run.
# By: Marko Petrovic
# Date: 2020-03-30

import gmxapi as gmx
import subprocess as sp
import numpy as np
from md_tools import *

class VIS:

    def __init__(self):
        self.mdp_files = {}
        self.topology_file = ""
        self.CV_keys = []
        self.CV_constraints = {}
        self.isReady = False
