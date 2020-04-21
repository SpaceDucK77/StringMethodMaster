# VIS collection classes for VISs for the same system
# By: Marko Petrovic
# Date: 2020-04-13

import numpy as np
import VIS

class VIS_string:
    """ One string of VISs representing one iteration of the string method. """
    def __init__(self, start, end, string, no_CVs):
        self.start = start
        self.end = end
        self.start_string = string
        self.cv_dim = no_CVs
        self.state = 0 # 0 = object created, 1 = string of origins created, 2 = swarms run
        self.SO = []

    def __str__(self):
        cvs = np.array([ swarm.averages() for swarm in self.SO])
        return "Swarm:\n" + str(cvs)

    def create_SO(self, redo = False):
        # Creates VIS_swarm objects (Swarm Origins (SOs))
        if self.state == 0 or redo:
            self.SO = [None] * len(self.start_string)
            for i in range(len(self.SO)):
                self.SO[i] = VIS_swarm(self.start_string[i])
        self.state = 1

    def run_swarms(self, redo = False):
        if self.state == 1 or (redo and self.state > 1):
            self.end_loc = np.zeros([len(self.SO),self.cv_dim])
            for i,swarm_or in enumerate(self.SO):
                swarm_or.run_swarm()
                self.end_loc[i,:] = swarm_or.averages()
            self.state = 2
        return self.end_loc




class VIS_swarm:
    """ One swarm of trajectories from a single origin. """
    def __init__(self, origin, trajectories = 20):
        self.origin = origin
        self.n_traj = trajectories
        self.trajs = [None] * trajectories
        for i in range(trajectories):
            self.trajs[i] = VIS.VIS.fresh_copy(origin, self)
        self.ran_swarm = False

    def run_swarm(self, parameters = None, run_override = False):
        if run_override or not self.ran_swarm:
            if parameters == None:
                parameters = {"nsteps": 50, "dt": 0.002}
            for vis in self.trajs:
                vis.swarm(parameters)
            self.ran_swarm = True

    def averages(self):
        n = len(self.trajs)
        cvs = []
        for vis in self.trajs:
            cvs.append(vis.get_CVs())
        cvs = np.array(cvs)
        return np.mean(cvs, axis = 0)

class VIS_collection:
    """ A collection of several iterations. """
    pass
