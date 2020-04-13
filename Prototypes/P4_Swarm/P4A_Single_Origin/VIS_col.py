# VIS collection classes for VISs for the same system
# By: Marko Petrovic
# Date: 2020-04-13

import numpy as np
import VIS

class VIS_string:
    """ One string of VISs representing one iteration of the string method. """
    pass

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
        return(cvs)

class VIS_collection:
    """ A collection of several iterations. """
    pass
