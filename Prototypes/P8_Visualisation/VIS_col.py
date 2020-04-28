# VIS collection classes for VISs for the same system
# By: Marko Petrovic
# Date: 2020-04-13

import numpy as np
import matplotlib.pyplot as plt
import VIS
from md_tools import *


class VIS_string:
    """ One string of VISs representing one iteration of the string method. """
    def __init__(self, start, end, string, no_CVs):
        self.start = start
        self.end = end
        self.start_string = string
        self.cv_dim = no_CVs
        self.state = 0 # 0 = object created, 1 = string of origins created, 2 = swarms run
        self.SO = []
        self.start_CVs = None
        self.end_loc = None
        self.new_CVs = None
        self.get_CVs()

    def __str__(self):
        cvs = np.array([ swarm.averages() for swarm in self.SO])
        return "Swarm:\n" + str(cvs)

    def create_SO(self, redo = False):
        # Creates VIS_swarm objects (Swarm Origins (SOs))
        if self.state == 0 or redo:
            print("Creating swarms from VISs")
            self.SO = [None] * len(self.start_string)
            for i in range(len(self.SO)):
                self.SO[i] = VIS_swarm(self.start_string[i])
            self.state = 1
        
            
    def get_CVs(self):
        #try:
        if self.new_CVs is not None:
            return self.new_CVs
        elif self.start_CVs is not None:
            return self.start_CVs
        #except AttributeError:
        #    self.new_CVs = None
        start_CVs = [vis.get_CVs() for vis in self.start_string]
        self.start_CVs = np.append(np.append([self.start.get_CVs()],
                                             start_CVs,
                                             axis = 0),
                                   [self.end.get_CVs()],
                                   axis =0 )
        return self.start_CVs

    ''' def get_CV_start(self):
        #try:
        if self.start_CVs is not None:
            return self.start_CVs
        #except AttributeError:
        #    self.new_CVs = None
        start_CVs = [vis.get_CVs() for vis in self.start_string]
        self.start_CVs = np.append(np.append([self.start.get_CVs()],
                                             start_CVs,
                                             axis = 0),
                                   [self.end.get_CVs()],
                                   axis =0 )
        return self.start_CVs'''

    def plot_CVs_2D(self, plotwindow, CV_index = (0,1), label = "A string"):
        #self.get_CV_start()
        p = plotwindow.plot(self.start_CVs[1:-1, CV_index[0]],
                            self.start_CVs[1:-1, CV_index[1]],"x")
        plotwindow.plot(self.start_CVs[:, CV_index[0]],
                        self.start_CVs[:, CV_index[1]],
                        p[0].get_color(),
                        label = label)

    def plot_spline_curve(self, plotwindow, CV_index = (0, 1), label = "A string"):
        if "spline_data" not in dir(self) and self.state > 1:
            state = self.state
            if state == 2:
                state = 3
            self.prep_new_CVs(redo = True)
            self.state = state
        elif self.state < 2:
            print("error label: " + label)
            return None
        splines = self.spline_data[0]
        mins = self.spline_data[1]
        deltas = self.spline_data[2]
        states = len(self.SO) + 1
        t = np.linspace(0, states, 200)
        xdim = CV_index[0]
        ydim = CV_index[1]
        x = splines[xdim](t) * deltas[xdim] + mins[xdim]
        y = splines[ydim](t) * deltas[ydim] + mins[ydim]
        p = plotwindow.plot(x,
                            y,
                            label = "spline: " + label)
        plotwindow.plot(self.drift_CVs[1:-1, xdim],
                        self.drift_CVs[1:-1, ydim], "o", 
                        color = p[0].get_color(),
                        label = "drift_CVs: " + label)
        plotwindow.plot(self.new_CVs[1:-1, xdim],
                        self.new_CVs[1:-1, ydim], "v",
                        color = p[0].get_color(),
                        label = "new_CVs: " + label)
        plotwindow.plot(self.start_CVs[:, xdim],
                        self.start_CVs[:, ydim],
                        label = "start_CVs: " + label)
    def prep_new_CVs(self, redo = False):
        if self.state == 2 or (redo and self.state > 2):
            print("Calculating CV values for new string")
            self.drift_CVs =  np.append(np.append([self.start.get_CVs()],
                                            self.end_loc,
                                            axis = 0),
                                  [self.end.get_CVs()],
                                  axis =0 )
            self.new_CVs, *self.spline_data = reparameterise(self.drift_CVs)
            self.state = 3

    def prep_new_string(self):
        if self.state == 3:
            old_CVs = [self.start.get_CVs()]
            new_string = [VIS.VIS.fresh_copy(old_state) for old_state in self.start_string]
            nsteps = 1000
            for i,fresh_copy in enumerate(new_string):
                parameters ={}
                parameters["nsteps"] = "{:10}; {}".format(nsteps, str(nsteps * 0.002) +" ps")
                old = self.start_string[i].get_CVs()
                old_CVs.append(old)
                for j in range(self.cv_dim):
                    delta = self.new_CVs[i+1,j] - old[j]
                    parameters["pull_coord" + str(j + 1) + "_rate"] = delta * 0.002 / nsteps
                    parameters["pull_coord" + str(j + 1) + "_k"] = 2000
                fresh_copy.steered(parameters)
            old_CVs.append(self.end.get_CVs())
            old_CVs = np.array(old_CVs)
            distance = np.linalg.norm(self.new_CVs - old_CVs)
        return distance, VIS_string(self.start, self.end, new_string, self.cv_dim)

    def run_swarms(self, redo = False):
        if self.state == 1 or (redo and self.state > 1):
            print("Running swarms")
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
                parameters = {"nsteps": 20, "dt": 0.002}
            for vis in self.trajs:
                vis.EM(restrain_cv = True)
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

    def __init__(self, start, end, topolgy, cv_def):
        self.start = start
        self.end = end
        self.topology = topology
        self.conf = cv_def
        self.strings = []
        self.state = 0

    def create_base_string():
        # self.extract_CV_end_point_values()
        # self.create_end_points()
        # create the string
        # self.state = 1
        pass

    def iterate_string():
        # self.current_string.create_SO()
        # self.current_string.run_swarms()
        
        pass

    def run_method():
        pass
