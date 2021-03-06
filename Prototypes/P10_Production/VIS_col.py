# VIS collection classes for VISs for the same system
# By: Marko Petrovic
# Date: 2020-04-13

import numpy as np
import matplotlib.pyplot as plt
import VIS
import shutil
from md_tools import *

# Dictionary associating parameter file key to attribute and type
P_FILE_DICT = {"name": ("name", str),
               "start_conf": ("start", str),
               "end_conf": ("end", str),
               "topology_file": ("topology", str),
               "struct_index": ("index", str),
               "solvated": ("solvated", int),
               "solvate": ("solvate", int),
               "max_iter": ("iterations", int),
               "max_conv": ("conv_limit", float),
               "beads": ("beads", int),
               "steered_run_time": ("steer_run", float),
               "swarm_size": ("swarm_size", int),
               "swarm_sim_steps": ("swarm_time", int)}


class Pull_group(frozenset):
# Pull group class

    # Constructor override
    def __new__(cls, atom_nos, pg_no = None):
        return super(Pull_group, cls).__new__(cls, atom_nos)

    # Constructor override
    def __init__(self, atom_nos, pg_no = None):
        self.COM_group_numbers = set()
        self.pg_no = pg_no
        self.name = ""
        self.__set_name()

    # adds COM groups to pull group
    def add_COM_group_number(self, number):
        self.COM_group_numbers.add(number)

    # name getter
    def get_name(self):
        return self.name

    # Generates a string for adding to index file
    def index_add(self):
        return "[ " + self.name + " ]\n" +\
            " ".join([str(i) for i in sorted(self)])

    # Name generator
    def __set_name(self):
        self.name = "grp_"
        if len(self) == 1:
            self.name += "a"+str(self.pg_no)
        else:
            self.name += "g"+str(self.pg_no)

# CV super class
class CV:

    #Constructor method
    def __init__(self, pull_groups):
        self.pull_groups = pull_groups
        self.name = None

    # Name setter
    def set_name(self, name):
        self.name = name

    # Name getter
    def get_name(self):
        return self.name


# Dihedral CV sublcass
class Dihedral(CV):

    # Constructor method
    def __init__(self, pull_groups):
        super(Dihedral, self).__init__(pull_groups)

    # Creates string for adding to index file
    def index_add(self):
        return " ".join([str(next(iter(i))) for i in self.pull_groups])

    # Returns which type of CV
    def geometry(self):
        return "dihedral"

    # creates string for .mdp file pull coordinate pul groups parameter
    def mdp_groups(self):
        return " ".join([str(pg) for pg in [self.pull_groups[cv].pg_no for cv in [0, 1, 1, 2, 2, 3]]])

# Distance CV subclass
class Distance(CV):

    # Constructor method
    def __init__(self, pull_groups):
        super(Distance, self).__init__(pull_groups)

    # Returns geometry type
    def geometry(self):
        return "distance"

    # creates string for .mdp file pull coordinate pul groups parameter
    def mdp_groups(self):
        return str(self.pull_groups[0].pg_no) + " " + str(self.pull_groups[1].pg_no)


class VIS_collection:
    """ A collection of several string iterations. """

    # Constructor method
    def __init__(self, parameter_file):
        self.parameter_file = parameter_file
        self.name = "String_method"
        self.start = None
        self.end = None
        self.topology = None
        self.conf = None
        self.index = None
        self.beads = None
        self.steer_run = 500
        self.swarm_size = 20
        self.swarm_time = 20
        self.strings = []
        self.state = 0
        self.curr_iter = 0
        self.progress = {} # saves replacement?
        self.pull_groups = {}
        self.CVs = {}
        self.end_pointCVs = []
        self.CV_vis = {}
        self.CV_2D = set()
        self.startVIS = None
        self.endVIS = None
        self.temp_dict = None

    # Generates initial string from end point VISs
    def create_base_string(self):
        if self.state == 2:
            self.CV_targets, self.delta, self.opposites = \
                linear_interpolation(self.end_pointCVs[0],
                                     self.end_pointCVs[1],
                                     parts = self.beads,
                                     no_dih = len(self.CVs["dihedrals"]))
            if self.temp_dict ==  None:
                self.temp_dict = {}
            #print(self.steer_run)
            #exit()
            self.beads_list = create_string(start = self.startVIS,
                                            end = self.endVIS,
                                            intermediaries = self.CV_targets,
                                            delta = self.delta,
                                            opposites = self.opposites,
                                            saves = self.temp_dict,
                                            collection = self,
                                            run_time = self.steer_run,
                                            name = self.name + "/Pre_String_0/")
            self.state = 2.5
            save(self)
        if self.state == 2.5:
            self.strings.append(VIS_string(start = self.startVIS,
                                           end = self.endVIS,
                                           string = self.beads_list,
                                           no_CVs = len(self.end_pointCVs[0]),
                                           name =  self.name + "/String_0/",
                                           trajs_per_swarm = self.swarm_size,
                                           swarm_time = self.swarm_time))
            self.state = 3
            save(self)

    # Parses CVs from .smg file content
    def parse_CVs(self, text):
        COM_group_no = {}
        current_pg_no = 1
        self.CVs={"dihedrals": [], "distances": []}
        divide = text.split("[ dihedrals ]")
        COM_group_text = divide[0]
        dihedrals, distances = divide[1].split("[ distances ]")
        for row in dihedrals.split("\n")[1:]:
            row = row.split("#")[0].strip()
            if len(row) == 0:
                continue
            current_pull_groups = []
            for atom in row.split():
                atom = frozenset([int(atom)])
                if atom not in self.pull_groups:
                    group = Pull_group(atom, current_pg_no)
                    self.pull_groups[group] = group
                    current_pg_no += 1
                current_pull_groups.append(self.pull_groups[atom])

            self.CVs["dihedrals"].append(Dihedral(current_pull_groups))

        curr_COM_group_no = 0
        for i,row in enumerate(COM_group_text.split("\n")[2:]):
            row = row.split("#")[0].strip()
            if len(row) == 0:
                continue
            row = [int(part.strip()) for part in row.split()]
            group = Pull_group(row, current_pg_no)
            if group not in self.pull_groups:
                self.pull_groups[group] = group
                current_pg_no += 1
            self.pull_groups[group].add_COM_group_number(curr_COM_group_no)
            curr_COM_group_no += 1
        for p_group in self.pull_groups:
            for g_no in p_group.COM_group_numbers:
                COM_group_no[g_no] = p_group
        for row in distances.split("\n")[1:]:
            row = row.split("#")[0].strip()
            if len(row) == 0:
                continue
            current_pull_groups = []
            for g_no in row.split():
                current_pull_groups.append(COM_group_no[int(g_no)])
            self.CVs["distances"].append(Distance(current_pull_groups))
        log("Pull groups")
        log("Dihedrals")
        for CV in self.CVs["dihedrals"]:
            log(str(CV.pull_groups))
        log("Distances")
        for CV in self.CVs["distances"]:
            log(str(CV.pull_groups))

    # Splits .smg file content into 3 parts. and calls their respective parser
    # Simulation parameters
    # CV parameters
    # Visualisation parameters
    def parse_parameters(self):
        if self.state > 0:
            return
        file = open(self.parameter_file)
        content = file.read().split('***')[0]
        self.parse_sim_parameters(content.split("---")[0])
        self.parse_CVs(content.split("---")[2])
        self.parse_visualisations(content.split("---")[4])
        try:
            os.mkdir(self.name)
        except:
            print("""Trying to run string method with same name as local directory.
Probably a previous string method run, this has been disallowed due to possible
naming coflicts\n""")
            log("""Trying to run string method with same name as local directory.
Probably a previous string method run, this has been disallowed due to possible
naming coflicts\n""")
            exit()
        self.state = 1
        save(self)

    # Parses simulation parameters from text
    def parse_sim_parameters(self, text):
        attributes = dictionarise(text)
        #print(attributes)
        for attribute in attributes:
            if attribute in P_FILE_DICT:
                setattr(self,
                        P_FILE_DICT[attribute][0],
                        P_FILE_DICT[attribute][1](attributes[attribute]))
            else:
                raise NameError("Parameter "+ attribute + " not recognised")
        #print(self.steer_run)
        #exit()

    # Parses visualisation parameters
    def parse_visualisations(self, text):
        for row in text.split("\n"):
            row = row.split("#")[0].strip()
            if len(row) == 0:
                continue
            plot = row[-1] == '*'
            row = row.split("=")
            name = row[0].strip()
            number = int(row[1].rstrip('*'))
            cv_list = self.CVs["dihedrals"]+self.CVs["distances"]
            cv_list[number].set_name(name)
            self.CV_vis[name] = number
            if plot:
                self.CV_2D.add(name)

    # Prepares end point VISs from configuration files
    def prepare_endpoints(self):
        if self.state == 1:
            backup_file(self.topology, True)
            update_topol_file(self.topology)
            make_index(self.start, self.index)
            update_index_file(self.index, self.CVs, self.pull_groups)
            dihedrals_exist = len(self.CVs["dihedrals"]) != 0
            startCVs = []
            endCVs = []
            if dihedrals_exist:
                startCVs = get_angle(self.start, self.index)
            for dist in self.CVs["distances"]:
                startCVs += get_distance (self.start, self.index, dist)
            if dihedrals_exist:
                endCVs =  get_angle(self.end, self.index)
            for dist in self.CVs["distances"]:
                endCVs += get_distance (self.end, self.index, dist)
            self.end_pointCVs = [startCVs, endCVs]
            self.startVIS = create_VIS(self.start,
                                       CV_definitions = self.CVs,
                                       topology_file = self.topology,
                                       index_file = self.index,
                                       pull_groups = self.pull_groups,
                                       solvated = self.solvated,
                                       path = self.name + "/startVIS/")
            self.endVIS = create_VIS(self.end,
                                     CV_definitions = self.CVs,
                                     topology_file = self.topology,
                                     index_file = self.index,
                                     pull_groups = self.pull_groups,
                                     solvated = self.solvated,
                                     path = self.name + "/endVIS/")
            if self.solvate and not self.solvated:
                self.startVIS.solvate()
                self.solvated = 1
                self.startVIS.ions()
                save(self)
            self.startVIS.EM(restrain_cv = True)
            if self.solvated:
                log("-------- NVT --------")
                self.startVIS.nvt()
                log("-------- NPT --------")
                self.startVIS.npt()
            self.state = 2
            save(self)


    # Iterates the string method
    def run_method(self):
        try:
            log_t = open("log.txt").read()
        except FileNotFoundError:
            log_t = ""
        if self.state == 3:
            dist = self.conv_limit + 1
            while dist > self.conv_limit and self.curr_iter < self.iterations:
                if str(self.curr_iter) + ": create_SO" not in log_t:
                    self.strings[-1].create_SO()
                    save(self)
                    log(str(self.curr_iter) + ": create_SO")
                print(str(self.curr_iter) + ": create_SO")
                if str(self.curr_iter) + ": run_swarms" not in log_t:
                    self.strings[-1].run_swarms()
                    save(self)
                    log(str(self.curr_iter) + ": run_swarms")
                print(str(self.curr_iter) + ": run_swarms")
                if str(self.curr_iter) + ": prep_new_CVs" not in log_t:
                    self.strings[-1].prep_new_CVs(self.opposites)
                    save(self)
                    log(str(self.curr_iter) + ": prep_new_CVs")
                print(str(self.curr_iter) + ": prep_new_CVs")
                if str(self.curr_iter) + ": prep_new_string" not in log_t:
                    dist, new_string = self.strings[-1].prep_new_string()
                    save(self)
                    #self.strings[-1].clear_string_storage_HD()
                    self.strings.append(new_string)
                    log(str(self.curr_iter) + ": prep_new_string")
                    save(self)
                print(str(self.curr_iter) + ": prep_new_string")
                self.curr_iter = len(self.strings) - 1
                save(self)
            self.state=4
            save(self)
        print("self.state:", self.state)

    # Generates visualisations
    def visualisations(self, lw):
        if self.state != 4:
            return
        if not os.path.isdir("plots"):
            os.mkdir("plots")
        sCV = self.startVIS.get_CVs()
        eCV = self.endVIS.get_CVs()
        CV_2D = list(self.CV_2D)
        CV_2D.sort(key = lambda x : self.CV_vis[x])
        for i,CV1 in enumerate(CV_2D):
            CV1index = self.CV_vis[CV1]
            for CV2 in CV_2D[i+1:]:
                CV2index = self.CV_vis[CV2]
                phie = np.array([sCV[CV1index], eCV[CV1index]])
                psie = np.array([sCV[CV2index], eCV[CV2index]])

                istates = np.array(self.CV_targets)
                CV_names = [CV1, CV2]
                print("VIS_coll CV_index: ", (CV1index, CV2index))
                print("VIS_coll CV_names: ", CV_names)

                plot_iterations_2D(phie,
                                   psie,
                                   istates,
                                   self.strings,
                                   CV_index = (CV1index, CV2index),
                                   spdim = sp_dim(len(self.strings) // 5 + 1),
                                   select = 3,
                                   savedir = "plots",
                                   CV_names = CV_names,
                                   lw = lw)
                plot_iter_splines_2D(phie,
                                     psie,
                                     istates,
                                     self.strings,
                                     CV_index = (CV1index, CV2index),
                                     select = list(range(len(self.strings))),
                                     ppf = 2,
                                     savedir = "plots",
                                     CV_names = CV_names,
                                     lw = lw)
                frame = create_panda(sCV,
                                     eCV,
                                     istates,
                                     self.strings,
                                     tuple([cv.get_name() for cv in self.CVs["dihedrals"]+self.CVs["distances"]]))
                plot_sim_par_coords(frame, "plots")





class VIS_string:
    """ One string of VISs representing one iteration of the string method. """

    # Constructor method
    def __init__(self, start, end, string, no_CVs, name, trajs_per_swarm = 20, swarm_time = 20):
        self.start = start
        self.end = end
        self.start_string = string
        self.cv_dim = no_CVs
        self.state = 0 # 0 = object created, 1 = string of origins created, 2 = swarms run
        self.SO = []
        self.start_CVs = None
        self.end_loc = None
        self.end_loc720 = None
        self.new_CVs = None
        self.target_drifts = None
        self.dih_no = len(self.start.get_CV_coll()["dihedral"])
        self.no_trajs = trajs_per_swarm
        self.swarm_time = swarm_time
        self.get_CVs()
        self.name = name
        try:
            os.mkdir(self.name)
        except FileExistsError:
            pass

    # String representation
    def __str__(self):
        cvs = np.array([ swarm.averages(origin_ok = True) for swarm in self.SO])
        return "Swarm:\n" + str(cvs)

    # Frees up disk space
    def clear_string_storage_HD(self):
        raise VIS.DeprecatedError("VIS_string.clear_string_storage_HD is no longer to be used.")
    """
        for swarm in self.SO:
            swarm.clear_swarm_storage_HD()"""

    # Creates VIS_swarm objects (Swarm Origins (SOs))
    def create_SO(self, redo = False):
        if self.state == 0 or redo:
            print("Creating swarms from VISs")
            self.SO = [None] * len(self.start_string)
            for i in range(len(self.SO)):
                self.SO[i] = VIS_swarm(self.start_string[i],
                                       trajectories = self.no_trajs,
                                       swarm_steps = self.swarm_time,
                                       name = self.name + "bead_" +
                                       str(i) + "/")
            self.state = 1

    # returns numpy matrix With the strings CV values
    def get_CVs(self):
        if self.new_CVs is not None:
            return self.new_CVs
        elif self.start_CVs is not None:
            return self.start_CVs
        start_CVs = [vis.get_CVs() for vis in self.start_string]
        self.start_CVs = np.append(np.append([self.start.get_CVs()],
                                             start_CVs,
                                             axis = 0),
                                   [self.end.get_CVs()],
                                   axis =0 )
        return self.start_CVs

    # plots itself on the given plot window
    def plot_CVs_2D(self, plotwindow, CV_index = (0,1), label = "A string", lw = 1):
        #self.get_CV_start()
        #print("VIS_string plot_CVs_2D, self.start_CVs", self.start_CVs.shape)
        p = plotwindow.plot(self.start_CVs[1:-1, CV_index[0]],
                            self.start_CVs[1:-1, CV_index[1]],"x")
        plotwindow.plot(self.start_CVs[:, CV_index[0]],
                        self.start_CVs[:, CV_index[1]],
                        p[0].get_color(),
                        label = label,
                        linewidth = lw)

    # Plots its full iteration to the given plot window
    def plot_spline_curve(self, plotwindow, CV_index = (0, 1), label = "A string", lw =1):
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
                            label = "spline: " + label,
                            linewidth = lw)
        plotwindow.plot(self.drift_CVs[1:-1, xdim],
                        self.drift_CVs[1:-1, ydim], "o",
                        color = p[0].get_color(),
                        label = "drift_CVs: " + label,
                        linewidth = lw)
        plotwindow.plot(self.new_CVs[1:-1, xdim],
                        self.new_CVs[1:-1, ydim], "v",
                        color = p[0].get_color(),
                        label = "new_CVs: " + label,
                        linewidth = lw)
        plotwindow.plot(self.start_CVs[:, xdim],
                        self.start_CVs[:, ydim],
                        label = "start_CVs: " + label,
                        linewidth = lw)

    # Calculates targets for the next string
    def prep_new_CVs(self, opposites, redo = False):
        if self.state == 2 or (redo and self.state > 2):
            print("Calculating CV values for new string")
            self.drift_CVs =  np.append(np.append([self.start.get_CVs()],
                                            self.end_loc,
                                            axis = 0),
                                  [self.end.get_CVs()],
                                  axis =0 )
            self.new_CVs, self.target_drifts, *self.spline_data =\
                    reparameterise(self.drift_CVs,
                                   opposites,
                                   len(self.start.get_CV_coll()["dihedral"]))
            self.state = 3

    # Creates a new string from its own origin and new targets
    def prep_new_string(self):
        if self.state == 3:
            old_CVs = [self.start.get_CVs()]
            new_string = [VIS.VIS.fresh_copy(old_state, new_name =
                                             self.name + "next_string_prep_" +
                                             str(i) + "/")
                          for i,old_state in enumerate(self.start_string)]
            nsteps = 1000
            sim_time = 0.002 * nsteps
            for i,fresh_copy in enumerate(new_string):
                parameters ={}
                parameters["nsteps"] = "{:10}; {}".format(nsteps, str(nsteps * 0.002) +" ps")
                old = self.start_string[i].get_CVs()
                old_CVs.append(old)
                for j in range(self.cv_dim):
                    if j < self.dih_no:
                        delta = delta_angle(self.new_CVs[i+1,j], old[j])
                    else:
                        delta = self.new_CVs[i+1,j] - old[j]
                    parameters["pull_coord" + str(j + 1) + "_rate"] = delta / sim_time
                    parameters["pull_coord" + str(j + 1) + "_k"] = 2000
                rem_path = fresh_copy.steered(parameters)
                if rem_path != None:
                    shutil.rmtree(rem_path)
            old_CVs.append(self.end.get_CVs())
            old_CVs = np.array(old_CVs)
            distance = np.linalg.norm(self.new_CVs - old_CVs)
        new_name = self.name.rsplit("_", 1)
        new_name = new_name[0]+"_" + str(int(new_name[1][:-1])+1)+"/"
        return distance, VIS_string(self.start,
                                    self.end,
                                    new_string,
                                    self.cv_dim,
                                    name = new_name,
                                    trajs_per_swarm = self.no_trajs)

    # Runs MD simulations for swarms of trajexctories for each bead
    def run_swarms(self, redo = False):
        if self.state == 1 or (redo and self.state > 1):
            print("Running swarms")
            self.end_loc = np.zeros([len(self.SO),self.cv_dim])
            self.end_loc720 = np.zeros([len(self.SO),self.cv_dim])
            for i,swarm_or in enumerate(self.SO):
                swarm_or.run_swarm()
                self.end_loc720[i,:], self.end_loc[i,:] = swarm_or.averages()
            self.state = 2
        return self.end_loc



class VIS_swarm:
    """ One swarm of trajectories from a single origin. """

    # Constructor method
    def __init__(self, origin, trajectories = 20, swarm_steps = 20, name = None):
        self.name = name
        try:
            os.mkdir(self.name)
        except FileExistsError:
            pass
        self.origin = VIS.VIS.fresh_copy(origin, self, new_name = name+"origin/")
        self.origin.EM(restrain_cv = True)
        self.origin.nvt()
        self.origin.npt()
        self.n_traj = trajectories
        self.steps = swarm_steps
        self.trajs = [None] * trajectories
        for i in range(trajectories):
            self.trajs[i] = VIS.VIS.fresh_copy(self.origin,
                                               new_name = self.name +
                                               "trajectory_"
                                               + str(i) + "/")
        self.ran_swarm = False
        self.origin_CVs = None
        self.CVs = None
        self.drift = None
        self.av_drift = None

    # Calculates swarms average drift
    def averages(self, origin_ok = False):
        if not self.ran_swarm and self.CVs is None and not origin_ok:
            raise ValueError("Trying to get average of non existing swarm")
        if self.ran_swarm == False and origin_ok:
            return np.array(self.origin.get_CVs())
        if self.CVs is None:
            n = len(self.trajs)
            CVs = []
            for vis in self.trajs:
                CVs.append(vis.get_CVs())
            self.CVs = np.array(CVs)
            self.drift = np.zeros(self.CVs.shape)
            self.origin_CVs = self.origin.get_CVs()
            CV_types = self.origin.get_CV_coll()
            CV_types = [len(CV_types[key]) for key in CV_types]
            col = 0
            for col in range(self.CVs.shape[1]):
                if col < CV_types[0]:
                    self.drift[:, col] = delta_angles(self.CVs[:, col],
                                                      self.origin_CVs[col])
                else:
                    self.drift[:, col] = self.CVs[:, col] - self.origin_CVs[col]
            self.av_drift = np.mean(self.drift, axis = 0)
        targets = np.zeros(self.av_drift.shape)
        targets720 = np.zeros(self.av_drift.shape)
        for i in range(self.av_drift.shape[0]):
            if i < CV_types[0]:
                targets[i] = delta_angle(self.av_drift[i], -self.origin_CVs[i])
            else:
                targets[i] = self.av_drift[i] + self.origin_CVs[i]
            targets720[i] = self.av_drift[i] + self.origin_CVs[i]
        return targets720, targets

    # Clears HD space
    def clear_swarm_storage_HD(self):
        raise VIS.DeprecatedError("VIS_swarm.clear_swarm_storage_HD is no longer to be used.")
    """
        for i,vis in enumerate(self.trajs):
            vis.delete_runs()
            #self.trajs[i] = VIS.VIS.fresh_copy(self.origin, self)
        self.ran_swarm = False"""

    # Runs all trajectories for this bead.
    def run_swarm(self, parameters = None, run_override = False):
        if run_override or not self.ran_swarm:
            if parameters == None:
                parameters = {"nsteps": self.steps, "dt": 0.002}
            for vis in self.trajs:
                vis.swarm(parameters)
            self.ran_swarm = True
        self.CVs = None


""" Do not use
if __name__ == "__main__":
    a = load()
    if a =={}:
        shutil.copyfile("topol_orig.top", "topol.top")
        a = VIS_collection("parameters.smg")
    a.parse_parameters()
    a.prepare_endpoints()
    a.create_base_string()
    a.run_method()
    a.visualisations()
"""
