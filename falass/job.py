import numpy as np

class Job(object):
    def __init__(self, files, layer_thickness, cut_off_size, flip = False):
        self.files = files
        self.layer_thickness = layer_thickness
        self.cut_off_size = cut_off_size
        self.flip = flip
        self.times = []

    def setRun(self, files = None, layer_thickness = None, cut_off_size = None, flip = None):
        if files:
            self.files = files
        if layer_thickness:
            self.layer_thickness = layer_thickness
        if cut_off_size:
            self.cut_off_size = cut_off_size
        if flip:
            self.flip = flip

    def setTimes(self, times = None):
        if times:
            self.times = np.arange(times[0], times[1], times[2])
        else:
            first_times = float(input(
                "Please define the first timestep to be analysed, the first in the pdb file was {} ps: ".format(
                    self.files.times[0])))
            while check_array(self.files.times, first_times) != True:
                first_times = float(input("TIMESTEP NOT FOUND. Please define the first timestep to be analysed, "
                                    "the first in the pdb file was {} ps: ".format(self.files.times[0])))
            last_times = float(input(
                "Please define the last timestep to be analysed, the last in the pdb file was {} ps: ".format(
                    self.files.times[-1])))
            while check_array(self.files.times, last_times) != True:
                last_times = float(input("TIMESTEP NOT FOUND. Please define the last timestep to be analysed, "
                                   "the last in the pdb file was {} ps: ".format(self.files.times[-1])))
            interval_times = float(input("Please define time interval for analysis, the smallest interval in the pdb "
                                   "file was {} ps: ".format(self.files.times[1] - self.files.times[0])))
            while interval_times > last_times:
                interval_times = float(input("NOT A VALID INTERVAL. Please define time interval for analysis, "
                                       "the smallest interval in the pdb file was {} ps: ".format(
                    self.files.times[1] - self.files.times[0])))
            self.times = np.arange(first_times, last_times, interval_times)


def check_array(array, check):
    for i in range(0, len(array)):
        if check == array[i]:
            return True
    return False