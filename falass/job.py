import numpy as np


class Job(object):
    def __init__(self, files, layer_thickness, cut_off_size):
        """
        This class is used for setting up the falass job -- and is generally a catch all for inputs that do not fit into
        other parts.
        :param files: falass.readwrite.Files
            a Files class item
        :param layer_thickness: float
            the thickness of the layers that the simulation cell should be sliced into
        :param cut_off_size: float
            the size of the simulation cell that should be ignored from the bottom -- this is to allow for the use
            of a vacuum gap at the bottom of the cell
        """
        self.files = files
        self.layer_thickness = layer_thickness
        self.cut_off_size = cut_off_size
        self.times = np.asarray(self.files.times)

    def set_run(self, files=None, layer_thickness=None, cut_off_size=None):
        """
        this allows parts of the class to be assigned after the initial assignment or changed
        :param files: falass.readwrite.Files
            a Files class item
        :param layer_thickness: float
            the thickness of the layers that the simulation cell should be sliced into
        :param cut_off_size: float
            the size of the simulation cell that should be ignored from the bottom -- this is to allow for the use
            of a vacuum gap at the bottom of the cell
        """
        if files:
            self.files = files
            self.times = np.asarray(self.files.times)
        if layer_thickness:
            self.layer_thickness = layer_thickness
        if cut_off_size:
            self.cut_off_size = cut_off_size

    def set_times(self, times=None):
        """
        the assignment of the simulation timesteps that should be analysed. if none are given all will be anlaysed
        :param times: array_like float
            the timesteps that should be analysed, in the unit of time that present in the pdbfile
        """
        if times:
            self.times = np.arange(float(times[0]), float(times[1]) + float(times[2]), float(times[2]))
        else:
            first_times = float(input(
                "Please define the first timestep to be analysed, the first in the pdb file was {} ps: ".format(
                    self.files.times[0])))
            while check_array(self.files.times, first_times) is not True:
                first_times = float(input("TIMESTEP NOT FOUND. Please define the first timestep to be analysed, "
                                    "the first in the pdb file was {} ps: ".format(self.files.times[0])))
            last_times = float(input(
                "Please define the last timestep to be analysed, the last in the pdb file was {} ps: ".format(
                    self.files.times[-1])))
            while check_array(self.files.times, last_times) is not True:
                last_times = float(input("TIMESTEP NOT FOUND. Please define the last timestep to be analysed, "
                                   "the last in the pdb file was {} ps: ".format(self.files.times[-1])))
            interval_times = float(input("Please define time interval for analysis, the smallest interval in the pdb "
                                   "file was {} ps: ".format(self.files.times[1] - self.files.times[0])))
            while interval_times > last_times:
                interval_times = float(input("NOT A VALID INTERVAL. Please define time interval for analysis, "
                                       "the smallest interval in the pdb file was {} ps: ".format(self.files.times[1] -
                                                                                                  self.files.times[0])))
            self.times = np.arange(first_times, last_times + interval_times, interval_times)


def check_array(array, check):
    """
    Checks if an item has already been added to an array.
    :param array: array-type`
        the array to check
    :param check: str
        the item to try and find
    :return: bool
        true if the item is already present in the scatlen type array, false if not
    """
    for i in range(0, len(array)):
        if check == array[i]:
            return True
    return False
