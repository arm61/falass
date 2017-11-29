import numpy as np
from falass import readwrite, dataformat
import os


class Job():
    """The catch all.

    This class is used for setting up the falass job -- and is generally a catch all for inputs that do not fit into
    other parts.

    Parameters
    ----------
    files: falass.readwrite.Files
        A Files class item.
    layer_thickness: float
        The thickness of the layers that the simulation cell should be sliced into.
    cut_off_size: float
        The size of the simulation cell that should be ignored from the bottom -- this is to allow for the use
        of a vacuum gap at the bottom of the cell.
    """
    def __init__(self, files, layer_thickness, cut_off_size):
        self.files = files
        self.layer_thickness = layer_thickness
        self.cut_off_size = cut_off_size
        self.times = np.asarray(self.files.times)
        self.new_file = False

    def set_run(self, files=None, layer_thickness=None, cut_off_size=None):
        """Edit job inputs.

        This allows parts of the class to be assigned after the initial assignment or changed

        Parameters
        ----------
        files: falass.readwrite.Files
            A Files class item.
        layer_thickness: float
            The thickness of the layers that the simulation cell should be sliced into.
        cut_off_size: float
            The size of the simulation cell that should be ignored from the bottom -- this is to allow for the use
            of a vacuum gap at the bottom of the cell.
        """
        if files:
            self.files = files
            self.times = np.asarray(self.files.times)
        if layer_thickness:
            self.layer_thickness = layer_thickness
        if cut_off_size:
            self.cut_off_size = cut_off_size

    def set_lgts(self):
        """Assign scattering lengths.

        Assigned the scattering lengths from the lgtfile to the different atom types. If no lgtfile is defined falass
        will help the user to build one by working through the atom types in the pdb file and requesting input of the
        real and imaginary scattering lengths. This will also occur if a atom type if found in the pdbfile but not in
        the given lgts file. falass will write the lgtfile to disk if atom types do not feature in the given lgtfile or
        one is written from scratch.
        """
        if self.files.lgtfile:
            path, extension = os.path.splitext(self.files.lgtfile)
            lgtfile_name = path + extension
            for i in range(0, len(self.files.atoms)):
                for j in range(0, len(self.files.atoms[i])):
                    duplicate = readwrite.check_duplicates(self.files.scat_lens, self.files.atoms[i][j].atom)
                    if not duplicate:
                        self.new_file = True
                        real_scat_len = input('The following atom type has no scattering length given '
                                              'in the lgt file {} \nPlease define a real scattering length for '
                                              'this atom type: '.format(self.files.atoms[i][j].atom))
                        imag_scat_len = input('\nPlease define a imaginary scattering length for '
                                              'this atom type: '.format(self.files.atoms[i][j].atom))
                        self.files.scat_lens.append(dataformat.ScatLens(self.files.atoms[i][j].atom, float(real_scat_len),
                                                                  float(imag_scat_len)))
        else:
            self.new_file = True
            print('There was no lgt file defined, falass will help you define one and save it for future use.')
            for i in range(0, len(self.files.atoms)):
                for j in range(0, len(self.files.atoms[i])):
                    duplicate = check_duplicates(self.files.scat_lens, self.files.atoms[i][j].atom)
                    if not duplicate:
                        real_scat_len = input('The following atom type has no scattering length given '
                                              'in the lgt file {} \nPlease define a real scattering length for '
                                              'this atom type: '.format(self.files.atoms[i][j].atom))
                        imag_scat_len = input('\nPlease define a imaginary scattering length for '
                                              'this atom type: '.format(self.files.atoms[i][j].atom))
                        self.files.scat_lens.append(dataformat.ScatLens(self.files.atoms[i][j].atom, float(real_scat_len),
                                                                  float(imag_scat_len)))
            lgtfile_name = input("What should the lgt file be named? ")
            path, extension = os.path.splitext(lgtfile_name)
            if extension != '.lgt':
                lgtfile_name = path + '.lgt'
        if self.new_file:
            i = 0
            while os.path.isfile(lgtfile_name):
                i+=1
                lgtfile_name = path + str(i) + '.lgt'

            lgtsf = open(lgtfile_name, 'w')
            for i in range(0, len(self.files.scat_lens)):
                lgtsf.write('{} {} {}\n'.format(self.files.scat_lens[i].atom, self.files.scat_lens[i].real * 1e5,
                                                self.files.scat_lens[i].imag * 1e5))
            print('A new lgtfile has been written with the name {}'.format(lgtfile_name))

    def set_times(self, times=None):
        """Assign times to analyse.

        The assignment of the simulation timesteps that should be analysed. if none are given all will be analysed.

        Parameters
        ----------
        times: array_like float
            The timesteps that should be analysed, in the unit of time that present in the pdbfile.
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
    """Checks if item is in array.

    Checks if an item has already been added to an array.

    Parameters
    ----------
    array: array-type
        The array to check.
    check: str
        The item to try and find.

    Returns
    -------
    bool
        true if the item is already present in the scatlen type array, false if not.
    """
    for i in range(0, len(array)):
        if check == array[i]:
            return True
    return False
