import numpy as np
import os
from falass import dataformat


class Files(object):
    def __init__(self, pdbfile, lgtfile=None, datfile=None, resolution=5., ierror=5., flip=False):
        """
        Definition of the files to be used in the analysis; such as the .pdb file, .lgt file (is one exists), .dat file
        (is one exists). Also defined at this time are the resolution (this is ignored if the .dat file has 4 columns)
        and the percentage error in the intensity (this is also ignored if the .dat file has >2 columns), and whether
        the simulation cell should be flipped in the xy-plane.
        :param pdbfile: str
            path and name of the .pdb file to be analysed
        :param lgtfile: str (optional)
            path and name of the .lgt file (which contains the scattering lengths of each of the atom types in the
            pdbfile). If a lgtfile is not defined falass will help in the creation of one.
            Currently the .lgt file style that is supported is a 3 column space separated txt file where the columns are
            atom_type, real_scattering_length, and imaginary_scattering_length respectively.
        :param datfile: str (optional)
            path and name of the .dat file (from which the analysis q vectors are drawn and also for subsequent
            comparison between theory and experiment). If a datfile is not defined falass will allow the user to
            define a range of q vectors to calculate the reflectometry over.
            Currently the .dat file style that is supported is a 2, 3, and 4 column space separated txt files where the
            columns are q, i, di, and dq respectively.
        :param resolution: float (optional)
            percentage of the q vector for the width of the resolution Gaussian function to be used in the data
            smearing, if a 4 column datfile is given this is ignored.
        :param ierror: float (optional)
            percentage error of the intensity to be assumed, if a >2 column datfile is used this is ignored.
        :param flip: bool (optional)
            false if the system should be read as is, true is the simulation cell should be rotated through the
            xy-plane -- note that falass treats the first side that the neutron or X-ray interacts with as that at z=0.
        """
        self.pdbfile = pdbfile
        self.cell = []
        self.atoms = []
        self.number_of_timesteps = 0
        self.times = []
        self.lgtfile = lgtfile
        self.scat_lens = []
        self.real_scat_lens = {}
        self.imag_scat_lens = {}
        self.elements = []
        self.datfile = datfile
        self.expdata = []
        self.ierror = ierror
        self.resolution = resolution
        self.flip = flip

    def set_file(self, pdbfile=None, lgtfile=None, datfile=None):
        """
        Let the subsequent definition, or redefinition of the pdbfile, lgtfile or datfile.
        :param pdbfile: str
            path and name of the .pdb file to be analysed
        :param lgtfile: str (optional)
            path and name of the .lgt file (which contains the scattering lengths of each of the atom types in the
            pdbfile). If a lgtfile is not defined falass will help in the creation of one.
            Currently the .lgt file style that is supported is a 3 column space separated txt file where the columns are
            atom_type, real_scattering_length, and imaginary_scattering_length respectively.
        :param datfile: str (optional)
            path and name of the .dat file (from which the analysis q vectors are drawn and also for subsequent
            comparison between theory and experiment). If a datfile is not defined falass will allow the user to
            define a range of q vectors to calculate the reflectometry over.
            Currently the .dat file style that is supported is a 2, 3, and 4 column space separated txt files where the
            columns are q, i, di, and dq respectively.
        """
        if pdbfile:
            self.pdbfile = pdbfile
        if lgtfile:
            self.lgtfile = lgtfile
        if datfile:
            self.datfile = datfile

    def read_pdb(self):
        """
        Reads the .pdb file into memory. Currently the atoms must have the title 'ATOM', the timestep time needs to
        be the last text in the 'TITLE' line, and the cell dimensions are taken from the 'CRYST1' line, and assumed to
        be orthorhomic. Non-orthorhomic cells are not necessarily supported.
        """
        lines = line_count(self.pdbfile)
        print("Reading PDB file \n[ 0 % ]")
        file = open(self.pdbfile, 'r')
        string = 0
        atoms_each_timestep = []
        for i, line in enumerate(file):
            string_new = np.floor(i / lines * 100)
            if string_new > string + 9:
                string = string_new
                print("[{} {} % ]".format('#' * int(string / 10), int(string/10)*10))
            if line[0:6] == "ATOM  ":
                if self.flip:
                    a = np.sqrt(np.square(self.cell[self.number_of_timesteps - 1][2] - float(line[46:54])))
                    atoms_each_timestep.append(dataformat.atompositions(line[12:16].strip(), a))
                else:
                    atoms_each_timestep.append(dataformat.atompositions(line[12:16].strip(), float(line[46:54])))
            if "TITLE  " in line:
                if self.number_of_timesteps == 0:
                    self.number_of_timesteps += 1
                    self.times.append(float(line.split()[-1]))
                else:
                    self.atoms.append(atoms_each_timestep)
                    atoms_each_timestep = []
                    self.number_of_timesteps += 1
                    self.times.append(float(line.split()[-1]))
            if "CRYST1  " in line:
                self.cell.append([float(line[6:15]), float(line[15:24]), float(line[24:33])])
        print("[{} {} % ]".format('#' * int(100 / 10), int(100)))

    def read_lgt(self):
        """
        Parses the lgtfile. If no lgtfile is defined falass will help the user to build one by working through the
        atom types in the pdb file and requesting input of the real and imaginary scattering lengths. This will also
        occur if a atom type if found in the pdbfile but not in the given lgts file. falass will write the lgtfile 
        to disk if atom types do not feature in the given lgtfile or one is written from scratch. 
        """
        if self.lgtfile:
            lines = line_count(self.lgtfile)
            lgtfile_name = self.lgtfile
            print("Reading LGT file \n[ 0 % ]")
            file = open(self.lgtfile, 'r')
            string = 0
            for i, line in enumerate(file):
                string_new = np.floor(i / lines * 100)
                if string_new > string + 9:
                    string = string_new
                    print("[{} {} % ]".format('#' * int(string / 10), int(string/10)*10))
                line_list = line.split()
                duplicate = check_duplicates(self.scat_lens, line_list[0])
                if not duplicate:
                    self.scat_lens.append(dataformat.scatlens(line_list[0], float(line_list[1]), float(line_list[2])))
            for i in range(0, len(self.atoms)):
                for j in range(0, len(self.atoms[i])):
                    duplicate = check_duplicates(self.scat_lens, self.atoms[i][j].atom)
                    if not duplicate:
                        real_scat_len = input('The following atom type has no scattering length given '
                                              'in the lgt file {} \nPlease define a real scattering length for '
                                              'this atom type: '.format(self.atoms[i][j].atom))
                        imag_scat_len = input('\nPlease define a imaginary scattering length for '
                                              'this atom type: '.format(self.atoms[i][j].atom))
                        self.scat_lens.append(dataformat.scatlens(self.atoms[i][j].atom, float(real_scat_len),
                                                                  float(imag_scat_len)))
            print("[{} {} %]".format('#' * int(100 / 10), int(100)))
        else:
            print('There was no lgt file defined, falass will help you define one and save it for future use.')
            for i in range(0, len(self.atoms)):
                for j in range(0, len(self.atoms[i])):
                    duplicate = check_duplicates(self.scat_lens, self.atoms[i][j].atom)
                    if not duplicate:
                        real_scat_len = input('The following atom type has no scattering length given '
                                              'in the lgt file {} \nPlease define a real scattering length for '
                                              'this atom type: '.format(self.atoms[i][j].atom))
                        imag_scat_len = input('\nPlease define a imaginary scattering length for '
                                              'this atom type: '.format(self.atoms[i][j].atom))
                        self.scat_lens.append(dataformat.scatlens(self.atoms[i][j].atom, float(real_scat_len),
                                                                  float(imag_scat_len)))
            lgtfile_name = input("What should the lgt file be named? ")
            path, extension = os.path.splitext(lgtfile_name)
            if extension != '.lgt':
                lgtfile_name = path + '.lgt'
            print("[{} {} %]".format('#' * int(100 / 10), int(100)))
        lgtsf = open(lgtfile_name, 'w')
        for i in range(0, len(self.scat_lens)):
            lgtsf.write('{} {} {}\n'.format(self.scat_lens[i].atom, self.scat_lens[i].real * 1e5,
                                            self.scat_lens[i].imag * 1e5))

    def read_dat(self):
        """
        Parses the .dat file, supporting 2, 3, and 4 column files consisting of q, i, di, and dq with comments in lines
        where the first character is a '#'. If there is no .dat file the get_qs() function should be used to generate
        q vectors to allow for the calculation of the reflectometry profile.
        """
        if self.datfile:
            lines = line_count(self.datfile)
            print("Reading DAT file \n[ 0 % ]")
            file = open(self.datfile, 'r')
            string = 0
            for i, line in enumerate(file):
                string_new = np.floor(i / lines * 100)
                if string_new > string + 9:
                    string = string_new
                    print("[{} {} % ]".format('#' * int(string / 10), int(string/10)*10))
                if line[0] != '#':
                    line_list = line.split()
                    if len(line_list) == 2:
                        self.expdata.append(dataformat.datastruct(float(line_list[0]), float(line_list[1]),
                                                                  float(line_list[1]) * (self.ierror / 100),
                                                                  float(line_list[0]) * (self.resolution / 100)))
                    if len(line_list) == 3:
                        self.expdata.append(dataformat.datastruct(float(line_list[0]), float(line_list[1]),
                                                                  float(line_list[2]),
                                                                  float(line_list[0]) * (self.resolution / 100)))
                    if len(line_list) == 4:
                        self.expdata.append(
                            dataformat.datastruct(float(line_list[0]), float(line_list[1]), float(line_list[2]),
                                                  float(line_list[3])))
            print("[{} {} %]".format('#' * int(100 / 10), int(100)))
        else:
            print('No DAT file has been given, therefore no comparison will be conducted, please use the get_qs '
                  'function. Alternatively the DAT file can be added using the setFile function.')

    def get_qs(self, start=0.005, end=0.5, number=50):
        """
        If no datfile exists this function should be used to generate a linear-spaced range of q-vector for the
        calculation of the reflectometry over.
        :param start: float (optional)
            the first q-vector for which the reflectometry should be calculated
        :param end: float (optional)
            the last q-vector for which the reflectometry should be calculated
        :param number: int (optional)
            the number of q-vectors
        """
        q_values = np.linspace(start, end, number)
        for i in range(0, len(q_values)):
            self.expdata.append(dataformat.datastruct(q_values[i], None, None, q_values[i] * (self.resolution / 100)))



def check_duplicates(array, check):
    """
    Checks if an atom type has already been added to an array.
    :param array: array-type of scatlens
        the array to check
    :param check: str
        the atom type to try and find
    :return: bool
        true if the atom type is already present in the scatlen type array, false if not
    """
    for i in range(0, len(array)):
        if array[i].atom == check:
            return True
    return False


def line_count(filename):
    """
    quickly counts the number of lines in a file
    :param filename: str
        name of the file that the number of lines is desired for
    :return: int
        number of lines in the file
    """
    i = 0
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
