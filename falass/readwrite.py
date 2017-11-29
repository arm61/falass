import numpy as np
from falass import dataformat


class Files:
    """File parsing.

    Definition of the files to be used in the analysis; such as the .pdb file, .lgt file (is one exists), .dat file
    (is one exists). Also defined at this time are the resolution (this is ignored if the .dat file has 4 columns)
    and the percentage error in the intensity (this is also ignored if the .dat file has >2 columns), and whether
    the simulation cell should be flipped in the xy-plane.

    Parameters
    ----------
    pdbfile: str
        Path and name of the .pdb file to be analysed.
    lgtfile: str, optional
        Path and name of the .lgt file (which contains the scattering lengths of each of the atom types in the
        pdbfile). If a lgtfile is not defined falass will help in the creation of one.
        Currently the .lgt file style that is supported is a 3 column space separated txt file where the columns are
        atom_type, real_scattering_length, and imaginary_scattering_length respectively.
    datfile: str, optional
        Path and name of the .dat file (from which the analysis q vectors are drawn and also for subsequent
        comparison between theory and experiment). If a datfile is not defined falass will allow the user to
        define a range of q vectors to calculate the reflectometry over.
        Currently the .dat file style that is supported is a 2, 3, and 4 column space separated txt files where the
        columns are q, i, di, and dq respectively.
    resolution: float, optional
        Percentage of the q vector for the width of the resolution Gaussian function to be used in the data
        smearing, if a 4 column datfile is given this is ignored.
    ierror: float, optional
        Percentage error of the intensity to be assumed, if a >2 column datfile is used this is ignored.
    flip: bool, optional
        False if the system should be read as is, true is the simulation cell should be rotated through the
        xy-plane -- note that falass treats the first side that the neutron or X-ray interacts with as that at z=0.
    """
    def __init__(self, pdbfile, lgtfile=None, datfile=None, resolution=5., ierror=5., flip=False):
        self.pdbfile = pdbfile
        self.cell = []
        self.atoms = []
        self.number_of_timesteps = 0
        self.times = []
        self.lgtfile = lgtfile
        self.scat_lens = []
        self.datfile = datfile
        self.expdata = []
        self.ierror = ierror
        self.resolution = resolution
        self.flip = flip
        return

    def set_file(self, pdbfile=None, lgtfile=None, datfile=None):
        """Edits files.

        Let the subsequent definition, or redefinition of the pdbfile, lgtfile or datfile.

        Parameters
        ----------
        pdbfile: str
            Path and name of the .pdb file to be analysed.
        lgtfile: str, optional
            Path and name of the .lgt file (which contains the scattering lengths of each of the atom types in the
            pdbfile). If a lgtfile is not defined falass will help in the creation of one.
            Currently the .lgt file style that is supported is a 3 column space separated txt file where the columns are
            atom_type, real_scattering_length, and imaginary_scattering_length respectively.
        datfile: str, optional
            Path and name of the .dat file (from which the analysis q vectors are drawn and also for subsequent
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
        return

    def read_pdb(self):
        """Parse .pdb.

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
                    atoms_each_timestep.append(dataformat.AtomPositions(line[12:16].strip(), a))
                else:
                    atoms_each_timestep.append(dataformat.AtomPositions(line[12:16].strip(), float(line[46:54])))
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
        self.atoms.append(atoms_each_timestep)
        print("[{} {} % ]".format('#' * int(100 / 10), int(100)))
        file.close()
        return

    def read_lgt(self):
        """Parses .lgt.

        Parses the lgtfile. If no lgtfile is defined falass will help the user to build one by working through the
        atom types in the pdb file and requesting input of the real and imaginary scattering lengths. This will also
        occur if a atom type if found in the pdbfile but not in the given lgts file. falass will write the lgtfile 
        to disk if atom types do not feature in the given lgtfile or one is written from scratch. 
        """
        if self.lgtfile:
            lines = line_count(self.lgtfile)
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
                    self.scat_lens.append(dataformat.ScatLens(line_list[0], float(line_list[1]), float(line_list[2])))
            print("[{} {} %]".format('#' * int(100 / 10), int(100)))
            file.close()
        else:
            raise ValueError("No lgtfile has been defined.")
        return

    def read_dat(self):
        """Parses .dat.

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
                        self.expdata.append(dataformat.QData(float(line_list[0]), float(line_list[1]),
                                                             float(line_list[1]) * (self.ierror / 100),
                                                             float(line_list[0]) * (self.resolution / 100)))
                    if len(line_list) == 3:
                        self.expdata.append(dataformat.QData(float(line_list[0]), float(line_list[1]),
                                                             float(line_list[2]),
                                                             float(line_list[0]) * (self.resolution / 100)))
                    if len(line_list) == 4:
                        self.expdata.append(
                            dataformat.QData(float(line_list[0]), float(line_list[1]), float(line_list[2]),
                                             float(line_list[3])))
            print("[{} {} %]".format('#' * int(100 / 10), int(100)))
            file.close()
        else:
            print('No DAT file has been given, therefore no comparison will be conducted, please use the get_qs '
                  'function. Alternatively the DAT file can be added using the setFile function.')
        return

    def get_qs(self, start=0.005, end=0.5, number=50):
        """Make custom q-vectors.

        If no datfile exists this function should be used to generate a linear-spaced range of q-vector for the
        calculation of the reflectometry over.

        Parameters
        ----------
        start: float, optional
            The first q-vector for which the reflectometry should be calculated.
        end: float, optional
            The last q-vector for which the reflectometry should be calculated.
        number: int, optional
            The number of q-vectors.
        """
        q_values = np.linspace(start, end, number)
        for i in range(0, len(q_values)):
            self.expdata.append(dataformat.QData(q_values[i], None, None, q_values[i] * (self.resolution / 100)))
        return


def check_duplicates(array, check):
    """Stops duplicate atom types.

    Checks if an atom type has already been added to an array.

    Parameters
    ----------
    array: array-type ScatLens
        The array to check.
    check: str
        The atom type to try and find.

    Returns
    -------
    bool
        True if the atom type is already present in the scatlen type array, false if not.
    """
    for i in range(0, len(array)):
        if array[i].atom == check:
            return True
    return False


def line_count(filename):
    """File length.

    Quickly counts the number of lines in a file

    Parameters
    ----------
    filename: str
        Name of the file that the number of lines is desired for.

    Returns
    -------
    int
        Number of lines in the file
    """
    i = 0
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
