import numpy as np
import os
from multiprocessing import Pool
from falass import dataformat

class Files(object):
    def __init__(self, pdbfile, lgtsfile = None, datfile = None, resolution = 5., ierror = 5.):
        """
        A class of the reading in of pdb files
        :param filename: str
            the pdb file name and path
        """
        self.pdbfile = pdbfile
        self.cell = []
        self.atoms = []
        self.number_of_timesteps = 0
        self.times = []
        self.lgtsfile = lgtsfile
        self.scat_lens = []
        self.real_scat_lens = {}
        self.imag_scat_lens = {}
        self.elements = []
        self.datfile = datfile
        self.expdata = []
        self.ierror = ierror
        self.resolution = resolution

    def setFile(self, pdbfile = None, lgtsfile = None, datfile = None):
        if pdbfile:
            self.pdbfile = pdbfile
        if lgtsfile:
            self.lgtsfile = lgtsfile
        if datfile:
            self.datfile = datfile

    def readPDB(self):
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
        print("[{} {} %]".format('#' * int(100 / 10), int(100)))

    def readLGT(self):
        lgtsfile_name = None
        if self.lgtsfile:
            lines = line_count(self.lgtsfile)
            lgtsfile_name = self.lgtsfile
            print("Reading LGT file \n[ 0 % ]")
            file = open(self.lgtsfile, 'r')
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
            lgtsfile_name = input("What should the lgt file be named? ")
            path, extension = os.path.splittext(output_name)
            if extension != '.lgt':
                lgtsfile_name = path + '.lgt'
        lgtsf = open(lgtsfile_name, 'w')
        for i in range(0, len(self.scat_lens)):
            lgtsf.write('{} {} {}\n'.format(self.scat_lens[i].atom, self.scat_lens[i].real * 1e5,
                                            self.scat_lens[i].imag * 1e5))


    def readDAT(self):
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
        else:
            print('No DAT file has been given, therefore no comparison will be conducted. '
                  'Alternatively the DAT file can be added using the setFile function.')

def check_duplicates(array, check):
    for i in range(0, len(array)):
        if array[i].atom == check:
            return True
    return False

def line_count(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1