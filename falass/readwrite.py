import numpy as np
import sys
from multiprocessing import Pool

class Files(object):
    def __init__(self, pdbfile, lgtsfile = None, datfile = None, resolution = 5., ierror = 5.):
        """
        A class of the reading in of pdb files
        :param filename: str
            the pdb file name and path
        """
        self.pdbfile = pdbfile
        self.cell = []
        self.zpos = []
        self.atoms = []
        self.number_of_timesteps = 0
        self.times = []
        self.lgtsfile = lgtsfile
        self.real_scat_lens = {}
        self.imag_scat_lens = {}
        self.elements = []
        self.datfile = datfile
        self.q = []
        self.i = []
        self.dq = []
        self.di = []
        self.ierror = ierror
        self.resolution = resolution

    def setFile(self, lgtsfile = None, datfile = None):
        self.lgtsfile = lgtsfile
        self.datfile = datfile

    def readPDB(self):
        lines = line_count(self.pdbfile)
        print("Reading PDB file \n[ 0 % ]")
        file = open(self.pdbfile, 'r')
        string = 0
        atoms_each_timestep = []
        zpos_each_timestep = []
        for i, line in enumerate(file):
            string_new = np.floor(i / lines * 100)
            if string_new > string + 9:
                string = string_new
                print("[{} {} % ]".format('#' * int(string / 10), int(string/10)*10))
            if line[0:6] == "ATOM  ":
                atoms_each_timestep.append(line[12:16].strip())
                zpos_each_timestep.append(float(line[46:54]))
            if "TITLE  " in line:
                if self.number_of_timesteps == 0:
                    self.number_of_timesteps += 1
                    self.times.append(float(line.split()[-1]))
                else:
                    self.atoms.append(atoms_each_timestep)
                    self.zpos.append(zpos_each_timestep)
                    atoms_each_timestep = []
                    zpos_each_timestep = []
                    self.number_of_timesteps += 1
                    self.times.append(float(line.split()[-1]))
            if "CRYST1  " in line:
                self.cell.append([float(line[6:15]), float(line[15:24]), float(line[24:33])])
        print("[{} {} %]".format('#' * int(100 / 10), int(100)))

    def readLGT(self):
        if self.lgtsfile:
            lines = line_count(self.lgtsfile)
            print("Reading LGT file \n[ 0 % ]")
            file = open(self.lgtsfile, 'r')
            string = 0
            for i, line in enumerate(file):
                string_new = np.floor(i / lines * 100)
                if string_new > string + 9:
                    string = string_new
                    print("[{} {} % ]".format('#' * int(string / 10), int(string/10)*10))
                line_list = line.split()
                self.real_scat_lens[line_list[0]] = float(line_list[1])
                self.imag_scat_lens[line_list[0]] = float(line_list[2])
            for i in range(0, len(self.atoms)):
                for j in range(0, len(self.atoms[i])):
                    if self.atoms[i][j] not in self.real_scat_lens:
                        real_scat_len = input('The following atom type has no scattering length given '
                                             'in the lgt file {} \nPlease define a real scattering length for '
                                              'this atom type: '.format(self.atoms[i][j]))
                        self.real_scat_lens[self.atoms[i][j]] = float(real_scat_len)
                    if self.atoms[i][j] not in self.imag_scat_lens:
                        imag_scat_len = input('The following atom type has no scattering length given '
                                              'in the lgt file {} \nPlease define a imaginary scattering length for '
                                              'this atom type: '.format(self.atoms[i][j]))
                        self.imag_scat_lens[self.atoms[i][j]] = float(imag_scat_len)
        else:
            print('There was no lgt file defined, falass will help you define one and save it for future use.')
            for i in range(0, len(self.atoms)):
                for j in range(0, len(self.atoms[i])):
                    if self.atoms[i][j] not in self.real_scat_lens:
                        real_scat_len = input('The following atom type has no scattering length given '
                                             'in the lgt file {} \nPlease define a real scattering length for '
                                              'this atom type: '.format(self.atoms[i][j]))
                        self.real_scat_lens[self.atoms[i][j]] = float(real_scat_len)
                    if self.atoms[i][j] not in self.imag_scat_lens:
                        imag_scat_len = input('The following atom type has no scattering length given '
                                              'in the lgt file {} \nPlease define a imaginary scattering length for '
                                              'this atom type: '.format(self.atoms[i][j]))
                        self.imag_scat_lens[self.atoms[i][j]] = float(imag_scat_len)

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
                        self.q.append(float(line_list[0]))
                        self.i.append(float(line_list[1]))
                        self.di.append(float(line_list[1]) * (self.ierror / 100))
                        self.dq.append(float(line_list[0]) * (self.resolution / 100))
                    if len(line_list) == 3:
                        self.q.append(float(line_list[0]))
                        self.i.append(float(line_list[1]))
                        self.di.append(float(line_list[2]))
                        self.dq.append(float(line_list[0]) * (self.resolution / 100))
                    if len(line_list) == 4:
                        self.q.append(float(line_list[0]))
                        self.i.append(float(line_list[1]))
                        self.di.append(float(line_list[2]))
                        self.dq.append(float(line_list[3]))
        else:
            print('No DAT file has been given, therefore no comparison will be conducted. '
                  'Alternatively the DAT file can be added using the setFile function.')



def line_count(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1