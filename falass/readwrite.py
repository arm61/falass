import numpy as np
import sys
from multiprocessing import Pool

class pdbFile(object):
    def __init__(self, filename):
        """
        A class of the reading in of pdb files
        :param filename: str
            the pdb file name and path
        """
        self.filename = filename
        self.cell = []
        self.zpos = []
        self.atoms = []
        self.numberoftimesteps = 0
        self.times = []

    def readPDB(self):
        lines = line_count(self.filename)
        print("Reading PDB file \n[ 0 %]")
        file = open(self.filename, 'r')
        string = 0
        atoms_each_timestep = []
        zpos_each_timestep = []
        for i, line in enumerate(file):
            string_new = np.floor(i / lines * 100)
            if string_new > string + 9:
                string = string_new
                print("[{} {} %]".format('#' * int(string / 10), int(string)))
            if line[0:6] == "ATOM  ":
                atoms_each_timestep.append(line[12:16])
                zpos_each_timestep.append(float(line[46:54]))
            if "TITLE  " in line:
                if self.numberoftimesteps == 0:
                    self.numberoftimesteps += 1
                    self.times.append(float(line.split()[-1]))
                else:
                    self.atoms.append(atoms_each_timestep)
                    self.zpos.append(zpos_each_timestep)
                    atoms_each_timestep = []
                    zpos_each_timestep = []
                    self.numberoftimesteps += 1
                    self.times.append(float(line.split()[-1]))
            if "CRYST1  " in line:
                self.cell.append([float(line[6:15]), float(line[15:24]), float(line[24:33])])
        print("[{} {} %]".format('#' * int(100 / 10), int(100)))


def line_count(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1