import numpy as np
from falass import dataformat, readwrite
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

class Analysis:
    """Trajectory analysis

    This is a special class for the analysis of simulation trajectories.

    Parameters
    ----------
    files: falass.readwrite.Files
        Contains information of the files for analysis.
    """
    def __init__(self, files):
        self.files = files
        self.cell = []
        self.atoms = []
        self.number_of_timesteps = 0
        self.times = []
        self.flip = files.flip
        self.number_of_molecules = 0
        self.chain_tilt_angles = []
        self.tail_thickness = []
        self.head_bin = []
        self.water_bin = []
        self.tail_bin = []
        self.bin_width = 1.
        self.all_head_bin = []
        self.all_water_bin = []
        self.all_tail_bin = []
        self.wph = []
        self.z = []
        self.average_wph = []
        self.tail_length = []
        self.head_thickness = []
        self.dda = []
        self.cca = []

    def read_full_pdb(self):
        """Parse 3 coordinate pdb.

        Similar to readwrite.read_pdb however instead of only reading a single z-position it reads the 3d position.
        """
        lines = readwrite.line_count(self.files.pdbfile)
        print("Reading PDB file \n")
        file = open(self.files.pdbfile, 'r')
        percentage = 0
        readwrite.print_update(percentage)
        atoms_each_timestep = []
        atoms_each_molecule = []
        for i, line in enumerate(file):
            percentage_new = np.floor(i / lines * 100)
            percentage = readwrite.check_update(percentage, percentage_new)
            if line[0:6] == "ATOM  ":
                if int(line[22:26]) == self.number_of_molecules + 1:
                    if self.number_of_molecules != 0:
                        atoms_each_timestep.append(atoms_each_molecule)
                        atoms_each_molecule = []
                    self.number_of_molecules += 1
                atoms_each_molecule.append(get_atom_position(self.cell[self.number_of_timesteps - 1][2], line,
                                                             self.flip))
            if "TITLE  " in line:
                self.number_of_molecules = 0
                if self.number_of_timesteps == 0:
                    self.number_of_timesteps, new_time = readwrite.iterate_time(self.number_of_timesteps, line)
                    self.times.append(new_time)
                else:
                    self.atoms.append(atoms_each_timestep)
                    atoms_each_timestep = []
                    self.number_of_timesteps, new_time = readwrite.iterate_time(self.number_of_timesteps, line)
                    self.times.append(new_time)
            if "CRYST1  " in line:
                self.cell.append(readwrite.get_cell_parameters(line))
        self.atoms.append(atoms_each_timestep)
        readwrite.print_update(100)
        file.close()
        return

    def get_number_density(self, labels):
        for i in range(0, len(self.atoms)):
            for j in range(0, len(self.atoms[i])):
                for k in range(0, len(labels)):
                    if self.atoms[i][j].atom == labels[k]:
                        print(self.atoms[i][j].atom)

    def get_chain_tilt(self, labels):
        self.chain_tilt_angles = []
        self.tail_thickness = []
        self.tail_length = []
        for i in range(0, len(self.atoms)):
            for j in range(0, len(self.atoms[i])):
                for l in range(0, len(labels)):
                    tail_atoms_x = []
                    tail_atoms_y = []
                    tail_atoms_z = []
                    for k in range(0, len(self.atoms[i][j])):
                        if self.atoms[i][j][k].atom in labels[l]:
                            tail_atoms_x.append(np.float(self.atoms[i][j][k].x))
                            tail_atoms_y.append(np.float(self.atoms[i][j][k].y))
                            tail_atoms_z.append(np.float(self.atoms[i][j][k].z))
                    if len(tail_atoms_z) != 0:
                        diffx = np.abs(tail_atoms_x[1] - tail_atoms_x[0])
                        diffy = np.abs(tail_atoms_y[1] - tail_atoms_y[0])
                        diffz = np.abs(tail_atoms_z[1] - tail_atoms_z[0])
                        if diffx > self.cell[i][0]/2:
                            diffx = np.abs(diffx - self.cell[i][0])
                        if diffy > self.cell[i][1]/2:
                            diffy = np.abs(diffy - self.cell[i][1])
                        if diffz > self.cell[i][2]/2:
                            diffz = np.abs(diffz - self.cell[i][2])
                        xy = np.sqrt(np.square(diffx) +
                                     np.square(diffy))
                        theta = np.arctan(xy/diffz)
                        dist = np.sqrt(np.square(diffx) +
                                     np.square(diffy) + np.square(diffz))
                        self.tail_length.append(dist)
                        self.tail_thickness.append(diffz)
                        self.chain_tilt_angles.append(np.rad2deg(theta))

    def get_number_density(self, head_labels, tail_labels, water_labels, bin_width, offset = 0):
        self.all_head_bin = []
        self.all_water_bin = []
        self.all_tail_bin = []
        for i in range(0, len(self.atoms)):
            water_bin = np.zeros(int(np.ceil(self.cell[0][2] / bin_width)))
            head_bin = np.zeros(int(np.ceil(self.cell[0][2] / bin_width)))
            tail_bin = np.zeros(int(np.ceil(self.cell[0][2] / bin_width)))
            for j in range(0, len(self.atoms[i])):
                for k in range(0, len(self.atoms[i][j])):
                    if self.atoms[i][j][k].atom in head_labels:
                        head_bin[int(self.atoms[i][j][k].z / bin_width)] += 1
                    elif self.atoms[i][j][k].atom in water_labels:
                        water_bin[int(self.atoms[i][j][k].z / bin_width)] += 1
                    elif self.atoms[i][j][k].atom in tail_labels:
                        tail_bin[int(self.atoms[i][j][k].z / bin_width)] += 1
            head_bin = head_bin / (len(head_labels) * self.cell[0][0] * self.cell[0][1] * bin_width)
            self.all_head_bin.append(head_bin)
            water_bin = water_bin / (len(water_labels) * self.cell[0][0] * self.cell[0][1] * bin_width)
            self.all_water_bin.append(water_bin)
            tail_bin = tail_bin / (len(tail_labels) * self.cell[0][0] * self.cell[0][1] * bin_width)
            self.all_tail_bin.append(tail_bin)
        self.head_bin = get_average_over_time(self.all_head_bin)
        self.water_bin = get_average_over_time(self.all_water_bin)
        self.tail_bin = get_average_over_time(self.all_tail_bin)
        self.z = np.arange(0, len(self.head_bin)) * bin_width - offset

    def get_wph(self, bin_width):
        self.head_thickness = []
        self.dda = []
        self.cca = []
        for i in range(0, len(self.all_head_bin)):
            summ = []
            summm = 0
            for j in range(0, len(self.all_head_bin[i])):
                summm += self.all_head_bin[i][j]
                summ.append(summm)
            aa = 0.05 * summm
            bb = 0.95 * summm
            cc = get_point(self.z, summ, aa)
            dd = get_point(self.z, summ, bb)
            ht = dd - cc
            self.dda.append(dd)
            self.cca.append(cc)
            self.head_thickness.append(ht)
        self.wph = np.zeros(len(self.all_head_bin))
        for i in range(0, len(self.dda)):
            a = 0
            b = 0
            gs = []
            for j in range(0, len(self.z)):
                if self.z[j] >= self.cca[i] and self.z[j] < self.dda[i] + bin_width:
                    a += self.all_head_bin[i][j]
                    b += self.all_water_bin[i][j]
            re = b / a
            self.wph[i] = np.sum(re)

    def plot_number_density(self, size, colors=['green', 'red', 'blue']):
        fig, ax1 = plt.subplots(figsize=(size[0], size[1]))
        ax1.plot(self.z, np.asarray(self.head_bin) * 10000, color = colors[0])
        ax1.plot(self.z, np.asarray(self.tail_bin) * 10000, color = colors[1])
        for i in range(0, len(self.all_water_bin)):
            ax1.plot(self.z, np.asarray(self.all_head_bin[i]) * 10000, color = colors[0], alpha = 0.0075)
            ax1.plot(self.z, np.asarray(self.all_tail_bin[i]) * 10000, color = colors[1], alpha=0.0075)
        a = [np.max(self.all_tail_bin), np.max(self.all_head_bin)]
        ax1.set_ylim([0., np.max(a) * 10000])
        ax2 = ax1.twinx()
        ax2.set_ylim([0., np.max(self.all_water_bin) * 1000])
        ax2.plot(self.z, np.asarray(self.water_bin) * 1000, color = colors[2])
        for i in range(0, len(self.all_water_bin)):
            ax2.plot(self.z, np.asarray(self.all_water_bin[i]) * 1000, color = colors[2], alpha = 0.0075)

        return fig, ax1, ax2


def get_atom_position(cell, line, flip):
    """Get the atom position and type

    Reads the text line from the pdb and assigns it to an falass.dataformat.AtomPosition type object.

    Parameters
    ----------
    cell: float
        z-cell dimension.
    line: str
        Line from pdb file.
    flip: bool
        Should the cell be flipped in the z-dimension.

    Returns
    -------
    falass.dataformat.AtomPositions
        Object with atom type and x, y, z-position from the given line.
    """
    if flip:
        new_zpos = readwrite.flip_zpos(cell, float(line[46:54]))
        return dataformat.Atom3Positions(line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), new_zpos)
    else:
        return dataformat.Atom3Positions(line[12:16].strip(), line[30:38].strip(), line[38:46].strip(),
                                        float(line[46:54]))

def get_average_over_time(bin_array):
    average_bin_array = []
    bin_array = np.transpose(bin_array)
    for time in range(0, len(bin_array)):
        average_bin_array.append(np.average(bin_array[time]))
    return average_bin_array

def get_point(x, y, target):
    for i in range(0, len(x) - 1):
        if y[i] <= target and y[i+1] > target:
            m = (y[i+1] - y[i]) / (x[i+1] - x[i])
            c = m * x[i] - y[i]
            x = (target + c) / m
            return x
