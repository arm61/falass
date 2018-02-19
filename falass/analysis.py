import numpy as np
from falass import dataformat, readwrite
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
                atoms_each_molecule.append(get_atom_position(self.cell[self.number_of_timesteps - 1][2], line,
                                                             self.flip))
                if int(line[22:26]) == self.number_of_molecules + 1:
                    if self.number_of_molecules != 0:
                        atoms_each_timestep.append(atoms_each_molecule)
                        atoms_each_molecule = []
                    self.number_of_molecules += 1
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
                        self.tail_thickness.append(diffz)
                        self.chain_tilt_angles.append(np.rad2deg(theta))


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
        Object with atom type and z-position from the given line.
    """
    if flip:
        new_zpos = readwrite.flip_zpos(cell, float(line[46:54]))
        return dataformat.Atom3Positions(line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), new_zpos)
    else:
        return dataformat.Atom3Positions(line[12:16].strip(), line[30:38].strip(), line[38:46].strip(),
                                        float(line[46:54]))