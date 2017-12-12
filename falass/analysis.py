import numpy as np
from falass import dataformat, readwrite

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
        for i, line in enumerate(file):
            percentage_new = np.floor(i / lines * 100)
            percentage = readwrite.check_update(percentage, percentage_new)
            if line[0:6] == "ATOM  ":
                atoms_each_timestep.append(get_atom_position(self.cell[self.number_of_timesteps - 1][2], line,
                                                             self.flip))
            if "TITLE  " in line:
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