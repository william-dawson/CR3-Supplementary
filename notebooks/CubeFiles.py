"""
This module contains some helper routines for dealing with cube files.

Note: this was only made for the very specific data/task herein.
"""

def cube_name(iorb, kpt=1, wf_type="wavefunction"):
    """
    Creates the name of the cubefile from the indices.
    
    Args:
        iorb (int): orbital index (1 based)
        kpt (int): k point index (1 based)
        wf_type (str): can be `wavefunction`, `virtuals`, or `minBasis`
    
    Returns:
        (str): the name of the file.
    """
    wfn = wf_type + "-"
    wfn += "k" + "{:03d}".format(kpt)
    wfn += "-NR."
    wfn += "b" + "{:06d}".format(iorb)
    
    return wfn

def generate_cube(log, iorb, kpt=1, temp_dir=".", skip=False):
    """
    Generate the cubefile for a given orbital and k point.
    
    Args:
        log (BigDFT.Logfiles.Logfile): log from the run
        iorb (int): orbital index (1 based)
        kpt (int): k point index (1 based)
        temp_dir (str): temporary data for storing the intermediary files, and the result.
        skip (logical): will skip if the file already exists
    Returns:
        (str): path to the generated cube file.
    """
    from os import system, environ, getcwd, chdir
    from os.path import join, basename, exists
    from futile.Utils import ensure_dir
    from shutil import copyfile
    
    # Generate the orbital file name
    wfn = cube_name(iorb, kpt)
    orbf = join(log.srcdir, log.data_directory, wfn)
    
    # Generate the output file name
    cubf = join(temp_dir, wfn) + ".cube"
    
    # Skip option
    if exists(cubf):
        return cubf
    
    # Ensure temporary directory
    ensure_dir(temp_dir)
    
    # Copy the input file to the temporary directory
    base = basename(log.label)
    copyfile(log.label, join(temp_dir, base))
    
    basedir = getcwd()
    try:
        chdir(temp_dir)
        
        # Generate Command
        cmd = environ.get("BIGDFT_MPIRUN", "")
        cmd += " " + join(environ.get("BIGDFT_ROOT", ""), "bigdft-tool")
        cmd += " -a export-wf"
        cmd += " " + orbf + " --name=" + basename(log.label).replace(".yaml", "")
        cmd += " >/dev/null 2>&1"
        
        # Run
        system(cmd)
    finally:
        chdir(basedir)
        
    return cubf

class CubeFile:
    """
    A class for encapsulating cubefiles.    
    """
    def __init__(self):
        # (BigDFT.Systems.System): The system associated with the points on teh grid.
        self.sys = None

    def read(self, ifile):
        """
        Read in from file.
        
        Args:
            ifile (TextIOBase): file to read from.
        """
        from numpy import zeros
        self.grid_points = []
        self.spacing = []
        
        self.comment1 = next(ifile).strip()
        self.comment2 = next(ifile).strip()
        split = next(ifile).split()
        self.natoms = int(split[0])
        self.origin = [float(x) for x in split[1:]]
        if self.natoms < 0:
            wf = True
        else:
            wf = False
        self.natoms = abs(self.natoms)
        for i in range(3):
            split = next(ifile).split()
            self.grid_points.append(int(split[0]))
            self.spacing.append([float(x) for x in split[1:]])
            
        self.data = zeros(self.grid_points)
        self._read_atoms(ifile)
        if wf:
            next(ifile)
        split = next(ifile).split()
        for i in range(self.grid_points[0]):
            for j in range(self.grid_points[1]):
                for k in range(self.grid_points[2]):
                    if len(split) == 0:
                        split = next(ifile).split()
                    self.data[i, j, k] = float(split[0])
                    split = split[1:]
                    
    def _read_atoms(self, ifile):
        from BigDFT.Systems import System
        from BigDFT.Fragments import Fragment
        from BigDFT.Atoms import Atom, number_to_symbol
        
        self.sys = System()
        for i in range(self.natoms):
            split = next(ifile).split()
            sym = number_to_symbol(int(split[0]))
            nel = float(split[1])
            pos = [float(x) for x in split[2:]]
            at = Atom({sym: pos, "units": "bohr", "nzion": nel})
            self.sys["FRA:"+str(i)] = Fragment([at])
            
    def dot(self, b):
        """
        Compute the dot product between two cubefiles.
        
        Args:
            b (CubeFiles.CubeFile): the second cube.
            
        Returns:
            (float): the dot product
        """
        from numpy import multiply, sum
        return sum(multiply(self.data, b.data), axis=(0, 1, 2))
        
    def write(self, ofile):
        """
        Write out to file.
        
        Args:
            ofile (TextIOBase): file to write to.
        """
        # Header
        ofile.write(self.comment1 + "\n")
        ofile.write(self.comment2 + "\n")
        
        # Box
        ofile.write(str(self.natoms) + "\t" + 
                    "\t".join([str(x) for x in self.origin]) + "\n")
        for i in range(3):
            ofile.write(str(self.grid_points[i]) + "\t")
            ofile.write("\t".join([str(x) for x in self.spacing[i]]))
            ofile.write("\n")
            
        # Write out the atoms
        for frag in self.sys.values():
            for at in frag:
                ofile.write(str(at.atomic_number) + "\t" + str(at.nel) + "\t" + 
                            "\t".join([str(x) for x in at.get_position("angstroem")]) + "\n")

        # Write out the data
        counter = 0
        for i in range(self.grid_points[0]):
            for j in range(self.grid_points[1]):
                for k in range(self.grid_points[2]):
                    ofile.write(str(self.data[i, j, k]) + " ")
                    if counter % 6 == 5:
                        ofile.write("\n")
                    counter += 1