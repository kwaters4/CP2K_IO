#! /usr/bin/python
import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
from pathlib import Path
import sys

def read_car_file(filename):
	'''
	Arguments: filename of the .car file to be read.

	Returns: pymatgen structure class

	Reads .in car file and parses the file (not tested for all formats) 
	'''
	with open(filename) as fp:
		coords = []
		species = []
		atoms = False
		for line in fp:
			if "end" in line:
				atoms = False
			if atoms == True:
				curr_line = line.split()
				specie = curr_line[0]
				xyz = [float(curr_line[1]),float(curr_line[2]),float(curr_line[3])]
				species.append(specie)
				coords.append(xyz)
			if "PBC " in line:
				curr_line = line.split()
				abc = [float(curr_line[1]),float(curr_line[2]),float(curr_line[3])]
				angles = [float(curr_line[4]),float(curr_line[5]),float(curr_line[6])]
				atoms = True
        lattice = Lattice.from_lengths_and_angles(abc,angles)
        structure = Structure(lattice, species, coords, coords_are_cartesian=True)
    return structure

# Usage:
# $python car_reader.py FILENAME.car
# Writes FILENAME.cif

# First Argument given to python script is filename
filename = sys.argv[1]
# reads file
structure = read_car_file(filename)
# Writes .cif file of converted structure
CifWriter(structure.get_sorted_structure()).write_file(Path(filename).with_suffix('.cif'))
