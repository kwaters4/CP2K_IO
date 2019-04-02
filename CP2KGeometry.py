import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser
import glob
import os

class Cp2k_output:

    def __init__(self, structures, energies, completed, converged):
        self.structures = structures
        self.energies = energies
        self.completed = completed
        self.converged = converged

    def final_energy(self):
        return self.energies[-1]

    def final_structure(self):
        return self.structures[-1]

def parse_output(file_path):
    with open(file_path) as fp:
        lattices = []
        energies = []
        completed = False
        converged = False
        geo_opt = False
        cell_opt = False
        single_point = False
        for line in fp:
            if "CELL| Vector a " in line:
                curr_line = line.split()
                a1 = [float(curr_line[4]),float(curr_line[5]),float(curr_line[6])]
            if "CELL| Vector b " in line:
                curr_line = line.split()
                a2 = [float(curr_line[4]),float(curr_line[5]),float(curr_line[6])]
            if "CELL| Vector c " in line:
                curr_line = line.split()
                a3 = [float(curr_line[4]),float(curr_line[5]),float(curr_line[6])]
                lattices.append(np.array([a1,a2,a3]))	
            if "Total energy:" in line:
                curr_line = line.split()
                energies.append(float(curr_line[2]))
            if "PROGRAM ENDED AT" in line:
                completed = True
                converged = True
            if "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in line:
                converged = False
            if "SCF run NOT converged" in line:
                converged = False 
            if "GEO_OPT" in line:
                geo_opt = True
            if "CELL_OPT" in line:
                cell_opt = True
            if "ENERGY_FORCE" in line:
                single_point = True
    if geo_opt == True:
        for i in energies:
            lattices.append(lattices[-1])

    # Removes duplicates of first and last structures from output file 
    if cell_opt == True:
        lattices.pop(0)
        lattices.pop(0)
#        lattices.pop(-1)
    if len(lattices) > 2:
        lattices.pop(-1)
        lattices.pop(0)
    output_list = [lattices, energies, completed, converged, single_point]
    return output_list

def geometry_parse(file_path):
    structures = []
    with open(file_path) as fp:
        atoms = int(fp.readline().split()[0])
        xyz = []
        species = []
        coords = []
        specie = []
        for line in fp:
            if "i = " in line:
                continue
            if len(line.split()) == 1:
                continue
            curr_line = line.split()
            specie.append(curr_line[0])
            coords.append([float(curr_line[1]),float(curr_line[2]),float(curr_line[3])])
            if len(specie) % atoms == 0:
                xyz.append(coords)
                species.append(specie)
                coords = []
                specie = []
    return species, xyz 

def parse_structures(output_path, xyz_path, coords_are_cartesian = True):

    output_list = parse_output(output_path)
    if output_list[-1] == True:
        # If Single Point Calculations Read Original Cif File
        directory = os.path.dirname(output_path)
        filename = glob.glob('{}/*.cif'.format(directory))[0]
        structure = CifParser(filename).get_structures(primitive=False)[-1]
        run_output = Cp2k_output(structure,output_list[1],output_list[2],output_list[3])
        return run_output
    else:
        species, xyz = geometry_parse(xyz_path)

    structures=[]
    # TEMP FIX 
    while len(output_list[0]) > len(xyz):
        output_list[0].pop(0)
    if len(output_list[0]) != len(xyz) and len(output_list[0]) > 2:
        output_list[0].pop(0)
        print(len(output_list[0]),len(xyz))
    for i in range(0,len(output_list[0])):
        structure = Structure(output_list[0][i], species[i], xyz[i], coords_are_cartesian = coords_are_cartesian)
        structures.append(structure)
    run_output = Cp2k_output(structures,output_list[1],output_list[2],output_list[3])
    return run_output

def cp2k_structure(structure, filename):
    
    cp2k_xyz = open(filename, 'w')
    cp2k_xyz.write("\t&CELL\n")
    cp2k_xyz.write(" ABC {0:.6f} {1:.6f} {2:.6f} \n".format(structure.lattice.abc[0], structure.lattice.abc[1], structure.lattice.abc[2]))
    cp2k_xyz.write(" ALPHA_BETA_GAMMA {0:.6f} {1:.6f} {2:.6f}\n".format(structure.lattice.angles[0],structure.lattice.angles[1], structure.lattice.angles[2]))
    cp2k_xyz.write("\t&END CELL\n")
    cp2k_xyz.write("\t&COORD\n")
    for atom in structure.sites:
        cp2k_xyz.write("\t {0:s} {1:.6f} {2:.6f} {3:.6f} \n".format(atom.specie,atom.coords[0], atom.coords[1], atom.coords[2]))
    cp2k_xyz.write("\t&END COORD\n")

#xyz_path ="./Si_bulk8-pos-1.xyz"
#output_path ="./test.out"
#
#run = parse_structures(output_path, xyz_path)

#structure = Structure(np.array([a1,a2,a3]),atom_type,atom_list, coords_are_cartesian = True)
