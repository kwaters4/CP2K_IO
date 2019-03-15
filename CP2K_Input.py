import CP2KGeometry as cp2kio
import datetime
import subprocess


def bash_command(cmd):
	subprocess.call(cmd, shell=True)
	return

class CP2K_Input_Deck:

	def __init__(self, structure, directory, name = "default", method = "QS", functional = "pbe", basis_set = "gth_basis_sets", potential = "potential", vdw_potential = "pair_potential"):

		self.name = name
		self.directory = directory
		self.structure = structure
		self.Global = Global(name)
		self.Motion = Motion(self.Global.run_type)
		self.Force_Eval = Force_Eval(structure, method, functional, basis_set, potential, vdw_potential)

	def write_file(self):

		filename = "{}.inp".format(self.name)
		input_deck = open("{}/{}.inp".format(self.directory, self.name), "w")
		# GLOBAL
		self.Global.write_to_file(input_deck)
		# MOTION
		self.Motion.write_to_file(input_deck)
		# FORCE_EVAL
		self.Force_Eval.write_to_file(input_deck)
		input_deck.close()
	
	def shell_script(self, queue = "debug", nodes = 1, time = 60):

		filename = "{}.sh".format(self.name)
		sub_script = open("{}/{}.sh".format(self.directory, self.name), "w")

		sub_script.write("#! /bin/bash \n\n")	
		sub_script.write("#PBS -A ARLAP02333700\n")
		sub_script.write("#PBS -q {}\n".format(queue))
		sub_script.write("#PBS -l select={}:ncpus={}:mpiprocs={}\n".format(nodes, nodes*48, nodes*48))
		sub_script.write("#PBS -l walltime={}\n\n".format(str(datetime.timedelta(minutes = time))))

		sub_script.write("module unload compiler/intel\n")
		sub_script.write("module load gcc\n")
		sub_script.write("module load mpt\n")
		sub_script.write("module load costinit\n")
		sub_script.write("module load fftw3-mpi/gnu/sgimpt/3.3.5\n\n")

		sub_script.write("INPUT={}\n".format(self.name))
		sub_script.write("cd {}\n\n".format(self.directory))

		sub_script.write('export CP2K_DATA_DIR="/app/ccm4/CP2K/cp2k_6_1/data"\n')
		sub_script.write('CP2KLOC="/app/ccm4/CP2K/cp2k_6_1"\n')
		sub_script.write('LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CP2KLOC/lib/Mustang-elpa-gf_mkl18_smm/popt\n')
		sub_script.write('LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/p/app/PET/pkgs/CP2K/lib/elpa-2016.05.003/lib\n')
		sub_script.write('CP2K="$CP2KLOC/exe/Mustang-elpa-gf_mkl18_smm"\n\n')
		
		sub_script.write("mpiexec_mpt -np {} $CP2K/cp2k.popt -o {}.out {}.inp > dump.txt".format(nodes*48, self.name, self.name))
		sub_script.close()

	def submit_job(self):
		print("{}/{}.sh".format(self.directory, self.name))
		bash_command("qsub {}/{}.sh".format(self.directory, self.name))
		

class Global:

	blacs_grid = "SQUARE" #COLUMN, ROW, SQUARE
	preferred_fft_library = "FFTW3"
	extended_fft_lengths = ".True."
	print_level = "MEDIUM"
	run_type = "CELL_OPT" # "GEO_OPT" # "ENERGY_FORCE"

	def __init__(self, name):
		self.name = name
	
	def print_options(self):
		print("--------------------------------------------")
		print("Global Class")
		print("--------------------------------------------")
		print("Name                   : {}".format(self.name))
		print("Run Type               : {}".format(self.run_type))
		print("Print Level            : {}".format(self.print_level))
		print("Preferred FFT Library  : {}".format(self.fft_library))
		print("Extended FFT Lengths   : {}".format(self.fft_lengths))
		print("Blacs grid             : {}".format(self.blacs_grid))
		print("--------------------------------------------")

	def asdict(self):
		return {"PROJECT_NAME" : self.name, 
			"RUN_TYPE" : self.run_type,
			"PRINT_LEVEL" : self.print_level,
			"BLACS_GRID" : self.blacs_grid,
			"PREFERRED_FFT_LIBRARY" : self.preferred_fft_library,
			"EXTENDED_FFT_LENGTHS" : self.extended_fft_lengths,
			}

	def write_to_file(self, input_deck):
		input_deck.write("&GLOBAL \n")
		for key, value in self.asdict().items():
			input_deck.write("\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("&END GLOBAL\n")

class Motion:
	
	
	def __init__(self, run_type):
		self.run_type = run_type
		if run_type == "CELL_OPT":
			self.Cell_Opt = Cell_Opt()
			self.Geo_Opt = Geo_Opt()
		if run_type == "GEO_OPT":
			self.Geo_Opt = Geo_Opt()

	def write_to_file(self, input_deck):
		input_deck.write("&MOTION\n")
		if self.run_type == "CELL_OPT":
			self.Cell_Opt.write_to_file(input_deck)
		if self.run_type == "CELL_OPT" or self.run_type == "GEO_OPT":
			self.Geo_Opt.write_to_file(input_deck)
		input_deck.write("&END MOTION\n")

class Cell_Opt:

	constraint = "None" # X, XY, XZ, Y, YZ, Z 
	external_pressure = 100 # Default 100 bar
	keep_angles = ".True."
	keep_symmetry = ".False."
	max_dr = 3E-3 
	max_force = 1.00E-04	
	max_iter = 1800
	optimizer = "BFGS" #BFGS, CG, LBFGS
	pressure_tol = 100
	rms_dr = 1.5E-3 # Default bohr
	rms_force = 3E-4 # Default bohr/hatree
	step_start_val = 1 #0
	type = "direct_cell_opt" # DIRECT_CELL_OPT, GEO_OPT, MD

	def asdict(self):
		return {"CONSTRAINT" : self.constraint,
			"EXTERNAL_PRESSURE" : self.external_pressure,
			"KEEP_ANGLES": self.keep_angles,
			"KEEP_SYMMETRY" : self.keep_symmetry,
			"MAX_DR" : self.max_dr,
			"MAX_FORCE" : self.max_force,
			"MAX_ITER" : self.max_iter,
			"OPTIMIZER" : self.optimizer, 
			"PRESSURE_TOLERANCE" : self.pressure_tol,
			"RMS_DR" : self.rms_dr,
			"RMS_FORCE" : self.rms_force,
			"STEP_START_VAL" : self.step_start_val,
			"TYPE" : self.type
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Cell_Opt Class")
		print("--------------------------------------------")
		print("Constraint     : {}".format(self.constraint))
		print("External Pres. : {}".format(self.external_pressure))
		print("Keep Angles    : {}".format(self.keep_angles))
		print("Keep Symmetry  : {}".format(self.keep_symmetry))
		print("Max DR.        : {}".format(self.max_dr))
		print("Max Iter.      : {}".format(self.max_iter))
		print("Max Force      : {}".format(self.max_force))
		print("Optimizer      : {}".format(self.optimizer))
		print("Pressure Tol.  : {}".format(self.pressure_tol))
		print("RMS DR         : {}".format(self.rms_dr))
		print("RMS Force      : {}".format(self.rms_force))
		print("Start Value    : {}".format(self.step_start_val))
		print("Type           : {}".format(self.type))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t&CELL_OPT\n")
		for key, value in self.asdict().items():
			input_deck.write("\t \t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t&END CELL_OPT\n")

class Geo_Opt:

	max_dr = 3E-3 
	max_force = 1.00E-04	
	max_iter = 1800
	optimizer = "BFGS" #BFGS, CG, LBFGS
	rms_dr = 1.5E-3 # Default bohr
	rms_force = 3E-4 # Default bohr/hatree
	step_start_val = 1 #0
	type = "minimization" # minimization, transition state

	def asdict(self):
		return {"MAX_DR" : self.max_dr,
			"MAX_FORCE" : self.max_force,
			"MAX_ITER" : self.max_iter,
			"RMS_DR" : self.rms_dr,
			"RMS_FORCE" : self.rms_force,
			"STEP_START_VAL" : self.step_start_val,
			"TYPE" : self.type
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Geo_Opt Class")
		print("--------------------------------------------")
		print("Max DR.        : {}".format(self.max_dr))
		print("Max Iter.      : {}".format(self.max_iter))
		print("Max Force      : {}".format(self.max_force))
		print("RMS DR         : {}".format(self.rms_dr))
		print("RMS Force      : {}".format(self.rms_force))
		print("Start Value    : {}".format(self.step_start_val))
		print("Type           : {}".format(self.type))
		print("--------------------------------------------")
	
	def write_to_file(self, input_deck):
		input_deck.write("\t&GEO_OPT\n")
		for key, value in self.asdict().items():
			input_deck.write("\t \t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t&END GEO_OPT\n")
	
class Force_Eval:

	stress_tensor = "Analytical" # ANALYTICAL, DIAGONAL_ANALYTICAL, DIAGONAL NUMERICAL, NONE, NUMERICAL 

	def __init__(self, structure, method, functional, basis_set, potential, vdw_potential):
		self.method = method
		self.Dft = Dft(method, functional, basis_set, potential, vdw_potential)
		self.Sub_Sys = Sub_Sys(structure, basis_set, potential)
	
	def asdict(self):
		return {"METHOD" : self.method, 
			"STRESS_TENSOR" : self.stress_tensor,
			}

	def write_to_file(self, input_deck):
		input_deck.write("&FORCE_EVAL\n")
		for key, value in self.asdict().items():
			input_deck.write("\t{} {}\n".format(key,str(value).upper()))
		self.Dft.write_to_file(input_deck)
		self.Sub_Sys.write_to_file(input_deck)
		input_deck.write("&END FORCE_EVAL\n")

class Dft:

	charge = 0
	excitations = "none" # NONE, TDDFPT, TDLR
	multiplicity = 0
	plus_u_method = "mulliken" 
	relax_multiplicity = 0
	roks = ".False."
	subcells = 2
	surface_dipole_correction = ".False."
	surf_dip_dir = "Z"
	uks = ".False."
	wfn_restart_file_name = ""

	def __init__(self, method, functional, basis_set, potential, vdw_potential):
		self.basis_set_file_name = basis_set
		self.potential_file_name = potential
		self.vdw_potential = vdw_potential
		self.method = method
		self.Scf = Scf()
		if self.method == "QS":
			self.Qs = Qs() 
		self.Mgrid = Mgrid()
		self.Xc = Xc(functional, vdw_potential)

	def asdict(self):
		return {"BASIS_SET_FILE_NAME" : self.basis_set_file_name, 
			"POTENTIAL_FILE_NAME" : self.potential_file_name,
			"CHARGE" : self.charge,
			"EXCITATIONS" : self.excitations,
			"MULTIPLICITY" : self.multiplicity,
			"PLUS_U_METHOD" : self.plus_u_method,
			"RELAX_MULTIPLICITY" : self.relax_multiplicity,
			"ROKS" : self.roks,
			"SUBCELLS" : self.subcells,
			"SURFACE_DIPOLE_CORRECTION" : self.surface_dipole_correction,
			"SURF_DIP_DIR" : self.surf_dip_dir,
			"UKS" : self.uks,
			"WFN_RESTART_FILE_NAME" : self.wfn_restart_file_name,
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Dft Class")
		print("--------------------------------------------")
		print("Basis Set File       : {}".format(self.basis_set_file_name))
		print("Potential File       : {}".format(self.potential_file_name))
		print("Charge               : {}".format(self.charge))
		print("Excitations          : {}".format(self.excitations))
		print("Multiplicity         : {}".format(self.multiplicity))
		print("Plus U Methods       : {}".format(self.plus_u_method))
		print("Relax Mult.          : {}".format(self.relax_multiplicity))
		print("Res. Open KS         : {}".format(self.roks))
		print("Subcells             : {}".format(self.subcells))
		print("Surface Dipole Corr. : {}".format(self.surface_dipole_correction))
		print("Surface Dipole Dir.  : {}".format(self.surf_dip_dir))
		print("Spin Polarized KS    : {}".format(self.uks))
		print("WFN Restart File     : {}".format(self.wfn_restart_file_name))
		print("--------------------------------------------")
	
	def write_to_file(self, input_deck):
		input_deck.write("\t&DFT\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t{} {}\n".format(key,str(value).upper()))
		self.Scf.write_to_file(input_deck)
		if self.method == "QS":
			self.Qs.write_to_file(input_deck)
		self.Mgrid.write_to_file(input_deck)
		self.Xc.write_to_file(input_deck)
		input_deck.write("\t&END DFT\n")
	
class Scf:

	max_scf = 200
	max_diis = 8
	eps_scf = 1E-6 
	eps_diis = 5E-1
	scf_guess = "atomic" # ATOMIC, CORE, HISTORY_RESTART, MOPOC, NONE, RANDOM, RESTART, SPARSE

	def __init__(self):
		self.Mixing = Mixing()

	def asdict(self):
		return {"MAX_SCF" : self.max_scf, 
			"MAX_DIIS" : self.max_diis,
			"EPS_SCF" : self.eps_scf,
			"EPS_DIIS" : self.eps_diis,
			"SCF_GUESS" : self.scf_guess,
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("SCF Class")
		print("--------------------------------------------")
		print("Max SCF              : {}".format(self.max_scf))
		print("Max DIIS             : {}".format(self.max_diis))
		print("EPS SCF              : {}".format(self.eps_scf))
		print("EPS DIIS             : {}".format(self.eps_diis))
		print("SCF Guess            : {}".format(self.scf_guess))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t\t&SCF\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t{} {}\n".format(key,str(value).upper()))
		self.Mixing.write_to_file(input_deck)
		input_deck.write("\t\t&END SCF\n")
	
		return
class Mixing:

	mixing = ".TRUE."
	alpha = 2E-1 

	def asdict(self):
		return {"&MIXING" : self.mixing,
			"ALPHA" : self.alpha, 
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Mixing Class")
		print("--------------------------------------------")
		print("Alpha                : {}".format(self.alpha))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t\t\t&END MIXING \n")

class Qs:

	eps_default = 9E-13
	map_consistent= ".True."
	extrapolation = "ps"
	extrapolation_order = 3

	def asdict(self):
		return {"EPS_DEFAULT" : self.eps_default, 
			"MAP_CONSISTENT" : self.map_consistent,
			"EXTRAPOLATION" : self.extrapolation,
			"EXTRAPOLATION_ORDER" : self.extrapolation_order,
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Mixing Class")
		print("--------------------------------------------")
		print("EPS Default           : {}".format(self.eps_default))
		print("Map Consistent        : {}".format(self.map_consistent))
		print("Extrapolation         : {}".format(self.extrapolation))
		print("Extrapolation Order   : {}".format(self.extrapolation_order))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t\t&QS\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t&END QS \n")

class Mgrid:

	ngrids = 4
	cutoff = 5E2 

	def asdict(self):
		return {"NGRIDS" : self.ngrids, 
			"CUTOFF" : self.cutoff,
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Mgrid Class")
		print("--------------------------------------------")
		print("NGrids           : {}".format(self.ngrids))
		print("Cutoff           : {}".format(self.cutoff))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t\t&MGRID\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t&END MGRID \n")

class Xc:

	density_cutoff = 1E-10
	gradient_cutoff = 1E-10
	tau_cutoff = 1E-10

	def __init__(self, functional, vdw_potential):
		self.Xc_Functional = Xc_Functional(functional)
		self.Vdw_Potential = Vdw_Potential(vdw_potential, functional)
	
	def asdict(self):
		return {"DENSITY_CUTOFF" : self.density_cutoff, 
			"GRADIENT_CUTOFF" : self.gradient_cutoff,
			"TAU_CUTOFF" : self.tau_cutoff,
			} 
	def print_options(self):
		print("--------------------------------------------")
		print("Xc Class")
		print("--------------------------------------------")
		print("Density Cutoff          : {}".format(self.density_cutoff))
		print("Gradient Cutoff         : {}".format(self.gradient_cutoff))
		print("Tau Cutoff              : {}".format(self.tau_cutoff))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t\t&XC\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t{} {}\n".format(key,str(value).upper()))
		self.Xc_Functional.write_to_file(input_deck)
		self.Vdw_Potential.write_to_file(input_deck)
		input_deck.write("\t\t&END XC \n")


class Xc_Functional:


	def __init__(self, functional):
		self.functional = functional
		self.xc_functional = "no_shortcut"

		if self.functional == "pbe":
			self.Functional = Pbe()

	def asdict(self):
		return {"&XC_FUNCTIONAL" : self.xc_functional, 
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Xc_Functional Class")
		print("--------------------------------------------")
		print("XC Functional       : {}".format(self.xc_functional))
		print("XC Functional Arg.  : {}".format(self.functional))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t\t{} {}\n".format(key,str(value).upper()))
		self.Functional.write_to_file(input_deck)
		input_deck.write("\t\t\t\t&END XC_FUNCTIONAL \n")

class Pbe:
	
	section_parameters = ".True."
	parameterization = "ORIG" # ORIG, PBESOL, REVPBE 
	scale_c = 1E0 
	scale_x = 1E0

	def asdict(self):
		return {"&PBE" : self.section_parameters,
			"PARAMETRIZATION" : self.parameterization,
			"SCALE_C" : self.scale_c,
			"SCALE_X" : self.scale_x, 
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Xc_Functional Class")
		print("--------------------------------------------")
		print("Section Parameters  : {}".format(self.section_parameters))
		print("Parameterization    : {}".format(self.parameterization))
		print("Scale C             : {}".format(self.scale_c))
		print("Scale X             : {}".format(self.scale_x))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t\t\t\t&END PBE \n")

class Vdw_Potential:

	def __init__(self, vdw_potential, functional):
		self.potential = vdw_potential
		self.functional = functional
		if self.potential == "pair_potential":
			self.Potential = Pair_Potential(self.functional)
		if self.potential == "none":
			return

	def asdict(self):
		return {"POTENTIAL_TYPE" : self.potential,
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Vdw_Potential Class")
		print("--------------------------------------------")
		print("Potential Type          : {}".format(self.potential_type))
		print("--------------------------------------------")
	
	def write_to_file(self, input_deck):
		input_deck.write("\t\t\t&VDW_POTENTIAL\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t\t{} {}\n".format(key,str(value).upper()))
		self.Potential.write_to_file(input_deck)
		input_deck.write("\t\t\t&END VDW_POTENTIAL \n")
	
class Pair_Potential:
	
	r_cutoff = 1E1 
	type = "dftd3"
	parameter_file_name = "dftd3.dat"
	calculate_c9_term = ".True."
	reference_c9_term = ".True."
	
	def __init__(self,functional):
		self.reference_functional = functional

	def asdict(self):
		return {"R_CUTOFF" : self.r_cutoff,
			"TYPE" : self.type,
			"PARAMETER_FILE_NAME" : self.parameter_file_name,
			"REFERENCE_FUNCTIONAL" : self.reference_functional,
			"CALCULATE_C9_TERM" : self.calculate_c9_term,
			"REFERENCE_C9_TERM" : self.reference_c9_term, 
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Pair_Potential Class")
		print("--------------------------------------------")
		print("R Cutoff             : {}".format(self.r_cutoff))
		print("Type                 : {}".format(self.type))
		print("Paramter File Name   : {}".format(self.parameter_file_name))
		print("Reference Functional : {}".format(self.reference_functional))
		print("Calculate C9 Term    : {}".format(self.calculate_c9_term))
		print("Reference C9 Term    : {}".format(self.reference_c9_term))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t\t\t\t&PAIR_POTENTIAL\n")
		for key, value in self.asdict().items():
			if key == "PARAMETER_FILE_NAME":
				input_deck.write("\t\t\t\t\t{} {}\n".format(key,value))
			else:
				input_deck.write("\t\t\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t\t\t&END PAIR_POTENTIAL \n")

class Sub_Sys:

	def __init__(self, structure, basis_set, potential):
		self.structure = structure
		self.basis_set = basis_set
		self.potentials = potential 
		self.Cell = Cell(self.structure)
		self.Coord = Coord(self.structure)
		self.kinds = []
		for atom in self.structure.composition.elements:
			self.kinds.append(Kind(atom))

	def write_to_file(self, input_deck):
		input_deck.write("\t&SUBSYS\n")
		self.Cell.write_to_file(input_deck)
		self.Coord.write_to_file(input_deck)
		for Atom in self.kinds:
			Atom.write_to_file(input_deck)
#		self.Kind.write_to_file(input_deck)
		input_deck.write("\t&END SUBSYS\n")

class Cell:
	
	periodic = "XYZ"
	multiple_unit_cell = "1 1 1"
	symmetry = "none"

	def __init__(self, structure):
		self.structure = structure 
		self.abc = "{0:.6f} {1:.6f} {2:.6f}".format(self.structure.lattice.abc[0], self.structure.lattice.abc[1], self.structure.lattice.abc[2])
		self.alpha_beta_gamma = "{0:.6f} {1:.6f} {2:.6f}".format(structure.lattice.angles[0], structure.lattice.angles[1], structure.lattice.angles[2])

	def asdict(self):
		return {"ABC" : self.abc,
			"ALPHA_BETA_GAMMA " : self.alpha_beta_gamma,
			"PERIODIC" : self.periodic,
			"SYMMETRY" : self.symmetry,
			"MULTIPLE_UNIT_CELL" : self.multiple_unit_cell,
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Cell Class")
		print("--------------------------------------------")
		print("ABC                  : {}".format(self.abc))
		print("Alpha Beta Gamma     : {}".format(self.alpha_beta_gamma))
		print("Periodic             : {}".format(self.periodic))
		print("Symmetry             : {}".format(self.symmetry))
		print("Multiple Unit Cell   : {}".format(self.multiple_unit_cell))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t\t\t&CELL\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t\t&END CELL \n")


class Coord:

	def __init__(self, structure):
		self.structure = structure

	def print_options(self):
		print("--------------------------------------------")
		print("Coord Class")
		print("--------------------------------------------")
		for atom in self.structure.sites:
			print("{0:s} {1:.6f} {2:.6f} {3:.6f}".format(atom.specie, atom.coords[0], atom.coords[1], atom.coords[2]))	
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t\t\t&COORD\n")
		for atom in self.structure.sites:
			input_deck.write("\t\t\t\t {0:s} {1:.6f} {2:.6f} {3:.6f} \n".format(atom.specie, atom.coords[0], atom.coords[1], atom.coords[2]))
		input_deck.write("\t\t\t&END COORD \n")

class Kind:
	
	# PASS FUNCTIONAL AND POTENTIAL THROUGHT (EDIT)#

	def __init__(self, element):
		self.element = element
		self.name = element.symbol
		self.basis_set = "TZVP-GTH"
		self.valence = self.valence_count()
		self.potential = "GTH-PBE-q{}".format(self.valence)

	def valence_count(self):
		valence = 0
		for shell in self.element.full_electronic_structure:
			if shell[0] == self.element.row:
				valence += shell[2]
		return valence

	def asdict(self):
		return {"BASIS_SET" : self.basis_set,
			"POTENTIAL" : self.potential,
			} 
	def write_to_file(self, input_deck):
		input_deck.write("\t\t\t&KIND {}\n".format(self.name))
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t\t{} {}\n".format(key,value))
		input_deck.write("\t\t\t&END KIND \n")
