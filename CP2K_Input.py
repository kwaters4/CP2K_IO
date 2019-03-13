import CP2KGeometry as cp2kio


class CP2K_Input_Deck:

	def __init__(self, structure, directory, name = "default", method = "QS", functional = "pbe", potential = "pair_potential"):

		self.name = name
		self.directory = directory
		self.structure = structure
		self.Global = Global(name)
		self.Motion = Motion(self.Global.run_type)
		self.Force_Eval = Force_Eval(method, functional, potential)

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
#		cp2kio.cp2k_structure(structure, filename)
		return

class Global:

	blacs_grid = "SQUARE" #COLUMN, ROW, SQUARE
	fft_library = "FFTW3"
	fft_lengths = ".True."
	print_level = "MEDIUM"
	run_type = "CELL_OPT" # "GEO_OPT" # "ENERGY_FORCE"

	def __init__(self, name):
		self.name = name
	
	def print_options(self):
		print("--------------------------------------------")
		print("Global Class")
		print("--------------------------------------------")
		print("Name        : {}".format(self.name))
		print("Run Type    : {}".format(self.run_type))
		print("Print Level : {}".format(self.print_level))
		print("FFT Library : {}".format(self.fft_library))
		print("FFT Lengths : {}".format(self.fft_lengths))
		print("Blacs grid  : {}".format(self.blacs_grid))
		print("--------------------------------------------")

	def asdict(self):
		return {"PROJECT_NAME" : self.name, 
			"RUN_TYPE" : self.run_type,
			"PRINT_LEVEL" : self.print_level,
			"BLACS_GRID" : self.blacs_grid,
			"FFT_LIBRARY" : self.fft_library,
			"FTT_LENGTHS" : self.fft_lengths,
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
	keep_symmerty = ".False."
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
			"KEEP_SYMMERTY" : self.keep_symmerty,
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
		print("Keep Symmetry  : {}".format(self.keep_symmerty))
		print("Start Value    : {}".format(self.step_start_val))
		print("Max DR.        : {}".format(self.max_dr))
		print("Max Iter.      : {}".format(self.max_iter))
		print("Max Force      : {}".format(self.max_force))
		print("Optimizer      : {}".format(self.optimizer))
		print("Pressure Tol.  : {}".format(self.pressure_tol))
		print("RMS DR         : {}".format(self.rms_dr))
		print("RMS Force      : {}".format(self.rms_force))
		print("Type           : {}".format(self.type))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t &CELL_OPT\n")
		for key, value in self.asdict().items():
			input_deck.write("\t \t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t &END CELL_OPT\n")

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
		input_deck.write("\t &GEO_OPT\n")
		for key, value in self.asdict().items():
			input_deck.write("\t \t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t &END GEO_OPT\n")
	
class Force_Eval:

	stress_tensor = "Analytical" # ANALYTICAL, DIAGONAL_ANALYTICAL, DIAGONAL NUMERICAL, NONE, NUMERICAL 

	def __init__(self, method, functional, potential):
#		self.method = method
		self.Dft = Dft(method, functional, potential)
#		self.Sub_Sys = Sub_Sys()

	def write_to_file(self, input_deck):
		input_deck.write("&FORCE_EVAL\n")
		self.Dft.write_to_file(input_deck)
		input_deck.write("&END FORCE_EVAL\n")

#class Sub_Sys(self):

class Dft:

	basis_set_file_name = "gth_basis_sets"
	potential_file_name = "potential"
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

	def __init__(self, method, functional, potential):
		self.method = method
		self.Scf = Scf()
		if method == "QS":
			self.Qs = Qs() 
		self.Mgrid = Mgrid()
		self.Xc = Xc(functional, potential)

	def asdict(self):
		return {"BASIS_SET_FILE_NAME" : self.basis_set_file_name, 
			"POTENIAL_SET_FILE_NAME" : self.potential_file_name,
			"CHARGE" : self.charge,
			"EXCITATIONS" : self.excitations,
			"MULTIPLICITY" : self.multiplicity,
			"PLUS_U_METHOD" : self.plus_u_method,
			"RELAX_MULTIPLICITY" : self.relax_multiplicity,
			"ROKS" : self.roks,
			"SUBCELLS" : self.subcells,
			"SURE_DIPOLE_CORRECTION" : self.surface_dipole_correction,
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
		input_deck.write("\t &DFT\n")
		for key, value in self.asdict().items():
			input_deck.write("\t \t{} {}\n".format(key,str(value).upper()))
		self.Scf.write_to_file(input_deck)
		if self.method == "QS":
			self.Qs.write_to_file(input_deck)
		self.Mgrid.write_to_file(input_deck)
		input_deck.write("\t &END DFT\n")
	
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
		input_deck.write("\t\t &SCF\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t{} {}\n".format(key,str(value).upper()))
		self.Mixing.write_to_file(input_deck)
		input_deck.write("\t\t &END SCF\n")
	
		return
class Mixing:

	alpha = 2E-1 

	def asdict(self):
		return {"ALPHA" : self.alpha, 
			} 

	def print_options(self):
		print("--------------------------------------------")
		print("Mixing Class")
		print("--------------------------------------------")
		print("Alpha                : {}".format(self.alpha))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t\t\t &MIXING\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t\t &END MIXING \n")

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
		input_deck.write("\t\t &QS\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t &END QS \n")

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
		input_deck.write("\t\t &MGRID\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t &END MGRID \n")

class Xc:

	density_cutoff = 1E-10
	gradient_cutoff = 1E-10
	tau_cutoff = 1E-10

	def __init__(self, functional, potential):
		self.Vdw_Potential = Vdw_Potential(potential)
		self.Functional = Functional(functional)
	
	def print_options(self):
		print("--------------------------------------------")
		print("Xc Class")
		print("--------------------------------------------")
		print("Density Cutoff          : {}".format(self.density_cutoff))
		print("Gradient Cutoff         : {}".format(self.gradient_cutoff))
		print("Tau Cutoff              : {}".format(self.tau_cutoff))
		print("--------------------------------------------")

	def write_to_file(self, input_deck):
		input_deck.write("\t\t &XC\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t &END XC \n")

	
#	section_parameters = "no_shortcut" 
#	functional = "PBE"
#
#class Pbe:

class Functional:

	def __init__(self, functional):
		func = "PBE"

class Vdw_Potential:

	def __init__(self, potential):
		if potential == "pair_potential":
			self.Pair_Poetential = Pair_Potential()
		if potential == "none":
			return

	def print_options(self):
		print("--------------------------------------------")
		print("Vdw_Potential Class")
		print("--------------------------------------------")
		print("Potential Type          : {}".format(self.potential_type))
		print("--------------------------------------------")
	
	def write_to_file(self, input_deck):
		input_deck.write("\t\t\t &VDW_POTENTIAL\n")
		for key, value in self.asdict().items():
			input_deck.write("\t\t\t\t{} {}\n".format(key,str(value).upper()))
		input_deck.write("\t\t\t &END VDW_POTENTIAL \n")
	
class Pair_Potential:
	
	test = 1

