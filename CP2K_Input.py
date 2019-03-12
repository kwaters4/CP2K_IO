import CP2KGeometry as cp2kio


class CP2K_Input_Deck:

	functional = "PBE"
	basis_set = "GTH"
	cell_opt_flag = True
	cell_angles = True
	run_type = "CELL_OPT" # "GEO_OPT" # "ENERGY_FORCE"

	def __init__(self, structure, name = "default"):
		self.Overview = Global(name)
		self.Motion = Motion()
		self.Force_Eval = Force_Eval()
		self.structure = structure

	def __repr__(self):
		return

	def write_file(self):
		filename = "{}.inp".format(self.name)
		# GLOBAL
		# COORDINATES
#		cp2kio.cp2k_structure(structure, filename)
		return

class Global:

	print_level = "MEDIUM"

	def __init__(self, name):
		self.name = name
		return
class Motion:
	if cell_opt_flag == True
		self.Cell_Opt = Cell_Opt()
	else:
		return
	return

class Cell_Opt:
	return

class Geo_Opt:
	return
	
	
class Force_Eval:
	return 


class Dft:
		return

class Scf:
		return

