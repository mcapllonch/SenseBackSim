"""
This module was created by Miguel Capllonch Juan on 13 November 2018.
It handles most aspects of the system's biophysics
"""

import os
import csv
import numpy as np
from neuron import h
from datetime import datetime
import planar as pl

import algebra as alg
import workspace as ws
import anatomy
import electrodes
import analysis as als
import tools as tls



class AxonGeneralModel():
	"""
	This is a class to put in all the possible properties of an 
	axon model type, including anatomical and biophysical general properties
	This class is only instantiated once per model type, so each 
	instance of a physical axon does not have to do everything 
	from scratch
	"""
	def __init__(self, model):
		self.model = model
		self.properties_ = ws.axonmodel_settings[model]

	def build_regressions(self):
		"""
		Build the necessary regression functions of the axon model from the 
		json files in the work space
		"""
		if 'regressions' in self.properties_.keys():
			regfile = os.path.join(ws.folders['models'], 
				self.properties_['regressions'])
			self.regressions = alg.parameter_regressions(regfile)


class Cable:
	"""
	This is the basic class for any cable, be this an axon or a NAELC
	"""
	def __init__(self):

		# Counter. Add one to the class instances counter and save it 
		# as an id number as well
		ws.counters['cables'] += 1
		self.cable_number = ws.counters['cables'] - 1

		# Flag indicating if this cable has recordings
		self.has_recordings = False

		# Basic geometrical properties
		if ws.settings["nerve model"] == "resistor network":
			self.x = ws.nvt.x[self.cable_number]
			self.y = ws.nvt.y[self.cable_number]
			self.r = ws.nvt.r[self.cable_number]
		elif ws.settings["nerve model"] == "simple":
			self.x = ws.nvt.x[self.cable_number + ws.nvt.nNAELC]
			self.y = ws.nvt.y[self.cable_number + ws.nvt.nNAELC]
			self.r = ws.nvt.r[self.cable_number + ws.nvt.nNAELC]

		# Create properties dictionary
		self.properties = {}

		# Find the fascicle to which this cable belongs
		self.properties["fascicle"] = ws.nvt.cables_fascicles[self.cable_number]
		
	def finish_build(self):
		"""
		Perform the last steps in building this cable model
		"""

		# Append some properties to the corresponding dictionaries
		i = self.cable_number

		# Append the cell to the list of cables in the work space
		ws.cables_hoc.append(self.cable_cell)

		# Sections and lengths
		secs = list(self.cable_cell.all)
		ws.sections[i] = secs[:] 
		ws.sectypes[i] = self.properties['sctyps']
		ws.seclens[i] = [sec.L for sec in secs]
		# print('created seclens:', i)
		
		# Segments
		for j, sec in enumerate(secs):
			ws.secsegs[(i, j)] = list(sec.allseg())[1:-1]
		segs = tls.flattenlist([list(sec.allseg())[1:-1] for sec in ws.sections[i]])
		ws.segments[i] = segs[:]
		
		# Longitudinal profile (z-profiles)
		ws.zprofs[i] = anatomy.zprofile(self.cable_cell)
		
		# Maximum number of segments
		if len(ws.zprofs[i]) > ws.maxnseg:
			ws.maxnseg = len(ws.zprofs[i])

		# Store attributes
		self.properties['zprofile'] = ws.zprofs[i]

		# Intracellular area
		try:
			self.intracellular_area = np.pi * (0.5 * self.properties["axonD"]) ** 2
		except KeyError:
			# This is a NAELC or is an axon not configured to have axonD (which will be tackled if it happens)
			self.intracellular_area = 0.

	def set_recordings(self):
		"""
		Prepare the recording vectors for the desired variables
		"""
		# ws.log("Setting recordings for cable %i"%self.cable_number)

		# Update the flag indicating if this cable has recordings
		self.has_recordings = True

		# Variables to record from the settings
		rvars = ws.settings["data records"]["variables"]

		# Dictionary with all the recordings:
		# One entry per variable. Each entry contains a list of 
		# recording identifiers (one for each variable and segment of  
		# the cable)
		self.recordings = {}

		# Which sections to record
		where_to_record = ws.settings["data records"]["record sections"]

		# If this cable is a NAELC, record only vext[0]
		if self.cabletype == 'NAELC':
			rvars = ['vext[0]']
			where_to_record = ["all"]

		# Record the variables
		timings = []
		for var in rvars:

			# Create the list
			self.recordings[var] = []

			# Create the recording vectors
			for i, seg in enumerate(ws.segments[self.cable_number]):

				record = False
				
				t0 = datetime.now()

				if where_to_record[0] == "all":
					record = True
				else:
					# for sectypename in where_to_record:
					# 	if self.properties["sctyps"][i] == sectypename:
							# record = True
					if self.properties["sctyps"][i] in where_to_record:
						record = True

				if record:
					self.recordings[var].append({
						"number": ws.counters['recordings'], 
						"segment": str(seg)
						})
					als.set_recording(seg, var)

				t1 = datetime.now()
				timings.append((t1 - t0).total_seconds())

		ws.log("\tAverage time needed by als.set_recording for cable %i: %f ms. Max: %f ms"%(self.cable_number, 1e3 * np.array(timings).mean(), 1e3 * np.array(timings).max()))


class Axon(Cable):
	"""
	This is a class to put in all the possible properties of an 
	axon model, including anatomy and biophysics
	This class inherits its properties anatomical and biophysical 
	general properties from THE INSTANCE of AxonGeneralModel which 
	describes the model type
	"""
	def __init__(self, model):

		# Inheritance
		Cable.__init__(self)
		# Type of cable
		self.cabletype = 'Axon'
		# print('***Creating axon with the model %s'%model)
		# Counter
		ws.counters['axons'] += 1
		# Specific identifier; specific for this type of cable
		self.specific_number = ws.counters['axons'] - 1
		# Properties
		self.__dict__.update(ws.ax_gral_models[model].__dict__.copy())
		# Put AxonGeneralModel.properties_ into self.properties
		self.properties.update(self.properties_)
		self.properties['fiber_radius'] = self.r
		self.properties['fiberD'] = 2. * self.properties['fiber_radius']
		# Axon model
		# self.properties['model'] = ws.nvt.models[self.cable_number]

	def build(self):
		"""
		General class that chooses between the build methods for 
		myelinated or unmyelinated axons
		"""
		build_methods = {
			'myelinated': self.build_myelinated, 
			'unmyelinated': self.build_unmyelinated
		}
		build_methods[self.properties['type']]()

		# Finish building this cable and saving its properties
		self.finish_build()


	def build_unmyelinated(self):
		"""
		Build the axon as unmyelinated
		"""
		nseg = anatomy.nseg_dlambda_rule(
			ws.length, 
			self.properties['fiber_radius'], 
			self.properties['rhoa'], 
			self.properties['cm']
			)
		self.properties['nseg'] = nseg
		self.cable_cell = h.UnMyelAxon(self.properties)
		self.properties['sctyps'] = np.zeros(nseg, dtype=np.int)

	def build_myelinated(self):
		"""
		Build the necessary properties of the axon model from the 
		json files in the work space
		"""

		# Variables that depend on a regression from some stored data
		iv = 'independent variable'
		rgs = self.regressions
		for k, func in rgs.items():
			if k != iv:
				x = self.properties[rgs[iv]]
				self.properties[k] = func(x)

		# Other properties:

		# Variable name shortenings
		deltax = self.properties['deltax']

		# 1. STIN_length
		# This is still model-specific hard-coding
		if self.model == 'MRG' or self.model == 'gaines_sensory':
			nodl = self.properties['nodelength']
			pl1 = self.properties['paralength1']
			pl2 = self.properties['paralength2']
			self.properties['STIN_length'] = deltax - nodl - 2. * (pl1 + pl2)
		elif self.model == 'MRG_STINonly':
			nodl = self.properties['nodelength']
			self.properties['STIN_length'] = deltax - nodl

		# Start point along the z-axis
		if ws.settings["nerve model"] == "resistor network":
			self.properties["start"] = deltax * ws.nvt.start_positions[self.cable_number]
		elif ws.settings["nerve model"] == "simple":
			self.properties["start"] = deltax * ws.nvt.start_positions[self.cable_number + ws.nvt.nNAELC]

		# 3. Total internodal length
		self.properties['in_l'] = deltax - self.properties['nodelength']

		# 4. Update the minimum length of one chunk of sections
		# and of the node to fiber ratio, g
		ws.deltax_min = min((ws.deltax_min, deltax))
		ws.g_min = min((ws.g_min, self.properties['g']))

		# Build the axon in hoc
		sctyps, section_counter, n_max, lengths = anatomy.all_vars(self)

		# Lengths to hoc
		lengths_hoc = {}
		for key, value in lengths.items():
			lengths_hoc[key] = h.Vector(value)

		# Number of sections of each type and section types
		nsecs = sum(section_counter.values())
		sectypes = h.List()
		for st in sctyps:
			sectypes.append(st)

		# Store variables as attributes inside the 'properties' dictionary
		self.properties['nsecs'] = nsecs
		self.properties['sectypes'] = sectypes
		self.properties['section_counter'] = section_counter
		self.properties['mx_n'] = n_max
		self.properties['lengths'] = lengths_hoc

		# Axon in hoc
		if self.model == 'MRG':
			# print(' ----------- Setting up %s'%self.model)
			self.cable_cell = h.MRGMyelAxon(self.properties)
		elif self.model == 'gaines_sensory':
			# print(' ----------- Setting up %s'%self.model)
			self.cable_cell = h.Gaines2018Sensory(self.properties)
		
		# z-profile
		zp = anatomy.zprofile(self.cable_cell)

		# If this is the smallest axon, save its longitudinal profile
		if self.properties['fiber_radius'] == ws.raxmin:
			ws.zlayers_min = zp

		# Number of segments on each section of this axon
		this_nseg = [sec.nseg for sec in self.cable_cell.all]

		# I need sctyps to reflect all the segments, not just the sections
		sctyps_str = [this_nseg[ii] * [sctyps[ii].sectype] for ii in range(len(sctyps))]
		sctyps_str = np.array(tls.flattenlist(sctyps_str))

		# Get the z-profiles of the nodes of Ranvier only
		zprofs_RN = zp[np.where(sctyps_str == 'node')[0]]
		ws.zprofs_RN[self.cable_number] = zprofs_RN

		# Store these properties
		self.properties['zprofs_RN'] = zprofs_RN
		self.properties['sctyps'] = sctyps_str


class NAELC(Cable):
	""" 
	This is a Non-Axonal Cable 
	"""
	def __init__(self, *dummyargs):
		# Inheritance
		Cable.__init__(self)
		# Type of cable
		self.cabletype = 'NAELC'
		# Counter
		ws.counters['NAELC'] += 1
		# Specific identifier; specific for this type of cable
		self.specific_number = ws.counters['NAELC'] - 1
		# Properties
		self.properties = self.properties.copy()
		self.properties.update({
					'diam': ws.cdiam, 
					'length': ws.length, 
					'extlayers': ws.NAELC_extlayers
				})

	def build(self):
		"""
		Build the cable into the hoc space
		"""
		self.properties['nseg'] = ws.NAELC_nseg
		self.properties['sctyps'] = np.zeros(ws.NAELC_nseg, dtype=np.int)
		self.cable_cell = h.Wire(self.properties)

		# Finish building this cable and saving its properties
		self.finish_build()

class ResistorNetwork():
	"""
	This is a class for the Resistor Network that describes the 
	interactions between all axons
	"""
	def __init__(self):
		pass
		
	def get_resistances(self):
		""" Calculate the values of the transverse resistances """

		self.resistors = {}
		self.res_lengs = {}
		for pair in ws.pairs:
			i, j = pair
			# Calculation of the resistor's value: 
			# shape factor of the resistor (adimensional)
			# and value in MOhm*cm
			l = ws.nvt.len_con[pair]
			if ws.EC["resistor network"]["interaxonal connections"] == "membranes":
				l -= (ws.nvt.r[i] + ws.nvt.r[j])
			elif ws.EC["resistor network"]["interaxonal connections"] == "centers":
				pass
			
			# Width
			w = ws.nvt.len_seg[pair]
			# Shape factor (adimensional (um / um))
			shape_factor = l / w
			# Large shape factors mean large resistances; not interested
			if shape_factor < 100:

				# Resistance value
				res = ws.rho['endo_trans'] * shape_factor

				# Perineurium
				resist_perineurium = 0.

				# If the pair is at opposite sides of the perineurium, 
				# add the perineurium's resistance p.u.l or 
				# resistivity (Ohm*cm)
				if ('NAELC' in ws.nvt.cables[i] and 'xon' in ws.nvt.cables[j]) \
				or ('NAELC' in ws.nvt.cables[j] and 'xon' in ws.nvt.cables[i]):
					resist_perineurium += ws.resist_perineurium

				# Also, that happens when two axons are in connected but 
				# in different fascicles
				# if 'xon' in ws.nvt.cables[i] and 'xon' in ws.nvt.cables[j]:
				if ws.cables[i].cabletype == "Axon" and ws.cables[j].cabletype == "Axon":
					# ws.log("Checking if there\'s double layer of perineurium for axons %i and %i:"%(ws.cables[i].specific_number, ws.cables[j].specific_number))
					# ws.log("%s, %s"%(ws.cables[i].properties["fascicle"], ws.cables[j].properties["fascicle"]))
					if ws.cables[i].properties["fascicle"] != ws.cables[j].properties["fascicle"]:
						resist_perineurium += 2. * ws.resist_perineurium
						# ws.log("\tAdding double perineurium between axons %i and %i"%(ws.cables[i].specific_number, ws.cables[j].specific_number))
						# ws.log("\tThese axons are placed at (%f, %f) um and (%f, %f) um, respectively"%(ws.cables[i].x, ws.cables[i].y, ws.cables[j].x, ws.cables[j].y))

				# Finally, add it up
				res += resist_perineurium / w
				# Store in the dictionary
				self.resistors[pair] = res

	def connect_resistors(self):
		""" Connect the transverse resistors with ephapses in NEURON """

		def locate_and_create(pair, positions):
			""" For a pair of cables and a set of resistor positions, 
			locate, create and connect those resistors to the cables """

			# Pair indices
			i, j = pair

			# Get resistors' lengths
			rl = get_rl(pair, positions)

			# Iterate over resistors positions
			for zr, l in zip(positions, rl):

				# For each cable, see where this falls

				# On cable i:
				i_sec, z_on_i = anatomy.locate(ws.seclens[i], zr)
				seci = ws.sections[i][i_sec]
				sgi = seci(z_on_i)

				# On cable j:
				j_sec, z_on_j = anatomy.locate(ws.seclens[j], zr)
				secj = ws.sections[j][j_sec]
				sgj = secj(z_on_j)

				# Value of the resistor (MOhm)
				ephres = res / l
				# Create the resistor
				create_eph(ephres, sgi, sgj, seci, secj, i, j, 
					zr)
				# ws.log("Created a transverse resistor at z = %f um with value %f MOhm"%(zr, ephres))

		ws.eph_counter = 0
		ws.nrn_res = {}
		for pair, res in self.resistors.items():
			i, j = pair

			# Cable
			cable_cell_i = ws.cables_hoc[i]
			cable_cell_j = ws.cables_hoc[j]

			# I have to connect the cables in a staggered manner

			# If I am only taking the nodes of Ranvier into account
			# and the connection is to be set between two axons:
			if ws.EC["resistor network"]["transverse resistor locations"] == "regular":

				# Separation between transverse resistors
				sep = float(ws.EC["resistor network"]["resistors separation"])

				# Positions of the resistors
				# res_positions = np.arange(0, ws.length, sep)
				res_positions = np.arange(0.5, ws.length, sep)

				# Create and connect the resistors
				locate_and_create(pair, res_positions)

			elif ws.EC["resistor network"]["transverse resistor locations"] == "only nodes of Ranvier" and \
			('xon' in ws.cables[i].cabletype and 'xon' in ws.cables[j].cabletype):

				# Get resistors' positions and lengths (um)

				# Merge the ws.zprofs_RN of the two axons
				# This array will contain the postions of the resistors 
				# over the z-axis
				zp_pair = np.append(ws.zprofs_RN[i], ws.zprofs_RN[j])
				zp_pair.sort()
				# Remove repeated values (can --will, actually, happen if 
				# there's aligned nodes)
				zp_pair = np.unique(zp_pair)

				# Create and connect the resistors
				locate_and_create(pair, zp_pair)

			else:
				# If I am taking all nodes into account, the compartments of only one
				# of the two cables (arbitrarily chosen; cable 'i' by default) 
				# are plugged to transverse resistors. These are then connected 
				# to the second cable wherever they must, regardless of the 
				# topology of this second cable.

				# Get resistors' lengths (um)
				# From what is explained above, in this case the z-profile of one
				# of the cables only is used as zp_pair

				rl = get_rl(pair, ws.zprofs[i])
				
				# Iterate over segments
				for ii, sgi in enumerate(ws.segments[i]):

					# Get the section and segment of cable j which sgi lies against
					zi = ws.zprofs[i][ii]
					j_sec, z_ioverj = anatomy.locate(ws.seclens[j], zi)
					seci = sgi.sec
					secj = ws.sections[j][j_sec]
					sgj = secj(z_ioverj)

					# Value of the resistor (MOhm)
					ephres = res / rl[ii]
					# Create the resistor
					create_eph(ephres, sgi, sgj, seci, secj, i, j, 
						zi)

########################################################################
# Functions

def finish_fascicles():
	"""
	Finish giving basic information to the fascicles once the cables have been set up and the resistivities have been adjusted (if they have)
	"""
	for k, fas in ws.nvt.fascicles.items():
		# Intracellular area
		ws.nvt.fascicles[k].get_intracellular_area()
		# Print properties
		fas.print_properties()
	# Now print the properties of the nerve
	ws.nvt.print_properties()

def get_rl(pair, zp_pair):
	"""
	Get the lengths (all over the z-axis) of all the transverse
	resistors for a pair of cables
	They will be used to obtain the values of each transverse resistor in MOhm.
	"""

	# Distances between nodes
	dz = zp_pair[1:] - zp_pair[:-1]
	# Get the lengths of the resistors (um)
	rl = 0.5 * (dz[1:] + dz[:-1])
	# Add one element on the left: first resistor
	# and one on the right: last resistor
	rl = np.append(np.array([zp_pair[0] + 0.5 * dz[0]]), rl)
	rl = np.append(rl, np.array([0.5 * dz[-1] + ws.length - zp_pair[-1]]))
	# Store its vale
	ws.rnet.res_lengs[pair] = rl
	return rl

# def create_eph(ephres, sgi, sgj, seci, secj, i, j, nrn_res, zr):
def create_eph(ephres, sgi, sgj, seci, secj, i, j, zr):
	"""
	Create a resistor that acts as an ephaptic connection
	"""

	# Choose the layer depending on the type of fiber or cable
	# Layer where to point to
	layeri = ws.cables[i].properties['extlayers']
	layerj = ws.cables[j].properties['extlayers']

	# If zr is not in a cuff, we will make the ephaptic connection
	# super resistive and hence disappear
	if not ws.EC["resistor network"]["outside cuffs"]:
		ephres = disconnect_under_cuff(zr, ephres)

	# If the connection is too large, do not implement it
	if ephres < 1.e2:

		eph = h.EphapDBA(2)
		eph.ex0_loc(0, sgi.x, layeri, sec=seci)
		eph.ex1_loc(1, sgj.x, layerj, sec=secj)
		eph.add_submatrix(0, 1, ephres)
		eph.install()
		# Disconnect segments from ground
		# This is redundant, but just in case
		sgi.xg[layeri - 1] = 0.
		sgj.xg[layerj - 1] = 0.

		# Add to dictionary
		ws.nrn_res = tls.append_items(ws.nrn_res, (i, j), [eph])
		ws.eph_counter += 1

def disconnect_under_cuff(zi, ephres):
	"""
	Check if a transverse resistor is not under the region of a cuff 
	electrode and if so, give it a highest value so it's not connected
	"""
	if electrodes.check_cuff_cover(zi) is None:
		ephres *= ws.rho['SuperResistive'] / ws.rho['endo_trans']
	return ephres

def build_cables():
	"""
	Build all the cables, NAELC and fibers
	"""

	ws.log("Building the cables in NEURON")
	t0 = datetime.now()

	# Dictionaries where to store information in the workspace
	ws.zprofs_RN = {}
	ws.sections = {}
	ws.sectypes = {}
	ws.secsegs = {}
	ws.seclens = {}
	ws.zprofs = {}
	ws.segments = {}
	cables = []
	ws.cables_hoc = h.List()

	# New building method
	# First, axon model properties
	# Axon general model(s)
	ws.ax_gral_models = {}
	# ws.ax_gral_models[ws.settings['axon model']] = AxonGeneralModel(ws.settings['axon model'])
	for m in ws.nvt.models_set:
		if 'NAELC' not in m:
			ws.ax_gral_models[m] = AxonGeneralModel(m)
			ws.ax_gral_models[m].build_regressions()

	# Classes of cables
	ws.cable_classes = {
		"NAELC": NAELC, 
		"Axon": Axon
	}

	# First, instantiate all the cables.
	# Then, build axons. From them, obtain the desired nseg for NAELC.
	# Finally, build the NAELC using this value

	# Instantiate everything
	nvt = ws.nvt
	for i in range(nvt.nc):

		# Create the cable(s)

		if ws.settings["nerve model"] == "resistor network":
			# Use NAELC
			# cable = ws.cable_classes[nvt.cables[i]](ws.settings['axon model'])
			cable = ws.cable_classes[nvt.cables[i]](nvt.models[i])
			# Build it only if it is an axon
			# The NAELC will be built later on
			if nvt.cables[i] == 'Axon':
				# print('building', i, nvt.cables[i])
				cable.build()
			# Append it
			cables.append(cable)

		elif ws.settings["nerve model"] == "simple":
			# Don't use NAELC
			if nvt.cables[i] == 'Axon':
				# cable = ws.cable_classes[nvt.cables[i]](ws.settings['axon model'])
				cable = ws.cable_classes[nvt.cables[i]](nvt.models[i])
				# Build the axon
				cable.build()
				# Append it
				cables.append(cable)
	
	# Obtain nseg for the NAELC if there are axons
	if ws.counters['axons'] > 0:
		# There's axons
		daxmin = 2. * ws.raxmin
		if ws.EC["resistor network"]["transverse resistor locations"] == "only nodes of Ranvier":
			nseg = int(ws.length / ws.deltax_min) + 1
		elif ws.EC["resistor network"]["transverse resistor locations"] == "regular":
			nseg = int(ws.length / float(ws.EC["resistor network"]["resistors separation"])) + 1
		else:
			# nseg for an internode
			nseg = anatomy.nseg_dlambda_rule(ws.deltax_min, ws.g_min * daxmin, 
				cables[ws.counters['NAELC']].properties['rhoa'], cm=cables[ws.counters['NAELC']].properties['cm'])
			# nseg for the whole nerve length
			nseg = int(ws.length / (ws.deltax_min / nseg)) + 1

		# Make sure it's an odd number
		# This adds 2 if it's odd, but it doesn't matter; it's even better
		nseg += nseg % 2 + 1

	else:
		# If we're here, there are no axons
		nseg = ws.nseg

	ws.NAELC_nseg = nseg
	ws.log("NAELC segment length = %f"%(ws.length / nseg))
	ws.log("nseg: %i"%nseg)
	if ws.counters['axons'] > 0:
		ws.log("Min. deltax: %f"%ws.deltax_min)

	# Now build the NAELC according to the value of nseg 
	# obtained from the axons
	for i in range(ws.counters['NAELC']):
		# print('building', i, cables[i])
		cables[i].build()

	# print(ws.nvt.cables)

	# Store the cables in the workspace
	ws.cables = cables

	# Maximum number of segments. This is a number to be found
	ws.maxnseg = 0

	t1 = datetime.now()
	ws.log("Created %i cables in %i seconds"%(ws.nvt.nc, (t1 - t0).total_seconds()))

def electrics():
	"""
	In this function, xraxial is set for every axons and all 
	connections to ground are made
	"""

	# The following loop over cables and segments and calculates:
	# - The cross-sectional area of each cable and therefore:
	#  - its value for xraxial
	#  - its connections to ground at the ends of the nerve

	for i, seglist in ws.segments.items():
		
		# Level or layer where to point to
		level = ws.cables[i].properties['extlayers'] - 1

		# Cross-sectional area and xraxial
		# Cross-sectional area assigned to the cable (cm2)
		carea = ws.nvt.free_areas[i] * 1.e-8

		# Attribute xraxial
		for j, segment in enumerate(seglist):
			# First and foremost, disconnect ALL segments from ground
			segment.xg[level] = 0.

			# If there's no EC outside the cuff-covered regions,
			# check if this point is under a cuff and connect it to ground
			if not ws.EC["resistor network"]["outside cuffs"]:
				cuff_cover = electrodes.check_cuff_cover(ws.zprofs[i][j])
				if cuff_cover is None:
					segment.xg[level] = 1.e9

			# Attribute xraxial (MOhm/cm; so I have to convert rho (in MOhm*um) 
			# to MOhm*cm)

			# Longitudinal resistivity for this cable as a weighted average between epineurium and endoneurium (the endoneurium is, in turn, a weighted average between all the fascicles intersecting the cross-sectional area of this cable)
			ctendo = ws.nvt.cables_tissues[i]['endoneurium']
			ctepi = ws.nvt.cables_tissues[i]['epineurium']
			fas_ = ws.nvt.fascicles
			resistivity = ws.rho['endo_longt'] * sum([ctendo[k] for k in fas_]) + ctepi * ws.rho['epi_longt']

			# Apply the resistivity in xraxial
			xrx =  1.e-4 * resistivity / carea
			segment.xraxial[level] = xrx

		# Connect to ground all cables at both ends of the nerve
		if ws.EC["resistor network"]["outside cuffs"]:
			
			segleft = seglist[0]
			segrght = seglist[-1]
			
			# Note: The areas of each segment need to be in cm2
			segleft.xg[level] = ws.xg_distal[i] / (segleft.area() * 1.e-8)
			segrght.xg[level] = ws.xg_proxim[i] / (segrght.area() * 1.e-8)
		else:
			# In this case, just ground the ends
			# This is regardless of the ends being cuffed or not; 
			# these are at the base of the cylinder
			seglist[0].sec(0.).xg[level] = 1.e9
			seglist[-1].sec(1.).xg[level] = 1.e9

	# Connect grounds to the nerve's walls (not the ends; these are already sorted)
	# The xg[level] in NEURON for each segment is in mho/cm2
	# All points in the contour are connected to a far ground
	for i in np.arange(ws.npc):
		# Level or layer where to point to
		level = ws.cables[i].properties['extlayers'] - 1
		for j, segment in enumerate(ws.segments[i]):

			# Convert area to cm2 (it's in um2)
			# Cable's cylindrical wall area
			area = segment.area() * 1.e-8

			# Set xg
			# This depends on whether a segment is under the cuff or not
			cuff_cover = electrodes.check_cuff_cover(ws.zprofs[i][j])
			if cuff_cover is not None:
				# xg must be in mho/cm2, so xg_rwall should be in mho; it is

				# If this segment is at one end of the nerve, it already 
				# has one connection to ground that needs to be added
				if segment in (ws.segments[i][0], ws.segments[i][-1]):
					segment.xg[level] += \
						ws.xg_rwall['with_cuff'][cuff_cover] / area
				else:
					segment.xg[level] = \
						ws.xg_rwall['with_cuff'][cuff_cover] / area
			else:
				if ws.EC["resistor network"]["outside cuffs"]:
					# If this segment is at one end of the nerve, it already 
					# has one connection to ground that needs to be added
					if segment in (ws.segments[i][0], ws.segments[i][-1]):
						segment.xg[level] += \
							ws.xg_rwall['wout_cuff'] / area
					else:
						segment.xg[level] = \
							ws.xg_rwall['wout_cuff'] / area
				else:
					# If this is not meant to be used, just ground it
					# (in a future version, I will simply crop these
					# bits of NAELC outside the cuffs)
					segment.xg[level] = 1.e9

def set_direct_extracellular_stimulation():
	""" Set up the stimulation over e_extracellular on the axons """

	# File to read the data from
	espath = os.path.join(ws.folders['data/load'], 'extstim.csv')

	# Create the time window
	# All ones for now
	# It will be cropped when the information for each pulse is fetched
	time_window = np.ones(ws.nt)

	# Time series of the extracellular field for the axons
	stim_tseries = []

	# Dictionaries with the arrays and more info
	pulses = {}

	# Open and read
	with open(espath, 'r') as f:
		reader = csv.reader(f)

		# Read line by line
		for row in reader:

			# First element of each row
			key = row[0]

			# New pulse
			if key == 'Pulse':
				k = int(row[1])
				pulses[k] = {}
				# Reset the time window
				time_window[:] = 1.

			else:

				# Parse the information in the line
				try:
					# See if this is an axon
					i = int(key)
				except ValueError:
					# If not, see if row[1] is a numerical value
					try:
						pulses[k][key] = float(row[1])
					except ValueError:
						# If not, it's a string; just save it as such 
						pulses[k][key] = row[1]
				else:
					# If it is an axon, save the array
					field = np.array([float(x) for x in row[1:]])

					# If we're here, this means we already have all the information
					# Give values to the time window (only once)
					if i == 0:
						if pulses[k]['Type'] == 'square pulse':
							time_window[np.where(ws.tarray < pulses[k]['Delay'])] = 0.
							time_window[np.where(ws.tarray >= pulses[k]['Delay'] + pulses[k]['Duration'])] = 0.

					# If this is the first pulse, create the time series for 
					# this axon
					if k == 0:
						stim_tseries.append(np.zeros((len(field), len(time_window))))

					# Now create the time series of the extracellular field 
					# for each segment of this axon and play it in NEURON
					for j, seg in enumerate(ws.segments[i]):
						stim_tseries[i][j] += field[j] * time_window
			
	# Now that each axon has its full stimulation time 
	# series, play it (fts: 'field time series')

	ws.stim_field_vectors = []

	for i in list(ws.segments.keys())[ws.counters['NAELC']:]:

		for j, seg in enumerate(ws.segments[i]):

			# Create and play the time series
			fts = h.Vector(stim_tseries[i][j])
			fts.play(seg._ref_e_extracellular, h.dt)

			# Add it to the list
			# This is necessary for this to work
			# Otherwise, fts is lost
			ws.stim_field_vectors.append(fts)
