"""
This module handles virtually all the necessary processes for setting 
up stimulation and recording from electrodes.
One more basic function (analysis.set_recording) that is used in this 
module is defined in the module analysis. This is so because this 
function is not necessarily used by an electrode
"""

import os
import numpy as np
from collections import OrderedDict
from neuron import h

import anatomy
import geometry as geo
import analysis as als
import workspace as ws



class CuffRing():
	"""
	Ring of a cuff electrode

	Input parameters:

	z: Position of the ring along the z-axis
	d: Diameter
	npads: Number of pads
	thts: Angular positions of the pads

	Note about the stimulation set-up: this class has no method 
	set_stimulation because some rings won't have any stimulation at 
	all if it is not defined in the configuration files
	"""
	def __init__(self, ring_number, z, center, d, npads, thts):

		self.ring_number = ring_number
		self.z = z
		self.center = center
		self.d = d
		self.npads = npads
		self.thts = thts
		# Dictionary with the information about all the injections
		self.injections = {}
		# Dictionary with the information about all the recordings
		self.recordings = {}

		# Locate the contact points with the nerve
		self.locate_pads_on_nerve()

	def locate_pads_on_nerve(self):
		""" Find the locations in the nerve corresponding 
		to the pads """

		# Iterate over the ring's pads
		self.contact_points = []

		for i_pad in range(self.npads):
			# The pad's angular position determines which cable cell(s) 
			# is (are) stimulated

			# Select which of the targets is in contact 
			# with a pad 'i_pad' according to its angular position
			angdiff = np.abs(ws.contour_angles - self.thts[i_pad])
			which_cables = np.where(angdiff == angdiff.min())[0]
			
			# Iterate over possible cables
			ncables_here = len(which_cables)
			for i_cable in which_cables:
				# For cable_cell i, see where this is
				try:
					i_sec, zpos = anatomy.locate(ws.seclens[i_cable], 
						self.z)
				except TypeError:
					msg = 'ERROR: A ring of an electrode could not be placed on the nerve'
					ws.log(msg)
					ws.terminate()
				# Add the point
				point = {
						'cable': i_cable, 
						'section': i_sec, 
						'z': zpos
					}
				try:
					self.contact_points[i_pad].append(point)
				except IndexError:
					self.contact_points.append([point])

	def assign_injections(self, pad, protocol):
		""" Assign current injections to the pads """
		
		# Dictionary with the information about all the injections
		self.injections[pad] = []

		# Iterate over the contact points of this pad
		# NO! Just one pad is allowed. Otherwise, the current with 
		# intensity 'amp', which is meant for only one injection point, 
		# would be added to all the self.contact_points, so we would 
		# have (len(self.contact_points) - 1) more current injections 
		# than we intended
		point = self.contact_points[pad][0]

		# IClamp object(s)
		# Only square pulses are supported for now
		if protocol['type'] == 'square pulses':
			square_pulses(self.injections[pad], point, protocol)
		if protocol['type'] == 'ground':
			ground_segment(point)

	def set_recordings(self):
		""" Assign the recording vectors to each pad """

		# Iterate over all pads on this ring
		for pad in range(self.npads):

			# Dictionary with the information about all the recordings
			self.recordings[pad] = {}

			# Only one is allowed for each pad in this case
			point = self.contact_points[pad][0]
			
			# Set the recording for this pad
			set_recording(self.recordings[pad], point, 'vext[0]')


class CuffElectrode():
	"""
	Cuff electrode with as many rings as one wants

	Input parameters:

	z: Position (center) along the z-axis
	nrings: Number of rings
	rng_sp: Separation between rings
	nppr: Number of pads per ring
	l: Length of the cuff
	d_in: Inner diameter
	thickness: Insulator thickness
	rho: Insulator resistivity (Ohm*cm)
	"""
	def __init__(self, name):

		self.name = name
		settings = ws.electrodes_settings[name]
		self.type = settings['type']
		self.role = settings['role']
		self.z = settings['z-position']
		self.nrings = settings['number of rings']
		self.rng_sp = settings['inter-ring distance']
		self.nppr = settings['pads per ring']
		self.length = settings['length']
		self.thickness = settings['thickness']
		self.rho  = settings['insulator resistivity']
		self.center = ws.nvt.circumcenter
		self.d_in = ws.nvt.circumdiameter
		self.d_out = self.d_in + 2. * self.thickness
		if self.role == 'stimulation':
			self.stimulation_protocol = \
				settings['stimulation protocol']

		# Necessary by-products

		# Ends of the cuff
		z = self.z
		l = self.length
		nrings = self.nrings
		left_end = z - 0.5 * l
		rght_end = z + 0.5 * l

		# Positions of the rings
		rings_halfspan = 0.5 * (nrings - 1) * self.rng_sp
		z_ring_left = z - rings_halfspan
		z_ring_rght = z + rings_halfspan
		self.rings_poss = np.linspace(z_ring_left, z_ring_rght, nrings)

		# Angular positions of the pads on the rings
		self.thts = settings['angular positions']
		# Consistency check
		if len(self.thts) != self.nppr:
			msg = 'ERROR: Different number of pads and angles for them'
			ws.log(msg)
			ws.terminate()

		# Warning if size mismatch
		if z_ring_left < left_end or z_ring_rght > rght_end:
			msg = 'ERROR: Cuff electrode rings outside the cuff.'
			ws.log(msg)
			ws.terminate()

		# Store variables into attributes
		self.left_end = left_end
		self.rght_end = rght_end
		self.z_ring_left = z_ring_left
		self.z_ring_rght = z_ring_rght

		# Create rings
		self.create_rings()

	def create_rings(self):
		""" Create the rings """
		self.rings = []
		for i in range(self.nrings):
			z = self.rings_poss[i]
			self.rings.append(CuffRing(i, z, self.center, self.d_in, 
				self.nppr, self.thts))

	def set_stimulation(self):
		""" Set up the current injections according to the assigned 
		protocol """
		sp = self.stimulation_protocol
		for pad, protocol in sp.items():
			ring_number, pad_no = [int(x) for x in pad.split(',')]
			self.rings[ring_number].assign_injections(pad_no, protocol)

	def set_recordings(self):
		""" Set up the recording vectors on the pads """
		for ring in self.rings:
			ring.set_recordings()


class PointElectrode():
	"""
	Extracellular point electrode for intracellular current injections, 
	recordings or connections to ground

	Input parameters:

	z: Position along the z-axis
	"""
	def __init__(self, name, index):

		self.name = name
		self.index = index
		self.settings = ws.electrodes_settings[name]
		self.type = self.settings['type']
		self.role = self.settings['role']
		
	def set_stim_info(self):
		""" Set up just a basic thing about stimulation info """

		# Stimulation protocol
		if self.role == 'stimulation':
			try:
				self.stimulation_protocol = \
					self.settings['stimulation protocol'][self.index]
			except KeyError:
				self.stimulation_protocol = \
					self.settings['stimulation protocol']
		if self.role == 'ground':
			self.stimulation_protocol = {
				'type': 'ground'
			}

	def set_stimulation(self):
		""" Set up the current injections according to the assigned 
		protocol """

		# Dictionary containing the information about all 
		# the current injections
		self.injections = {}

		# Injections list
		# There's only one 'pad' for an intracellular electrode
		self.injections[0] = []

		# Fetch protocol
		protocol = self.stimulation_protocol

		# Create IClamp object(s)
		# Only square pulses are supported for now
		if protocol['type'] == 'square pulses':
			ws.log('Setting up injection on %s'%str(self.point))
			square_pulses(self.injections[0], self.point, protocol)
		if protocol['type'] == 'ground':
			ws.log('Setting up ground on %s'%str(self.point))
			ground_segment(self.point)

	def set_recordings(self):
		""" Set up the recording vector """
		self.recordings = {}
		self.recordings[0] = {}
		set_recording(self.recordings[0], self.point, 'v')


class IntracellularElectrode(PointElectrode):
	"""
	Intracellular electrode for intracellular current injections or 
	recordings

	Input parameters:

	z: Position along the z-axis
	"""
	def __init__(self, name, index=0):

		PointElectrode.__init__(self, name, 0)
		settings = self.settings

		try:
			self.z = settings['z-position']
		except KeyError:
			self.axon = settings['axon']
			try:
				self.node = settings['node of Ranvier']
			except KeyError:
				self.node = settings['internode']
				self.position = settings['position']
			else:
				# For a node of Ranvier, inject current in the middle
				self.position = 0.5
		else:
			# Warning if position mismatch
			if self.z < 0. or self.z > ws.length:
				ws.log('ERROR: Electrode is outside the limits of the nerve.')
				ws.terminate()
		# Note: the warnings or errors for targetting a non-existing 
		# node will be automatically raised by NEURON

		# Setting up this point is necessary to set up stimulation and 
		# recording
		# Besides, it gives the process uniformity across 
		# different types of electrodes
		cable_i = ws.nNAELC + self.axon
		k = 0
		if "node of Ranvier" in settings.keys():
			for j, sec in enumerate(ws.sections[cable_i]):
				secname = sec.name()
				if "NODE" in secname or "node" in secname:
					if k == self.node:
						break
					k += 1

		self.point = {
				'cable': cable_i, 
				'section': j, 
				'z': self.position
			}

		self.set_stim_info()


class ExtracellularPointElectrode(PointElectrode):
	"""
	Extracellular point electrode for intracellular current injections, 
	recordings or connections to ground

	Input parameters:

	z: Position along the z-axis
	"""
	def __init__(self, name, index=0):

		PointElectrode.__init__(self, name, index)
		settings = self.settings

		try:
			self.xyz = settings['point']
		except KeyError:
			self.xyz = settings['points'][self.index]
		self.x, self.y, self.z = self.xyz
		
		# Find the cable this corresponds to
		dists = geo.dist((self.x, self.y), (ws.nvt.x, ws.nvt.y))
		these = np.where(dists == dists.min())[0]
		# print('Points:', these)
		# Choose only one cable
		cable_i = these[0]

		# print('Cable:', ws.nvt.cables[cable_i])

		# Locate the section and segment (z value) of the cable
		j_sec, z_on_j = anatomy.locate(ws.seclens[cable_i], self.z)

		# print('Section and segment:', j_sec, z_on_j)

		# Cable and point of the cable to work on
		self.point = {
				'cable': cable_i, 
				'section': j_sec, 
				'z': z_on_j
			}

		self.set_stim_info()


class ExtracellularPointElectrodeSet():
	"""
	Set of ExtracellularPointElectrode instances
	"""
	def __init__(self, name):
		self.name = name
		self.role = ws.electrodes_settings[name]["role"]
		self.type = "extracellular points set"
		self.electrodes = []

	def set_stimulation(self):
		""" Create the list of extracellular point electrodes """
		for i, (x, y, z) in ws.electrodes_settings[self.name]["points"].items():
			# Create the ExtracellularPointElectrode instance and store 
			# it
			pe = ExtracellularPointElectrode(self.name, i)
			pe.set_stimulation()
			self.electrodes.append(pe)

class GroundPointSet():
	"""
	Set of points which are connected to ground
	"""
	def __init__(self, name):
		self.name = name
		self.role = "ground"
		self.type = "ground points set"

	def set_stimulation(self):
		""" Connect grounds """

		# Go straigth to connect the wanted points to ground
		for (x, y, z) in ws.electrodes_settings[self.name]["points"]:

			# Find the cable this corresponds to
			dists = geo.dist((x, y), (ws.nvt.x, ws.nvt.y))
			these = np.where(dists == dists.min())[0]
			# print('Points:', these)
			# Choose only one cable
			cable_i = these[0]

			# print('Cable:', ws.nvt.cables[cable_i])

			# Locate the section and segment (z value) of the cable
			j_sec, z_on_j = anatomy.locate(ws.seclens[cable_i], z)

			# Cable and point of the cable to work on
			point = {
					'cable': cable_i, 
					'section': j_sec, 
					'z': z_on_j
				}

			# Ground it
			ground_segment(point)

		
########################################################################
# Once the classes are defined, define the types of electrodes
electrode_types = {
	"cuff": CuffElectrode, 
	"intracellular": IntracellularElectrode, 
	"extracellular point": ExtracellularPointElectrode, 
	"extracellular points set": ExtracellularPointElectrodeSet, 
	"ground points set": GroundPointSet
}

########################################################################
# Functions

def create_electrodes():
	""" Create the electrodes which are specified in the 
	electrodes.json file """
	ws.electrodes = {}	
	for name, el in ws.electrodes_settings.items():
		print('Electrode: %s'%name)
		ws.electrodes[name] = electrode_types[el['type']](name)

def set_stimulation():
	""" Set up the stimulation protocols described in the json files. 
	This information is already in the electrodes attributes """
	for electrode in ws.electrodes.values():
		if electrode.role == "stimulation":
			electrode.set_stimulation()

def set_recordings():
	""" Set up the stimulation protocols described in the json files. 
	This information is already in the electrodes attributes """

	# # Time
	# set_time_recording()

	# Electrode recordings
	for electrode in ws.electrodes.values():
		if electrode.role == 'recording':
			# Tell ws that there are electrode recordings
			ws.settings["electrode recordings"] = True
			# Set the recordings
			electrode.set_recordings()

def check_cuff_cover(zi):
	"""
	Check if a value over the z-axis is under the region of a cuff 
	electrode
	"""
	for name, electrode in ws.electrodes.items():
		if electrode.type == 'cuff':
			if electrode.left_end <= zi <= electrode.rght_end:
				return name
	return None

def prepare_ground_paths():
	""" 
	Calculate the variables that need to be used to connect the outer 
	axons to the distant ground, through the cuff insulators and the 
	container 
	"""

	# Shortennings
	nvt = ws.nvt
	container_rho = ws.container_settings['resistivity']
	container_diam = ws.container_settings['diameter']
	# Convert the resistivity to Ohm*um
	container_rho *= 1.e4

	# Bases of the cylinder: conductivity to ground
	# Conductances p.u.area in both cases (mho/cm2)
	# Just direct grounding
	ws.xg_distal = 1.e9
	ws.xg_proxim = 1.e9
	# Actually, I want to think that there's medium at the bases
	# This base will have the cross-sectional area of the container, 
	# not the nerve, so I need to use that to compute rg
	# ASSUMPTION *1
	# Let's assume its thickness is the same as in the radial direction
	ws.xg_distal = []
	ws.xg_proxim = []
	# Iterate over cables
	for i in range(ws.counters['cables']):
		area = nvt.free_areas[i]
		# Area that corresponds to this cable for the base of the 
		# container
		# NOTE: This is just an approximation that would be ideal 
		# should the nerve be circular in its cross-section
		# If a nerve is not, it's still OK, because I use its 
		# circumcircle
		area_augmented = area * (container_diam / nvt.circumdiameter) ** 2
		# Conductance in mho
		# Note that the 'depth' (distance from segment to ground along 
		# the z-axis, which means through the saline bath) is the 
		# container's radius, according to our assumption above (ASSUMPTION *1)
		ws.xg_distal.append(area_augmented / (container_rho * 0.5 * container_diam))
		ws.xg_proxim.append(ws.xg_distal[i])

	# Thickness of the containing medium. There's one for each electrode
	ws.container_settings['thickness'] = {}

	# Round wall of the cylinder

	# Area that the current from one segment to ground has to cross (um2)
	# Length of a segment times its portion of the perimeter
	area = (ws.length / ws.NAELC_nseg) \
	* np.pi * 0.5 * (nvt.circumdiameter + container_diam) / ws.npc
	
	# Resistance to ground p.u.area (Ohm*um2)
	# The distances must be in cm (they have to be converted from um)
	# No, it now must be in um again, becuase I've put container_rho in Ohm*um
	rg = {}
	rg['wout_cuff'] = container_rho * 0.5 * (container_diam - nvt.circumdiameter)
	rg['with_cuff'] = {}

	# Conductance to ground in mho (p.u.area (mho/um2) * area (um2))
	ws.xg_rwall = {}

	# Without cuff

	ws.xg_rwall['wout_cuff'] = area / rg['wout_cuff']

	# With cuff
	ws.xg_rwall['with_cuff'] = {}

	# Iterate over cuff electrodes
	for name, el in ws.electrodes.items():
		if el.type == 'cuff':

			# Add the thickness from the electrodes thickness
			ws.container_settings['thickness'][name] = \
				0.5 * container_diam - (nvt.circumradius \
				 + ws.electrodes_settings[name]['thickness'])

			# Round wall of the cylinder

			# Resistance to ground p.u.area (Ohm*um2)
			# The resistivity has to be converted from Ohm*cm to Ohm*um
			rg['with_cuff'][name] = \
				1.e4 * ws.electrodes_settings[name]['insulator resistivity'] \
				* ws.electrodes_settings[name]['thickness'] \
				+ container_rho * ws.container_settings['thickness'][name]

			# Conductance to ground in mho (p.u.area (mho/um2) * area (um2))

			# Conductance
			ws.xg_rwall['with_cuff'][name] = area / rg['with_cuff'][name]

def square_pulses(injections_, point, protocol):
	""" Assign current injections to the pads """
	
	# Injection segment in hoc
	seg = ws.sections[point['cable']][point['section']](point['z'])

	# IClamp object(s)
	for i, amp in enumerate(protocol['currents']):
		# Gather the data
		dur = protocol['pulse durations'][i]
		delay = protocol['pulse onset times'][i]
		# Create IClamp object and add its properties
		injection = h.IClamp(seg)
		# Intensity from uA to nA
		injection.amp = amp * 1.e3
		injection.dur = dur
		injection.delay = delay
		# Add this information as an attribute
		injections_.append({
			'injection id': ws.counters['current injections'], 
			'point': point, 
			'amp': amp, 
			'dur': dur, 
			'delay': delay
		})
		# Update the counter
		ws.counters['current injections'] += 1
		# Add the duration and delay to ws.totaldur
		ws.totaldur = max(ws.totaldur, delay + dur)
		# Append it to ws
		ws.injections.append(injection)
		ws.injections_info = ws.injections_info + injections_
		ws.log('Injection added to: %s, %s, %s, %s'%(str(point), str(amp), str(dur), str(delay)))

def ground_segment(point):
	""" This connects a segment directly to ground. It is obviously not 
	a current injection protocol, but still is a necessary kind of 
	protocol in case it's needed, I think. And this is a good place to 
	define it.
	Because the paths or connections to ground have already been 
	defined in this module previously and declared in h in 
	biophysics.electrics, this replaces whatever path to ground this 
	segment had, if any. If it didn't have it, it creates it. """

	# Corresponding segment in hoc
	seg = ws.sections[point['cable']][point['section']](point['z'])

	# Level or layer where to point to
	level = ws.cables[point['cable']].properties['extlayers'] - 1

	# Ground it directly
	seg.xg[level] = 1.e9



def set_recording(recordings_, point, var):
	""" Assign the recording vectors to each pad """

	# This is the contact segment in hoc
	seg = ws.sections[point['cable']][point['section']](point['z'])

	# Set up the recording information
	recordings_.update({
			'recording id': ws.counters['recordings'], 
			'point': point, 
			'variable': var
		})

	# Set up the recording object
	als.set_recording(seg, var)
