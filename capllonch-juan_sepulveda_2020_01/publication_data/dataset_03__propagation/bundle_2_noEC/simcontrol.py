"""
This module contains the functions that set up a simulation
"""

import os
import sys
import numpy as np
import json
from neuron import h
from datetime import datetime
import time

import electrodes
import anatomy
import analysis as als
import biophysics as bio
import workspace as ws
import tools as tls



def get_paths():
	""" Get the paths to the different folders used by the code. 
	Note: In the future, I will prepare this function to get all the 
	paths in a recursive way. For now, the relevant sub-paths are hard-coded """

	# Dictionary with folders
	ws.folders = {}

	# Current Working Directory
	ws.folders["cwd"] = os.getcwd()

	# Get all the folders in cwd and sub-folders

	_, subdirs = tls.explore_folder(ws.folders["cwd"])

	# for item in os.listdir(ws.folders["cwd"]):
	# 	if os.path.isdir(item):
	for subdir in subdirs:

		# Add this folder
		path = os.path.join(ws.folders["cwd"], subdir)
		ws.folders[subdir] = path

		# Check if there's no more sub-folders here
		_, subsubdirs = tls.explore_folder(path)
		wholeend = len(subsubdirs) == 0

		# If it's not a 'whole end', dig on
		for ssd in subsubdirs:
			path = os.path.join(ws.folders["cwd"], subdir)
			# Do nothing for now
			# Later on, I'll finish this function

	# Hard-coded part: sub-folders
	ws.folders["data/load"] = os.path.join(ws.folders["data"], "load")
	ws.folders["data/saved"] = os.path.join(ws.folders["data"], "saved")
	ws.folders["data/results"] = os.path.join(ws.folders["data"], "results")
	ws.folders["src/hoc"] = os.path.join(ws.folders["cwd"], "src/hoc")

def create_log_file():
	""" Create a log file for the simulation """

	# Create the name
	ws.logfile = os.path.join(ws.folders["data"], "log.txt")

	# Create the file
	f = open(ws.logfile, "w")
	f.close()

def sanity_checks():
	""" Perform some checks to make sure that everything is defined 
	correctly and do whatever is necessary if not """

	# Using the Resistor Network model implies...
	if ws.settings["nerve model"] == "resistor network":

		# that there is ephaptic coupling
		if not ws.EC["presence"]:
			msg = "WARNING: Using the Resistor Network Model implies "\
				 + "the presence of ephaptic coupling\nEphaptic "\
				 + "coupling was activated"
			ws.log(msg)
			ws.EC["presence"] = True
			ws.EC["use"] = "resistor network"

		# that the stimulation forcing can't be "pre-computed"
		if ws.settings["stimulation"]["method"] == "pre-computed":
			msg = "ERROR: Using the Resistor Network Model implies "\
				 + "that the stimulation method can't be pre-computed\n"\
				 + "Please check your settings"
			ws.log(msg)
			ws.terminate()

def read_settings():
	""" Read all the settings from the json files """

	# General settings
	# Open json file containing the information about the chosen axon model
	file = os.path.join(ws.folders["settings"], "settings.json")
	with open(file, "r") as f:
		ws.settings = json.load(f)

	for k, v in ws.settings["ephaptic coupling"].items():
		try:
			ws.settings["ephaptic coupling"][k] = tls.str2bool(v)
		except TypeError:
			# It's probably a dictionary
			for kk, vv in v.items():
				ws.settings["ephaptic coupling"][k][kk] = tls.str2bool(vv)

	# Shortening for the settings for ephaptic coupling
	ws.EC = ws.settings["ephaptic coupling"]

	# String to bool for the rest of settings
	for k, v in ws.settings.items():
		try:
			v_is_str = "False" in v or "True" in v
		except TypeError:
			# Not a string, dictionary or list
			v_is_str = False

		# If it was a string, then parse it
		if v_is_str:
			ws.settings[k] = tls.str2bool(v)

	# Axon model
	# Open json file containing the information about the chosen axon model
	axon_files = {}
	ws.axonmodel_settings = {}
	for m in ws.settings['axon models']:
		af = m + '.json'
		axon_files[m] = af
		with open(os.path.join(ws.folders["models"], af), "r") as f:
			ws.axonmodel_settings[m] = json.load(f)

	###############################################################################
	# PHYSICAL PROPERTIES AND GEOMETRY

	# Anatomy
	# Open json file containing the information about anatomy
	file = os.path.join(ws.folders["settings"], "anatomy.json")
	with open(file, "r") as f:
		ws.anatomy_settings = json.load(f)

	# Electrodes
	file = os.path.join(ws.folders["settings"], "electrodes.json")
	with open(file, "r") as f:
		ws.electrodes_settings = json.load(f)

	# Container
	file = os.path.join(ws.folders["settings"], "container.json")
	with open(file, "r") as f:
		ws.container_settings = json.load(f)

	# Tissues
	file = os.path.join(ws.folders["settings"], "tissues.json")
	with open(file, "r") as f:
		ws.tissues_settings = json.load(f)

	# Check that everything is fine
	sanity_checks()

def physics():
	""" Set up some physical properties of the system from the settings """

	# Length (cm)
	ws.length = ws.anatomy_settings['length']
	# Internodal length to fiberD ratio
	# This will rarely be used in some models
	ws.ilD = ws.anatomy_settings['axons']['myelinated']['internodal length to fiberD ratio']

	# Temperature
	h.celsius = ws.settings['temperature']

	# Resistivities and thicknesses
	# Units:
	# Resistivity: Ohm*cm
	# Distance units: um

	# Resistivities
	ws.rho = {}
	# Thicknesses
	ws.th = {}

	# Nominal thickness of the perineurium (um)
	ws.th['perineurium'] = ws.tissues_settings['thickness']['perineurium']


	# Resistivities of the tissues. They must be in MOhm*um 
	# (convert from Ohm*cm)
	ws.rho['perineurium'] = ws.tissues_settings['resistivities']['perineurium'] * 1.e-6 * 1.e4

	ws.rho['epi_trans'] = ws.tissues_settings['resistivities']['epineurium']['transverse'] * 1.e-6 * 1.e4
	ws.rho['epi_longt'] = ws.tissues_settings['resistivities']['epineurium']['longitudinal'] * 1.e-6 * 1.e4
	ws.rho['endo_trans'] = ws.tissues_settings['resistivities']['endoneurium']['transverse'] * 1.e-6 * 1.e4
	ws.rho['endo_longt'] = ws.tissues_settings['resistivities']['endoneurium']['longitudinal'] * 1.e-6 * 1.e4

	# Perineurium resistance (MOhm*um2)
	ws.resist_perineurium = ws.rho['perineurium'] * ws.th['perineurium']

	ws.log("Perineurium resistance %f MOhm*um2:"%ws.resist_perineurium)


	# INFINITE RESISTIVITY OF AN IDEALLY TOTALLY RESISTIVE MEDIUM
	ws.rho['SuperResistive'] = 1.e9

def miscellanea():
	""" Miscellaneous variables. It's OK to have hard-coded 
	variables here """

	# DISCRETISATION OF THE VOLUME CONDUCTOR

	# Cables
	# NAELC. All this will need to be inside a json file
	# Number of segments for the cables
	ws.nseg = ws.settings["miscellanea"]["NAELC def nseg"]
	# Diameter (um)
	ws.cdiam = 1.
	# Extracellular layers
	ws.NAELC_extlayers = 1

	# Other variables
	ws.raxmin = 1e99
	ws.deltax_min = 1e99
	ws.g_min = 1e99
	# Counters
	ws.counters = {
		"cables": 0, 
		"axons": 0, 
		"NAELC": 0,
		"current injections": 0,
		"recordings": 0
	}
	ws.maxnseg = 0
	ws.totaldur = 0.
	ws.injections = []
	ws.injections_info = []
	ws.recordings = []

def setup_simcontrol():
	""" Try to set up the simulation control if there's a simulation 
	meant to happen """

	try:
		scsettings = ws.settings['simulation control']
		h.v_init = scsettings['v_init']
		ws.dt = h.dt = scsettings['dt']
		ws.nt = scsettings['nt']
		ws.tstop = scsettings['nt'] * h.dt
		h.tstop = ws.tstop
		ws.tarray = np.arange(0., h.tstop, h.dt)
	except LookupError:
		# There's not supposed to be any simulation
		pass

def prepare_workspace(remove_previous=True):
	""" Prepare all the necessary information in the workspace to be 
	able to do the necessary work, including running a simulation """

	# Data management
	# Prepare the folders that the code will use
	get_paths()
	# Remove results from previous simulation(s)
	if remove_previous:
		als.remove_previous_results()
	# Create the log file
	create_log_file()
	# Pass a couple of functions to ws
	ws.log = log
	ws.terminate = terminate
	# Read all the settings from the json files
	read_settings()
	now = datetime.now()
	ws.log("Execution started on: %s"%now.strftime("%d %h %Y at %H:%M:%S"))
	# Set up some physical properties of the system from the settings
	physics()
	# Set up a simulation
	setup_simcontrol()
	# Miscellaneous variables
	miscellanea()

def prepare_simulation(prep_workspace=True):
	""" Do the whole process of preparing a simulation """

	# Prepare the necessary variables from the configuration files
	if prep_workspace:
		prepare_workspace()

	# Shortenings
	stimmethod = ws.settings["stimulation"]["method"]
	nervemodel = ws.settings["nerve model"]

	# Create the NERVE TOPOLOGY. Cross-section only
	anatomy.create_nerve()

	# BIOPHYSICS

	# LOAD NEURON HOC CODE AND DEFINE FUNCTIONS

	# Load NEURON code
	h.load_file(os.path.join(ws.folders["src/hoc"], "ephap.hoc"))
	h.load_file(os.path.join(ws.folders["src/hoc"], "wire.hoc"))
	for m in ws.axonmodel_settings.values():
		h.load_file(m['hoc file'])

	# BUILD ALL THE WIRES with all their length
	bio.build_cables()

	# Finish setting up information for the fascicles
	if nervemodel == "resistor network":
		bio.finish_fascicles()

	# Create electrodes
	# This needs to be after the creation of the cables and before 
	# connecting the ephapses
	if stimmethod == "from electrodes":
		ws.log("Creating the electrodes")
		electrodes.create_electrodes()
	# elif stimmethod == "pre-computed":
	# 	pass

	# EPHAPSES
	if ws.EC["presence"]:
		if nervemodel == "resistor network":

			# Create an instance of the Resistor Network
			ws.log("Creating the Resistor Network")
			ws.rnet = bio.ResistorNetwork()

			# Get the resistances of the connections for each pair
			ws.log("\tGetting resistor values")
			ws.rnet.get_resistances()

			# Connect resistors with ephapses in NEURON
			ws.log("\tConnecting resistors")
			ws.rnet.connect_resistors()

	# OTHER CONNECTIONS

	ws.log("Preparing connections to ground on the boundaries and xraxial")
	# Prepare connections to ground in the system
	if stimmethod == "from electrodes":
		electrodes.prepare_ground_paths()

	# Set xraxial and make the connections to ground
	if nervemodel == "resistor network":
		bio.electrics()

	# STIMULATION
	ws.log("Setting up stimulation")

	# Inject currents from the cuff electrode

	# Add injected currents, delays and durations
	if stimmethod == "from electrodes":
		# if nervemodel == "resistor network":
		electrodes.set_stimulation()
	elif stimmethod == "pre-computed":
		bio.set_direct_extracellular_stimulation()

	# RECORDINGS
	ws.log("Setting up recordings")

	# Time
	als.set_time_recording()

	# From electrodes
	# Flag: presence of electrode recordings
	# Set it to False by default; it will become True if a recording 
	# electrode is created
	ws.settings["electrode recordings"] = False
	if stimmethod == "from electrodes":
		if nervemodel == "resistor network":
			electrodes.set_recordings()

	# Directly on the cables/axons
	als.record_cables()

def log(msg):
	""" Write a message into the log file """
	with open(ws.logfile, "a") as f:
		f.write(msg + "\n")
	# Print if wanted
	if ws.settings["verbose"]:
		print(msg)

def terminate():
	""" Terminate the whole program """
	now = datetime.now()
	log("THE PROGRAM WAS TERMINATED ON: %s"%now.strftime("%d %h %Y at %H:%M:%S"))
	sys.exit()

def run():
	""" Run a simulation """

	if h.tstop <= ws.totaldur:
		msg = 'WARNING: The simulation tstop is shorter than the total '\
		+ 'stimulation time. This will cause trouble when processing '\
		+ 'the recordings.'
		ws.log(msg)

	def initialize():
		""" Initialize the simulation for the first time step """
		h.finitialize(h.v_init)
		h.fcurrent()

	def go_by_steps():
		""" Main driver for the simulation """
		initialize()
		# Loop over time steps
		for i_t, t in enumerate(ws.tarray):
			h.fadvance()
			# Rest for 500 ms to allow the computer to cool down
			if False:
				time.sleep(0.5)
			# Time control
			t1 = datetime.now()
			msg = "Time step %i of %i. Total elapsed time: %i seconds"%(i_t, ws.nt, (t1 - t0).total_seconds())
			# Log the time step into the simulation log file
			ws.log(msg)

	def go_full():
		""" Run the simulation in one go with h.run() """
		initialize()
		h.run()

	# Run
	t0 = datetime.now()
	# log('Running the simulation in NEURON')
	now = datetime.now()
	ws.log("Starting the simulation in NEURON on: %s"%now.strftime("%d %h %Y at %H:%M:%S"))
	go_by_steps()
	now = datetime.now()
	ws.log("NEURON simulation finished on: %s"%now.strftime("%d %h %Y at %H:%M:%S"))
	log('Whole simulation time: %i seconds'%(now - t0).total_seconds())

