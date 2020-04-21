"""
Get the extracellular fields over the axons given a certain stimulation 
protocol
"""

import os
import csv
from neuron import h, gui
from collections import OrderedDict

import workspace as ws
import simcontrol
import analysis as als
import biophysics as bio



# Prepare the necessary variables from the configuration files
simcontrol.prepare_workspace()

# Change the number of time steps and refresh the settings for the 
# simulation control
ws.settings['simulation control']['nt'] = 1

# Manually modify some options in ws.settings so that this uses the 
# Resistor Network model and the electrodes defined in electrodes.json
ws.settings["nerve model"] = "resistor network"
ws.settings["ephaptic coupling"]["presence"] = "True"
ws.settings["stimulation"] = \
	{
		"method": "from electrodes", 
		"file": "None"
	}

# Pulses dictionary and counter
pulses = {}
j = 0
# Iterate over stimulating electrodes
for k, el in ws.electrodes_settings.items():
	if el["role"] == "stimulation":
		# print(k)

		# Iterate over active pads
		for pad, protocol in el["stimulation protocol"].items():
			# print("\t", pad, protocol["type"])

			# Iterate over pulses on this pad
			if protocol["type"] == "square pulses":
				for i, (o, d) in enumerate(zip(
								protocol["pulse onset times"], 
								protocol["pulse durations"])):

					# Add its properties to the dictionary
					pulse = {}
					pulse["electrode name"] = k
					pulse["pad"] = pad
					pulse["type"] = "square pulse"
					pulse["delay"] = o
					pulse["duration"] = d
					pulse["pulse number for this pad"] = i
					pulses[j] = pulse

					# Update the pulse counter
					j += 1

					# For each pulse, trick the model into thinking that 
					# the pulse virtually won't happen
					# This is necessary to isolate all the pulses from 
					# each other
					ws.electrodes_settings[k]["stimulation protocol"] \
						[pad]["pulse onset times"][i] = 2. * ws.tstop

# Save the original tstop
tstop_original = ws.tstop

# File to save the data in
espath = os.path.join(ws.folders['data/load'], 'extstim.csv')


# Now iterate over pulses and compute the field over 
# the axons for each pulse
for i, (k, p) in enumerate(pulses.items()):
	# print(k, p)

	# Empty the current injections list and reset the counters
	simcontrol.miscellanea()

	# Turn the delay for this pulse to zero so the stimulation 
	# starts inmediately
	ws.electrodes_settings[p["electrode name"]]["stimulation protocol"] \
		[p["pad"]]["pulse onset times"] \
		[p["pulse number for this pad"]] = 0.

	# Prepare the necessary stuff for the simulation
	simcontrol.prepare_simulation(prep_workspace=False)
	simcontrol.setup_simcontrol()

	# Run the simulation
	simcontrol.run()

	# Save the values into the file

	# If this is the first pulse, open the file as a new one
	mode = "a"
	if k == 0:
		mode = "w"
	with open(espath, mode) as f:

		# Write the number and type of pulse
		fw = csv.writer(f)
		fw.writerow(["Pulse", i])
		fw.writerow(["Type", p["type"]])
		fw.writerow(["Delay", p["delay"]])
		fw.writerow(["Duration", p["duration"]])

		# Get the fields for each axon
		nvt = ws.nvt
		for i, j in enumerate(range(nvt.nNAELC, nvt.nc, 1)):
			exst = [seg.vext[ws.axonmodel_settings[nvt.models[j]]['extlayers']-1] \
				for seg in ws.segments[j]]
			fw.writerow([i] + exst)

	# Once, finished, turn the delay for this pulse back to later than 
	# tstop so it doesn't interfere with the following pulses
	ws.electrodes_settings[p["electrode name"]]["stimulation protocol"] \
		[p["pad"]]["pulse onset times"] \
		[p["pulse number for this pad"]] = 2. * tstop_original