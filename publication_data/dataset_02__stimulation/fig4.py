"""
Compare two datasets, or two simulations.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from collections import OrderedDict
from gi.repository import Gtk
import zipfile
import math
import json
import csv
import sys
import os



def round_up(x, order=1):
	""" Round a number up to the given order of magnitude """
	return order * math.ceil(float(x) / order)
def round_down(x, order=1):
	""" Round a number down to the given order of magnitude """
	return order * math.floor(float(x) / order)

# Cwd
cwd = os.getcwd()

# Items
cwd_items = os.listdir(cwd)

# List zip files
zipflist = [
		"revs_round_2_stim_set01_current_02000nA_EC.zip", 
		"revs_round_2_stim_set01_current_02000nA_noEC.zip"
	]

print(zipflist)

# Datasets where to store all data
datasets = OrderedDict()

# Initialize some variables that are needed for the graphs
vmin = 1e99
vmax = -1e99

# Dictionary that identifies the znames with the ephaptic coupling
znames = {}

# Iterate over zip files
for zf in zipflist:
	print(zf)
	zname = zf.replace(".zip", "")

	# Create a directory for it if it doesn't exist
	isincwd = zname in cwd_items

	if not isincwd:
		# Create the directory
		os.mkdir(zname)
		# Reference the folder
		zfolder = os.path.join(cwd, zname)
		# Extract everything in it
		zip_ref = zipfile.ZipFile(zf, 'r')
		zip_ref.extractall(zfolder)
		zip_ref.close()

	# Reference the folder
	zfolder = os.path.join(cwd, zname)

	# Now process data
	datasets[zname] = OrderedDict()

	# Axons' activity

	# Results folder
	results_folder = os.path.join(zfolder, "data/results")
	# List all the items in the results folder
	items = sorted(os.listdir(results_folder))
	# Select the csv files
	items = [item for item in items if ".csv" in item]
	# Select the axon recording files
	items = [item for item in items if item[:4] == "Axon"]

	# Array for the axons' activity maxima
	maxima = []

	# Time
	dt = 0.005

	# AP peak times
	appt_ = {}

	# AP peak places
	APpeakplace = {}

	# Flags indicating which axons did fire and which not
	hasAP = []

	# Voltage data
	data = {}

	# Iterate over the files and read them
	i = 0
	for filename in items:
		data[i] = {}
		with open(os.path.join(results_folder, filename), "r") as f:
			fr = csv.reader(f)
			for row in fr:
				if "NODE" in row[0]:
					data[i][key].append([float(x) for x in row[1:]])
				elif len(row) == 3:
					pass
				elif len(row) == 1:
					try:
						data[i][key] = np.array(data[i][key])
					except NameError:
						# There's no key yet
						pass
					key = row[0]
					data[i][key] = []

		# When the last key is read, don't forget storing it
		data[i][key] = np.array(data[i][key])
		del key
		# Store data in the dictionary
		i += 1

	# Check maxima and relevant stuff
	for i, data_ in data.items():
		axondata = data_["v"]
		# Check if the maximum is an AP
		maximum = axondata.max()
		maxima.append(maximum)
		if maximum > 0:
			# Regions where v is greater than 15
			whereAPs = np.where(axondata > 15)
			# Time when the first AP is fired (v rises above 15mV)
			when1stAP = whereAPs[1].min()
			where1stAP = whereAPs[0][np.where(whereAPs[1] == when1stAP)][0]
			segment_maxima = argrelextrema(axondata[where1stAP], np.greater)[0]
			# Local maxima
			local_maxima = axondata[where1stAP][segment_maxima]
			# Local maxima greater than 15 mV
			# IMPORTANT: I named the following variable when1stAP_ just so 
			# it doesn't overwrite when1stAP, but I can make it overwrite 
			# it if I want
			when1stAP_ = segment_maxima[np.where(local_maxima > 15)][0]
			APpeaktime = dt * (when1stAP - 1)
			APpeakplace[i] = where1stAP
			appt_[i] = APpeaktime
			hasAP.append(True)
		else:
			hasAP.append(False)
		i += 1

	# Maxima to array
	maxima = np.array(maxima)

	# AP peak times to array
	# And subtract the pulse delay from them
	appt_values = np.array(list(appt_.values())) - 0.01

	if len(appt_values) == 0:
		print("No axon fired an AP")
		sys.exit()

	# Delete the key for tidyness
	del key

	# Store everything in the datasets dictionary
	datasets[zname]["APpeakplace"] = APpeakplace
	datasets[zname]["appt_"] = appt_
	datasets[zname]["appt_values"] = appt_values
	datasets[zname]["hasAP"] = hasAP
	datasets[zname]["figrow"] = zipflist.index(zf)
	datasets[zname]["data"] = data
	if "_EC" in zname:
		datasets[zname]["title"] = "With Ephaptic \nCoupling (sEC)"
		datasets[zname]["histcolor"] = "r"
		znames["EC"] = zname
	elif "noEC" in zname:
		datasets[zname]["title"] = "Without Ephaptic \nCoupling (snoEC)"
		datasets[zname]["histcolor"] = "b"
		znames["noEC"] = zname

	# Common variables
	vmin = min(vmin, appt_values.min())
	vmax = max(vmax, appt_values.max())

########################################################################
# Process vext[1]

# Get the initial value (for it = 0 or 1), which is what is used for 
# extstim
extstim = {}
vext_EC = {}
vext_noEC = {}
vext_noEC_mod = {}
vext_diffs = {}
vm_EC = {}
vm_noEC = {}
# Numbers of fired axons
nfa_sEC = 0
nfa_snoEC = 0
# Newly fired axons
newlyfa_total = 0



# Choose the node in the middle, which is roughly the closest to the 
# stimulating pad (not really, because that is on the first of four 
# rings)

dataset = datasets[znames["EC"]]
data = dataset["data"]
APpeakplace = dataset["APpeakplace"]
APpeakplace_noEC = datasets[znames["noEC"]]["APpeakplace"]

# Get the fields
for i, variables in data.items():
	field_EC = variables['vext[1]']
	field_noEC = datasets[znames["noEC"]]["data"][i]['vext[1]']
	transmv_EC = variables['v']
	transmv_noEC = datasets[znames["noEC"]]["data"][i]['v']

	# Chose the node of Ranvier that lies closest to the source
	try:
		n = APpeakplace[i]
	except KeyError:
		# This axon has no AP in any simulation
		n = field_EC.shape[0] // 2
	else:
		# It has one in sEC, so now check if there is an AP in snoEC
		nfa_sEC += 1
		try:
			_ = APpeakplace_noEC[i]
		except KeyError:
			# This axon has no AP in snoEC
			# Save the axon
			newlyfa_total += 1

	# Count the number of fired axons in snoEC
	try:
		_ = APpeakplace_noEC[i]
	except KeyError:
		pass
	else:
		nfa_snoEC += 1

	# Store the data into the dictionaries
	vext_EC[i] = field_EC[n]
	vext_noEC[i] = field_noEC[n]
	# Modify field_noEC to chop the last time step of the pulse from it
	# This is a way to compute the difference between the fields from the two 
	# simulations without getting spurious values
	field_noEC_mod = np.zeros_like(field_noEC[n])
	where = np.where(field_noEC[n] < -5.e-1)[0][:-1]
	field_noEC_mod[where] = field_noEC[n][where]
	vext_noEC_mod[i] = field_noEC_mod

	# Go on with the others
	vm_EC[i] = transmv_EC[n]
	vm_noEC[i] = transmv_noEC[n]
	extstim[i] = vext_noEC[i][3]
	vext_diffs[i] = vext_EC[i] - vext_noEC_mod[i]

print("Fired axons in sEC: %i"%nfa_sEC)
print("Fired axons in snoEC: %i"%nfa_snoEC)
print("Total number of newly fired axons: %i"%newlyfa_total)
print("Sanity check: %i"%(nfa_sEC - nfa_snoEC))

# Compound variables

# vext_diffs from all the axons grouped in one array
vext_diffs_array = np.array([vd for vd in vext_diffs.values()])
# Limits
vext_diffs_min = vext_diffs_array.min()
vext_diffs_max = vext_diffs_array.max()
# Mean and standard deviation
vda_mean = vext_diffs_array.mean(axis=0)
vda_std = vext_diffs_array.std(axis=0)
# Percentiles
vda_perc25 = np.percentile(vext_diffs_array, 25, axis=0)
vda_perc75 = np.percentile(vext_diffs_array, 75, axis=0)
vda_perc10 = np.percentile(vext_diffs_array, 10, axis=0)
vda_perc90 = np.percentile(vext_diffs_array, 90, axis=0)

print("Peak of averaged fields strengths: %f mV"%vda_mean.min())
print("Peak of all fields strengths: %f mV"%vext_diffs_array.min())
########################################################################
# Provisional figures
plt.close('all')

# Cut the time domain in half
nt = len(vext_EC[0]) // 2
time = dt * np.arange(nt)

# Choose n*n random axons
n = 4
these_axons = np.random.choice(list(extstim.keys()), n * n)
print('randomly selected axons:')
print(these_axons)
# Choose one specific axon
chosen_axon = 507

########################################################################
# Show figure

plt.style.use('./PaperFig_1textwidth_v5.mplstyle')

fig, ax = plt.subplots(1, 3)
ax = ax.flatten()

# Potentials
ax[0].plot(time, vm_noEC[chosen_axon][:nt], c='b', alpha=0.5)
ax[0].plot(time, vm_EC[chosen_axon][:nt], c='k')
ax[1].plot(time, vext_noEC_mod[chosen_axon][:nt], c='b', alpha=0.5)
ax[1].plot(time, vext_EC[chosen_axon][:nt], c='k')
ax[0].set_ylabel(r"$\rm v_{m}$ $\rm (mV)$")
ax[1].set_ylabel(r"$\rm v_{E}$ $\rm (mV)$")
ax[0].set_yticks(np.arange(-80, 80, 40))
ytick_separation = 50
yticks = np.arange(
		round_down(vext_EC[chosen_axon][:nt].min(), order=10), 
		round_up(vext_EC[chosen_axon][:nt].max(), order=10) + ytick_separation, 
		ytick_separation
	)

ax[1].set_yticks(yticks)
ax[1].set_ylim(yticks.min(), yticks.max())

# Plot all the vext_diff vs time
for vd in vext_diffs.values():
	ax[2].plot(time, vd[:nt], c='k', alpha=0.2)
	ax[2].plot(time, vda_mean[:nt], c='r', lw=1., alpha=1.)
	ax[2].plot(time, vda_mean[:nt] - vda_std[:nt], c='r', lw=0.2, ls='-', alpha=0.3)
	ax[2].plot(time, vda_mean[:nt] + vda_std[:nt], c='r', lw=0.2, ls='-', alpha=0.3)
	ax[2].axvline(time[where[0]], ls='-', c='k', lw=0.2, alpha=0.3)
	ax[2].axvline(time[where[-1]], ls='-', c='k', lw=0.2, alpha=0.3)

ax[1].set_xlabel("Time " + r"$\rm (ms)$")
ax[2].set_ylabel(r"$\rm v_{E}$ Diff. $\rm (mV)$")
ytick_separation = 30
yticks = np.arange(
		round_down(vext_diffs_min, order=10), 
		round_up(vext_diffs_max, order=10) + ytick_separation, 
		ytick_separation
	)
ax[2].set_yticks(yticks)
xticks = np.arange(time.min(), time.max() + 0.1, 0.1)
xticklabels = ["%0.1f"%x for x in xticks]
ax[0].set_xticks(xticks)
ax[1].set_xticks(xticks)
ax[2].set_xticks(xticks)
ax[0].set_xticklabels(xticklabels)
ax[1].set_xticklabels(xticklabels)
ax[2].set_xticklabels(xticklabels)
ax[0].set_xlim(time.min(), time.max())
ax[1].set_xlim(time.min(), time.max())
ax[2].set_xlim(time.min(), time.max())

########################################################################
# Save and show

fig.savefig("fig4.png", dpi=300)
plt.show()