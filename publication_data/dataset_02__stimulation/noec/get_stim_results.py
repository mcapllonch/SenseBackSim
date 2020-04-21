"""
Print the relevant results of the stimulation simulations
"""

import numpy as np
from scipy.signal import argrelextrema
from collections import OrderedDict
import json
import csv
import sys
import os

import workspace as ws
import geometry as geo


def process_results(path, ec_key, dataset_name):
	""" Read and process the results from this dataset"""
	print(path, dataset_name, ec_key)

	# Results folder
	results_folder = os.path.join(path, 'data/results')
	# List all the items in the results folder
	items = sorted(os.listdir(results_folder))
	# Select the csv files
	items = [item for item in items if ".csv" in item]
	# Select the axon recording files
	items = [item for item in items if item[:4] == "Axon"]

	# Array for the axons' activity maxima
	maxima = []

	# AP peak times
	appt_ = {}
	# AP latency times
	aplt_ = {}

	# Flags indicating which axons did fire and which not
	hasAP = {}

	# Voltage data
	data = {}

	# Geometrical properties
	xx_ = []
	yy_ = []
	rr_ = []


	# Iterate over the files and read them
	for filename in items:
		# Actually, i is taken from the file name
		i = int(filename.replace('Axon', '').replace('.csv', ''))
		data[i] = {}
		with open(os.path.join(results_folder, filename), "r") as f:
			fr = csv.reader(f)
			for row in fr:
				r0 = row[0]
				if ("NODE" in r0) or ("node" in r0):
					data[i][key].append([float(x) for x in row[1:]])
				elif len(row) == 3:
					xx_.append(float(r0))
					yy_.append(float(row[1]))
					rr_.append(float(row[2]))
				elif len(row) == 1:
					try:
						# print(key)
						data[i][key] = np.array(data[i][key])
					except NameError:
						# There's no key yet
						pass
					key = r0
					data[i][key] = []

		# When the last key is read, don't forget storing it
		data[i][key] = np.array(data[i][key])
		del key

	# Check maxima and relevant stuff
	vcrit = 15.
	for i, data_ in data.items():
		axondata = data_["v"]
		# Check if the maximum is an AP
		# print(axondata)
		maximum = axondata.max()
		maxima.append(maximum)
		if maximum > 0:
			# Regions where v is greater than vcrit
			whereAPs = np.where(axondata > vcrit)
			# Time when the first AP is fired (v rises above vcrit mV)
			when1stAP = whereAPs[1].min()
			where1stAP = whereAPs[0][np.where(whereAPs[1] == when1stAP)][0]
			segment_maxima = argrelextrema(axondata[where1stAP], np.greater)[0]
			# Local maxima
			local_maxima = axondata[where1stAP][segment_maxima]
			# Local maxima greater than vcrit mV
			# IMPORTANT: I named the following variable when1stAP_ just so 
			# it doesn't overwrite when1stAP, but I can make it overwrite 
			# it if I want
			when1stAP_ = segment_maxima[np.where(local_maxima > vcrit)][0]
			if True:
				APpeaktime = dt * (when1stAP - 1)
				appt_[i] = APpeaktime
				aplt_[i] = APpeaktime - 0.01
			hasAP[i] = True
		else:
			hasAP[i] = False
			aplt_[i] = 'nan'
		i += 1

	# Maxima to array
	maxima = np.array(maxima)

	# AP peak times to array
	# And subtract the pulse delay from them
	appt_values = np.array(list(appt_.values())) - 0.01

	if len(appt_values) == 0:
		# print("No axon fired an AP")
		pass
		# sys.exit()

	# Geometrical properties to array
	xx_ = np.array(xx_)
	yy_ = np.array(yy_)
	rr_ = np.array(rr_)

	# Topology

	# Open and read topology file
	topo_path = os.path.join(path, 'data/load/created_nerve_internal_topology.json')
	with open(topo_path, 'r') as f:
		topology = json.load(f)

	# Open and read the contours file
	contours_file = os.path.join(path, 'data/load/created_nerve_contour.csv')
	contours = {}
	with open(contours_file, "r") as f:
		fr = csv.reader(f)
		for row in fr:
			try:
				# Try to get numbers
				x = float(row[0])
			except ValueError:
				# It's a string
				key = row[0]
				contours[key] = []
			else:
				y = float(row[1])
				# Append the point to the corresponding contour
				contours[key].append([x, y])

	# Delete the key for tidyness
	del key

	# Polygons for each fascicle and the nerve
	polygons = {}
	for k, c in contours.items():
		pol = geo.Polygon(c)
		polygons[k] = pol.plpol

	# Fired axons per fascicle and in total in the nerve
	fascicle_ap_counter = {}
	for k in contours:
		if 'ascicle' in k:
			fascicle_ap_counter[k] = 0

	# Find them
	for i, b in hasAP.items():
		if b:
			# Find fascicle of this axon
			for k, p in polygons.items():
				if 'ascicle' in k:
					if p.contains_point((xx_[i], yy_[i])):
						# print('Axon %i is in %s'%(i, k))
						fascicle_ap_counter[k] += 1
						break

	# Read electrode settings
	settings_path = os.path.join(path, 'settings/electrodes.json')
	with open(settings_path, 'r') as f:
		stim_settings = json.load(f)

	current = list(list(stim_settings.values())[0]['stimulation protocol'].values())[0]['currents'][0]

	total_number_APs = sum(list(fascicle_ap_counter.values()))

	# Dictionary to gather all the important data
	data_final = OrderedDict()
	data_final['dataset_name'] = dataset_name
	data_final['current'] = current
	data_final['fascicle_ap_counter'] = fascicle_ap_counter
	data_final['total_number_APs'] = total_number_APs
	data_final['AP times'] = aplt_

	# Save the data in the 'all data' dictionary
	data_all[dataset_name] = data_final

	# Save data into the data_recruitment dictionary
	data_recruitment['currents'].append(current)
	data_recruitment['recruitment'][ec_key]['nerve'].append(total_number_APs)
	for k, n in fascicle_ap_counter.items():
		data_recruitment['recruitment'][ec_key][k].append(n)

	# Save results in a json file
	with open('stim_results_%s%s'%(dataset_name, '.json'), 'w') as f:
		json.dump(data_final, f, indent=4)


# Time
dt = 0.005
	
# Cwd
cwd = os.getcwd()

# Items
cwd_items = sorted(os.listdir(cwd))

# All data dictionary
data_all = {}

# Dictionary for all the recruitment data for the article
data_recruitment = {
	'currents': [], 
	'recruitment': 
		{}
}

for s in ('EC', 'noEC'):
	data_recruitment['recruitment'][s] = {
		'nerve': [], 
		'Fascicle_00': [], 
		'Fascicle_01': [], 
		'Fascicle_02': [], 
		'Fascicle_03': [], 
		'Fascicle_04': [], 
		'Fascicle_05': [], 
		'Fascicle_06': []
	}

# Explore and find the relevant directories
dataset_dirs = []
for item in cwd_items:
	sd_ = list(os.walk(os.path.join(cwd, item)))
	for x in sd_:
		x = x[0]
		xsp = x.split('EC')
		if len(xsp[-1]) <= 0:
			dataset_dirs.append(x)
			# Print the dataset name
			path_elements = [s for s in x.split('/') if ('current' in s) or ('EC' in s)]
			ec_key = [s for s in path_elements if 'EC' in s][0]
			dataset_name = '_'.join(path_elements)
			process_results(x, ec_key, dataset_name)

# Cwd
cwd = os.getcwd()

# Items
cwd_items = os.listdir(cwd)

# # Save the 'all data' and recruitment dictionaries
# with open('recruitment_data.json', 'w') as f:
# 	json.dump(data_recruitment, f, indent=4)

# Print the recruitment data on screen
data_recr_file = open('recruitment_data.txt', 'w')

print('currents:')
data_recr_file.write('currents:\n')
for c in data_recruitment['currents']:
	print(c)
	data_recr_file.write('%0.1f\n'%c)
print('')
data_recr_file.write('\n')
print('recruitment:')
data_recr_file.write('recruitment:\n')
for k, l in data_recruitment['recruitment']['noEC'].items():
	print('')
	data_recr_file.write('\n')
	print(k)
	data_recr_file.write('%s\n'%k)
	for x in l:
		print(x)
		data_recr_file.write('%i\n'%x)

data_recr_file.close()