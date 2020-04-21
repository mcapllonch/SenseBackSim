"""
Open the results csv files, read their information and return it
"""

import csv
import os
import numpy as np
from collections import OrderedDict



possible_secnames = [
	"node", 
	"NODE", 
	"FLUT", 
	"STIN", 
	"MYSA", 
	"Wire"
]

def read_results(folder):
	""" Just read all the results """

	# Select the csv files
	items = sorted(os.listdir(folder))
	items = [item for item in items if ".csv" in item]

	# Data storage variables	
	# Voltage data
	data = {}
	# Geometrical properties
	xx_ = []
	yy_ = []
	rr_ = []
	# z-axis anatomy
	zprofile = {}

	# Iterate over the files and read them
	i = 0
	for filename in items:
		j = 0
		data[i] = {}
		zprofile[i] = OrderedDict()
		with open(os.path.join(folder, filename), "r") as f:
			fr = csv.reader(f)
			for row in fr:
				isasection = np.array([x in row[0] for x in possible_secnames]).any()
				if isasection:
					data[i][key].append([float(x) for x in row[1:]])
				elif len(row) == 1:
					try:
						data[i][key] = np.array(data[i][key])
					except NameError:
						# There's no key yet
						pass
					key = row[0]
					data[i][key] = []
				elif "segment positions" in row[0]:
					for item in row[1:]:
						try:
							value = float(item)
						except ValueError:
							segname = item
						else:
							zprofile[i][segname] = value
				elif len(row) == 3:
					try:
						xx_.append(float(row[0]))
						yy_.append(float(row[1]))
						rr_.append(float(row[2]))
					except ValueError:
						pass

		# When the last key is read, don't forget storing it
		data[i][key] = np.array(data[i][key])
		# wherenodes[i] = np.unique(np.array(wherenodes[i]))
		del key
		# axondata = np.array(axondata)
		# Store data in the dictionary
		i += 1

	return data, zprofile, (xx_, yy_, rr_)