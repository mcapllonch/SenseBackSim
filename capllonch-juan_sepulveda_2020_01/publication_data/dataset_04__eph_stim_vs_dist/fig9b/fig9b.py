"""
Plot the ephaptic influence vs. inter-axonal distance
"""

import numpy as np
from scipy.signal import argrelextrema
from matplotlib.gridspec import GridSpec
from collections import OrderedDict
import matplotlib.pyplot as plt
import json
import csv
import sys
import os


# Load style file
plt.style.use('../PaperFig_0.45_textwidth.mplstyle')

# Results folder
cwd = os.getcwd()
results_folder = os.path.join(cwd, "code/data/results")
# List all the items in the results folder
items = sorted(os.listdir(results_folder))
# Select the csv files
items = [item for item in items if ".csv" in item]
# Select the axon recording files
items = [item for item in items if item[:4] == "Axon"]

# Array for the axons' activity maxima and minima
maxima = {}
minima = {}
vminimum = 1e99

# Flags indicating which axons did fire and which not
hasAP = {}

# Voltage data
data = {}

# Geometrical properties
xx_ = []
yy_ = []
rr_ = []

# z-axis anatomy
# Nodes of Ranvier
zprofile_RN = {}

# Diamters
diameters = OrderedDict()

# Iterate over the files and read them
i = 0
for filename in items:
	data[i] = {}
	zprofile_RN[i] = []
	with open(os.path.join(results_folder, filename), "r") as f:
		# data[i][key] = []
		fr = csv.reader(f)
		for row in fr:
			if "NODE" in row[0]:
				data[i][key].append([float(x) for x in row[1:]])
			elif len(row) == 3:
				xx_.append(float(row[0]))
				yy_.append(float(row[1]))
				rr_.append(float(row[2]))
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
						if "NODE" in segname:
							zprofile_RN[i].append(value)
				# Store zprofile_RN[i] as an array
				zprofile_RN[i] = np.array(zprofile_RN[i])
	# When the last key is read, don't forget storing it
	data[i][key] = np.array(data[i][key])
	del key
	# Store data in the dictionary
	i += 1

# Geometrical properties to array
xx_ = np.array(xx_)
yy_ = np.array(yy_)
rr_ = np.array(rr_)
# Distances of the axons from the center of the bundle
dists_center = np.hypot(xx_, yy_)
# Scale it
dists_center_scaled = dists_center / dists_center.max()

vrest = -80

# Avoid the margins to avoid boundary effects
nodes = np.arange(len(zprofile_RN[0]))[2:-2]

# deflections will be an array, but now it needs to be a list
deflections = []

for node in nodes:
	maxima[node] = []
	minima[node] = []
	for i, data_ in data.items():
		axondata = data_["v"]

		# Store maximum and minimum value of v
		maxima[node].append(axondata[node].max())
		minima[node].append(axondata[node].min())
		# Store overall minimum
		vminimum = min(vminimum, axondata.min())

		# Check if the maximum is an AP
		if axondata.max() > 15:
			hasAP[i] = True
			stimulated_axon = i
		else:
			try:
				hasAP[i]
			except KeyError:
				hasAP[i] = False

	# Deflections from -80 mV
	maxima[node] = np.array(maxima[node])
	minima[node] = np.array(minima[node])
	deflections_item = maxima[node] - vrest

	# Remove the highest value, which corresponds to the source axon
	deflections_item = np.ma.masked_where(deflections_item == deflections_item.max(), deflections_item)
	# Append deflections_item to the list deflections
	deflections.append(deflections_item)

# deflections to array
deflections = np.ma.array(deflections)

# Mean of deflections per node
defmeanpernode = deflections.mean(axis=1)

# Max of deflections per axon
defmaxperaxon = deflections.max(axis=0)

# Figure
plt.close('all')
fig, ax = plt.subplots()

# Distances from the stimulated axon
distances = np.hypot(xx_ - xx_[stimulated_axon], yy_ - yy_[stimulated_axon])

ax.scatter(distances, defmaxperaxon, c="k", s=1)

ax.set_xlabel(r"Distance ($\rm \mu m$)")
ax.set_ylabel(r"Rise above $v_{r}$ $\rm (mV)$")

ax.set_xticks(np.arange(0, 500, 100))

plt.show()

# Save results in a csv
with open('fig9b.csv', 'w') as f:
	fw = csv.writer(f)
	fw.writerow(['distance (um)', 'rise (mv)'])
	for d, r in zip(distances, defmaxperaxon):
		fw.writerow([d, r])