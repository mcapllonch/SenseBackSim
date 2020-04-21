"""
Show vm vs t of all the axons for a certain node
This is for simulations where I stimulated only one axon.
I want to show the (maximum) influence on vm that its AP has on other 
axons, and plot this magnitude vs. the distance to this axon
IMPORTANT: All the axons must have the same diameter and be perfectly 
aligned
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
# plt.style.use('../../../../PaperFig_0.45_textwidth.mplstyle')

# Results folder
cwd = os.getcwd()
results_folder = os.path.join(cwd, "data/results")
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

nodes = np.arange(len(zprofile_RN[0]))
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
		# print('Axon %i. minimum[node%i] + 80 mV = %f mV'%(i, node, minima[-1] + 80.))
		# print(i, axondata.max(), axondata.shape)

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
	# print("Total Minimum: %f"%vminimum)
	maxima[node] = np.array(maxima[node])
	minima[node] = np.array(minima[node])
	if True:
		fn_append = ''
		deflections_item = maxima[node] - vrest
	else:
		fn_append = '_v2'
		deflections_item = maxima[node] - minima[node]
	# Remove the highest value
	deflections_item = np.ma.masked_where(deflections_item == deflections_item.max(), deflections_item)
	# Append deflections_item to the list deflections
	deflections.append(deflections_item)

# deflections to array
deflections = np.ma.array(deflections)

# Mean of deflections per node
defmeanpernode = deflections.mean(axis=1)

# Time domain size
nt = len(axondata[node])
# Time array
dt = 0.005
tarray = np.arange(0., nt * dt, dt)

# Figure
plt.close('all')
fig, ax = plt.subplots()

# Distances from the stimulated axon
distances = np.hypot(xx_ - xx_[stimulated_axon], yy_ - yy_[stimulated_axon])
# Scale it
distances_scaled = distances / distances.max()

# Choose one node of Ranvier to sample
node_ = 1

for i, data_ in data.items():
	if not hasAP[i]:
		# ax.plot(nodes, deflections[:, i], c=plt.cm.jet(dists_center_scaled[i]))
		# ax.plot(nodes, deflections[:, i], c=plt.cm.jet(distances_scaled[i]))
		# ax.plot(tarray, data_["v"][node_], c=plt.cm.jet(distances_scaled[i]))
		# ax.plot(tarray, data_["v"][node_], c=plt.cm.jet(dists_center_scaled[i]))
		ax.plot(tarray, data_["v"][node_], c=plt.cm.jet(dists_center_scaled[i]))
		# ax.plot(nodes, deflections[:, i] - defmeanpernode, c=plt.cm.jet(dists_center_scaled[i]))
		# ax.plot(nodes, deflections[:, i] - defmeanpernode, c=plt.cm.jet(distances_scaled[i]))
# sc = ax.scatter(xx_, yy_, c=deflections, s=rr_)
# cb = plt.colorbar(sc, ax=ax)

# Stimulated axon
ax.plot(tarray, data[stimulated_axon]["v"][node_], c='k')

# Arrange the figure
# ax.set_xlabel(r"z ($\rm \mu m$)"))
ax.set_title('deflections minus mean')
# ax.set_yscale('log')
ax.set_xlabel(r"Time (ms)")
ax.set_ylabel(r"Rise above $v_{r}$ $\rm (mV)$")
# ax.set_xlim(0.1, 0.2)
# ax.set_ylim(-80.06, -80.04)

# ax.set_ylabel(r"Rise above $v_{r}$ $\rm (mV)$")
# ax.set_xticks(np.arange(0, 500, 100))
# ax.set_xticks(np.arange(0, 250, 50))
# fig.tight_layout()

model = 5
number = 1
if False:
	fig.savefig("../../../../figures_for_publication/eph_stim_vs_distance_model%i_node%03i_small.pdf"%(model, node_), dpi=300)
	fig.savefig("../../../../figures_for_publication/eph_stim_vs_distance_model%i_node%03i_small.png"%(model, node_), dpi=300)
	fig.savefig("../../../../figures_for_publication/eph_stim_vs_distance_model%i_node%03i_small.tiff"%(model, node_), dpi=300)
fig.savefig("vm_vs_t_model%i_node%03i_small%s_%s.png"%(model, node_, fn_append, number), dpi=300)
plt.show()