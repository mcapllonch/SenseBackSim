"""
Compare the AP trajectories of two simulations: SNOEC vs SEC.
Bundle 2.
"""

import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import scipy.signal as scs
import zipfile
import json
import csv
import os


# Time
dt = 0.005

# Cwd
cwd = os.getcwd()

# Items
cwd_items = os.listdir(cwd)

# List zip files
zipflist = [
	"bundle_2_noEC.zip", 
	"bundle_2_nominalEC.zip"
]

# Datasets where to store all data
datasets = OrderedDict()

# Method to select the maxima
method = "AP peaks"
method = "15 mV"

# Load style file
plt.style.use('./PaperFig_1textwidth_v3.mplstyle')

# Figure
fig = plt.figure()

grid = AxesGrid(fig, 111, 
			nrows_ncols=(1, len(zipflist)),
			direction="row",
			axes_pad=0.1,
			add_all=True,
			label_mode="L",
			share_all=True,
			cbar_location="right",
			cbar_mode="single",
			cbar_size="5%",
			cbar_pad=0.05,
			aspect=False
			)

# Arrange the figure
grid[0].set_ylabel("z (cm)")

# Select the titles
titles = {
	zipflist[0].replace(".zip", ""): "SNOEC", 
	zipflist[1].replace(".zip", ""): "SEC", 
}
# Iterate over zip files
for iz, zf in enumerate(zipflist):

	# Name of the dataset
	zname = zf.replace(".zip", "")

	# Dictionary to store data from this dataset
	datasets[zname] = {}

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

	# Axons' activity

	# Results folder
	results_folder = os.path.join(cwd, zfolder, "data/results")
	print(results_folder)
	# List all the items in the results folder
	items = sorted(os.listdir(results_folder))
	# Select the csv files
	items = [item for item in items if ".csv" in item]
	# Select the axon recording files
	items = [item for item in items if item[:4] == "Axon"]

	# Flags indicating which axons did fire and which not
	hasAP = []

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

	# Minima and maxima for the plot limits
	tmin = 0.
	tmax = -1e99
	zmin = 1e99
	zmax = -1e99

	# Iterate over the files and read them
	i = 0
	for filename in items:
		data[i] = {}
		zprofile_RN[i] = []
		with open(os.path.join(results_folder, filename), "r") as f:
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
		i += 1

	trajectories = {}

	for i, data_ in data.items():
		axondata = data_["v"]

		# Check if the maximum is an AP
		if axondata.max() > 15:

			# Time and position along the z-axis of the AP maxima
			ptime_ = []
			pzpos_ = []

			# Iterate over compartments
			for j, zj in enumerate(zprofile_RN[i]):

				if method == "AP peaks":
					# Peaks along time
					peak_time_indices = argrelextrema(axondata[j], np.greater)[0]
					peak_times = dt * peak_time_indices
					# Select only those which are greater than 15 mV
					local_maxima = axondata[j][peak_time_indices]
					peak_times = peak_times[np.where(local_maxima > 15)]
				elif method == "15 mV":
					# First time v reaches 15 mV
					tinds = np.where(axondata[j] >= 15)[0]
					if len(tinds) > 0:
						peak_time_indices = np.array([np.where(axondata[j] >= 15)[0].min()])
						peak_times = dt * peak_time_indices
						if j == 0:
							print('%s; v_max: %i; number of APs: %i'%(titles[zname], axondata.max(), len(peak_times)))
					else:
						peak_times = []
				# Add the data if there are valid peaks
				if len(peak_times) > 0:
					ptime_.append(peak_times)
					pzpos_.append(zj)
			data[i]["ptime_"] = np.array(ptime_)
			data[i]["pzpos_"] = np.array(pzpos_)

			# Minima and maxima
			tmax = max(tmax, max(ptime_))
			zmin = min(zmin, min(pzpos_))
			zmax = max(zmax, max(pzpos_))

			hasAP.append(True)
		else:
			hasAP.append(False)

	# Geometrical properties to array
	xx_ = np.array(xx_)
	yy_ = np.array(yy_)
	rr_ = np.array(rr_)

	# Colors for the diameters
	diameters_ = 2 * rr_
	dmin = diameters_.min()
	dmax = diameters_.max()
	diameters_scaled = (diameters_ - dmin) / diameters_.ptp()
	cmap_diams = plt.cm.gnuplot
	cmap_diams = plt.cm.jet
	diamcolors = [cmap_diams(x) for x in diameters_scaled]

	for i, data_ in data.items():
		if hasAP[i]:
			ptime_ = data_["ptime_"]
			pzpos_ = data_["pzpos_"]
			# Select the time of the first AP in each case
			times = np.array([t[0] for t in ptime_])

			# Plot
			grid[iz].plot(times, 1.e-4 * pzpos_, color=diamcolors[i], 
				markersize=1, lw=1, label="Axon %i"%i)

			# Store the dictionary that will be saved in a file
			datasets[zname][i] = {
				'times': times, 
				'pzpos_': pzpos_
			}

	# Arrange figure
	grid[iz].set_xlabel("Time (ms)")
	grid[iz].set_xticks(np.arange(tmin, 3., 0.5))
	grid[iz].set_xlim(tmin, tmax)
	grid[iz].set_ylim(zmin * 1.e-4, zmax * 1.e-4)
	grid[iz].set_title(titles[zname])

# Dummy scatters for the plots' colorbars
try:
	xdummy = times.mean() * np.ones_like(diameters_)
	ydummy = pzpos_.mean() * np.ones_like(diameters_)
except NameError:
	xdummy = np.zeros_like(diameters_)
	ydummy = np.zeros_like(diameters_)
scdummy = grid[iz].scatter(xdummy, ydummy, s=0, c=diameters_, cmap=cmap_diams, vmin=int(dmin), vmax=dmax)
# Colorbars
grid[-1].cax.colorbar(scdummy)
grid[-1].cax.axis["right"].toggle(ticklabels=True,label=True)
grid[-1].cax.set_ylabel(r"Fiber Diameter $(\rm \mu m)$")
grid[-1].cax.set_yticks(np.arange(3., 20., 3.))

plt.draw()
plt.show()

# Store data in a json file

# Custom json encoder
class NumpyEncoder(json.JSONEncoder):
	def default(self, obj):
		if isinstance(obj, np.ndarray):
			return obj.tolist()
		return json.JSONEncoder.default(self, obj)

# Write data to file
with open("fig8.json", "w") as f:
	json.dump(datasets, f, indent=4, cls=NumpyEncoder)