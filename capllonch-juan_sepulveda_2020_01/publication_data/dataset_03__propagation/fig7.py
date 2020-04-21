"""
Show the CVs of "EC" scaled with the CVs of "No EC"
"""

import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import scipy.signal as scs
import scipy.stats as sst
import zipfile
import json
import csv
import os


def get_slope(x, y):
	"""
	Obtain the slope of the linear regression of a number of points
	"""
	return sst.linregress(x, y)[0]

def draw_subfigure(ax, zipflist):
	""" Draw a subfigure """

	# Method to select the maxima
	method = "AP peaks"
	method = "15 mV"

	# Velocities
	velocities = {}

	# Iterate over zip files
	for iz, zf in enumerate(zipflist):

		# Name of the dataset
		zname = zf.replace(".zip", "")

		# Dictionary to store data from this dataset
		velocities[zname] = {}

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
		diameters_scaled = (diameters_ - 9.) / 2.
		cmap_diams = plt.cm.gnuplot
		cmap_diams = plt.cm.jet
		diamcolors = [cmap_diams(x) for x in diameters_scaled]

		for i, data_ in data.items():
			if hasAP[i]:
				ptime_ = data_["ptime_"]
				pzpos_ = data_["pzpos_"]
				# Select the time of the first AP in each case
				times = np.array([t[0] for t in ptime_])

				# Minima and maxima
				tmax = max(tmax, times.max())

				# Conduction Velocities

				# Avoid boundaries
				margin = 3
				times_cropped = times[margin:-margin]
				pzpos_cropped = pzpos_[margin:-margin]

				# Velocities: LINEAR REGRESSIONS
				vels = []

				# Window: ALWAYS ODD!
				window_lr = 11
				lr_margin = window_lr // 2

				for j in range(len(times_cropped) - window_lr + 1):
					vels.append(
						get_slope(
								times_cropped[j : j + window_lr], 
								pzpos_cropped[j : j + window_lr]
							)
						)
				# Velocities to array
				vels = np.array(vels)

				# Save data
				velocities[zname][i] = {
					"time": times_cropped[lr_margin:-lr_margin], 
					"velocities": vels
				}

	# Dummy scatter for the plots' colorbar
	try:
		xdummy = times.mean() * np.ones_like(diameters_)
		ydummy = pzpos_.mean() * np.ones_like(diameters_)
	except NameError:
		xdummy = np.zeros_like(diameters_)
		ydummy = np.zeros_like(diameters_)
	scdummy = ax.scatter(xdummy, ydummy, s=0, c=diameters_, cmap=cmap_diams, vmin=int(dmin), vmax=11)

	# Averages for noEC
	avgnoEC = np.zeros_like(diameters_)
	for zname, data_ in velocities.items():
		if "noEC" in zname:
			for i, vels in data_.items():
				vels = vels["velocities"]
				avgnoEC[i] = vels.mean()

	# Plot velocities
	for zname, data_ in velocities.items():
		if "nominalEC" in zname:
			for i, vels in data_.items():
				time = vels["time"]
				vels = vels["velocities"]
				# Scale the velocities and plot
				ax.plot(time, vels / avgnoEC[i], color=diamcolors[i])

	# Return scdummy to use it for the colorbar
	return velocities, scdummy

# Time
dt = 0.005

# Cwd
cwd = os.getcwd()

# Items
cwd_items = os.listdir(cwd)

# Load style file
plt.style.use('./PaperFig_0.8_textwidth_v6.mplstyle')

# Figure
plt.close('all')

fig, ax = plt.subplots()

# List zip files containing the outputs of the simulations
zipflist = [
	"bundle_1_noEC.zip", 
	"bundle_1_nominalEC.zip"
]

# Draw the subfigures
velocities, scdummy = draw_subfigure(ax, zipflist)

# Arrange the figure
ax.set_xlim(0.5, 3.0)
ax.set_xlabel("Time (ms)")
ax.set_ylim(0.3, 1.)
ax.set_ylabel("Scaled CV")

# Colorbar
ticks = np.arange(9., 11., 0.3)
cb = plt.colorbar(scdummy, ax=ax, label=r"Fiber Diameter ($\rm \mu m$)", ticks=ticks)

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
with open("fig7.json", "w") as f:
	json.dump(velocities, f, indent=4, cls=NumpyEncoder)