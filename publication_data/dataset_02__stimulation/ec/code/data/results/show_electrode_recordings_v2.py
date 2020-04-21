import os
import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.signal as scs




# With or without artefact
withart = 'no'
withart = 'with'

# Dictionary for all the recordings time series
recordings = {}

# Colors for the plots
colors = ["b", "g", "r", "purple", "k", "brown"]

# Maxima of rings and pads
maxrings = 0
maxpads = 0

# Choose if to visualise only one pad or not
only1pad = True


# Choose the time limit on the left
tmin = 0.21

# Read data
files = sorted(os.listdir(os.getcwd()))
for filename in files:
	time = []
	recs = []
	if 'recordings' in filename and '.csv' in filename:
		print(filename)
		with open(filename, 'r') as f:
			frl = list(csv.reader(f))
			for row in frl[1:]:
				time.append(float(row[0]))
				recs.append(float(row[1]))
		key = filename.replace('.csv', '').replace('recordings_E_', '')
		key = key.replace('_noartefacts', '')
		key = key.replace('_withartefacts', '')
		key = key.split('_')[-1]
		key = key.replace('R', '')
		ring, pad = key.split('P')
		ring = int(ring)
		pad = int(pad)
		data = {
			"time": time, 
			"recs": recs, 
			"color": colors[ring]
		}
		recordings[(ring, pad)] = data
		# print(ring, pad)
		maxrings = max(maxrings, ring)
		maxpads = max(maxpads, pad)

# Add 1 to the maxima to account for the Python indexing convention
maxrings += 1
maxpads += 1

# Force to plot only for one pad of each ring
if only1pad:
	maxpads = 1

# Create the figure and plot stuff in it
fig, ax = plt.subplots(maxpads, 2, figsize=(10, 5*maxpads))

# for (ring, pad), data in recordings.items():
for ring in range(maxrings):
	for pad in range(maxpads):

		# This conditional is here only in case I want to show the 
		# recordings until a certain pad
		if pad < maxpads:

			s0 = np.s_[pad, 0]
			s1 = np.s_[pad, 1]
			if maxpads == 1:
				s0 = np.s_[0]
				s1 = np.s_[1]

			# Gather data
			data = recordings[(ring, pad)]
			time = data["time"]
			recs = data["recs"]
			color = data["color"]

			# Show recordings
			
			ax[s0].plot(time, recs, color=color)
			ax[s0].set_ylabel('Vext (mV)')
			ax[s0].set_xlim(left=tmin)

			# Compute PSD
			dt = 0.005 # ms
			fs = 1. / (dt * 1.e-3) # Hz
			f, Pxx_den = scs.welch(recs, fs)

			# Show PSD only until 10 kHz
			xmax = 10000
			window = np.where(f <= xmax)
			
			ax[s1].semilogy(f[window], Pxx_den[window], label="Ring %i"%ring, color=color)
			# plt.ylim(0.5e-3)
			ax[s1].set_xlim(0, xmax)
			ax[s1].set_ylabel('PSD')

for i in range(maxpads):
	s0 = np.s_[i, 0]
	s1 = np.s_[i, 1]
	if maxpads == 1:
		s0 = np.s_[0]
		s1 = np.s_[1]
	ax[s0].set_ylim(-7, 7)
if maxpads > 1:
	ax[0, 0].set_title('Electrode recordings')
	ax[maxpads-1, 0].set_xlabel('Time (ms)')
	ax[0, 1].set_title('PSD of electrode recordings')
	ax[maxpads-1, 1].set_xlabel('frequency (Hz)')
	ax[0, 1].legend(fontsize=10)
else:
	ax[0].set_title('Electrode recordings')
	ax[0].set_xlabel('Time (ms)')
	ax[1].set_title('PSD of electrode recordings')
	ax[1].set_xlabel('frequency (Hz)')
	ax[1].legend(fontsize=10)

fig.tight_layout()
plt.show()
fig.savefig('./images/%s'%'recordings_%sartefact.png'%withart)