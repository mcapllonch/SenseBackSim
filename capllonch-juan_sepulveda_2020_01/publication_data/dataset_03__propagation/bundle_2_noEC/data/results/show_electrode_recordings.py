import os
import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.signal as scs




# With or without artefact
withart = 'no'
withart = 'with'

# Read data
files = sorted(os.listdir(os.getcwd()))
for filename in files:
	time = []
	recs = []
	if 'recordings' in filename and '.csv' in filename:
		print filename
		with open(filename, 'r') as f:
			frl = list(csv.reader(f))
			for row in frl[1:]:
				time.append(float(row[0]))
				recs.append(float(row[1]))

		# Show recordings
		fig, ax = plt.subplots()
		ax.plot(time, recs)
		ax.set_title('Electrode recordings')
		ax.set_xlabel('time (ms)')
		ax.set_ylabel('Vext (mV)')
		newfilename = filename.replace('.csv', '.png')
		fig.savefig('./images/%s'%newfilename)

		# Compute PSD
		dt = 0.005 # ms
		fs = 1. / (dt * 1.e-3) # Hz
		f, Pxx_den = scs.welch(recs, fs)

		# Save PSD in file

		# Show PSD
		fig, ax = plt.subplots()
		ax.semilogy(f, Pxx_den)
		# plt.ylim(0.5e-3)
		ax.set_title('PSD of electrode recordings')
		ax.set_xlabel('frequency (Hz)')
		ax.set_ylabel('PSD')
		fig.savefig('./images/PSD%s'%newfilename)

		# Show PSD only until 10 kHz
		xmax = 10000
		window = np.where(f <= xmax)
		fig, ax = plt.subplots()
		ax.semilogy(f[window], Pxx_den[window])
		# plt.ylim(0.5e-3)
		ax.set_xlim(0, xmax)
		ax.set_title('PSD of electrode recordings')
		ax.set_xlabel('frequency (Hz)')
		ax.set_ylabel('PSD')
		fig.savefig('./images/PSD_until_10kHz%s'%newfilename)
		# Close all figures to save memory
		plt.close('all')