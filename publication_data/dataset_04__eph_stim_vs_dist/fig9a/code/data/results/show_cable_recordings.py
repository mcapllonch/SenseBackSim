"""
Show the recordings on cables (which can be either NAELC or cells)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import csv


# Read data
files = sorted(os.listdir(os.getcwd()))
for filename in files:
	if 'Axon' in filename or 'NAELC' in filename and '.csv' in filename:
		recs = []
		# print filename
		with open(filename, 'r') as f:
			frl = list(csv.reader(f))
			for row in frl:
				recs.append([float(x) for x in row])
		recs = np.array(recs).T

		# Coordinates:
		# Time and segment number
		time = 0.005 * np.arange(recs.shape[0])
		segments = np.arange(recs.shape[1])
		# Mesh them up
		# time, segments = np.meshgrid(time, segments)
		time, segments = np.meshgrid(segments, time)
		# print 'shapes:'
		# print time.shape, segments.shape, recs.shape

		# Show recordings
		fig, ax = plt.subplots()
		p = ax.pcolormesh(time, segments, recs, vmin=recs.min(), vmax=recs.max())
		cb = plt.colorbar(p, ax=ax)
		ax.set_title(filename)
		ax.set_xlabel('Segment number')
		ax.set_ylabel('time (ms)')
		newfilename = filename.replace('.csv', '.png')
		fig.savefig('./images/%s'%newfilename)
		plt.show()

		# Close all figures to save memory
		plt.close('all')