"""
Reduce the size of the data files
"""

import os
import re
import csv
import numpy as np



# Keys to be removed
removekeys = ["STIN", "MYSA", "FLUT"]
# removekeys = iter(removekeys)

# cwd
cwd = os.getcwd()

# List all the simulation folders
files = sorted(os.listdir(cwd))

for filename in files:

	if "Axon" in filename:

		# Find the file
		filename = os.path.join(cwd, filename)

		# File rows
		rows_ = []

		# Open and read it
		with open(filename, "r") as f:
			# Fetch data
			frl = list(csv.reader(f))
			for row in frl:
				skip = False
				# while not skip:
				# 	try:
				# 		skip = next(removekeys) in row[0]
				# 	except StopIteration:
				# 		skip = 
				# 	if next(removekeys) in row[0]:
				# 		skip = True
				# if not skip:
				# 	rows_.append(row[:])
				for rk in removekeys:
					if rk in row[0]:
						skip = True
						break
				if not skip:
					rows_.append(row[:])


		# Finally, overwrite the old file with just the necessary data 
		# to save disk space
		with open(filename, "w") as f:
			fw = csv.writer(f)
			for row in rows_:
				fw.writerow(row)
