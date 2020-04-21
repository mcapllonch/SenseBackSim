"""
Create folders for the new sets of results
"""
import os
import numpy as np
import shutil
import json



# List of numbers
values = np.around(np.arange(0.2, 4.2, 0.2), decimals=1)

# cwd
cwd = os.getcwd()

# Identify the folder that I want to replicate
mother = 'code'

# Get the contents of the file we want
with open(mother + '/EC/settings/electrodes.json', 'r') as f:
	info = json.load(f)

# Function to return the name of a folder given a number
fname = lambda x: 'current_%05inA'%(x*1000)

# Iterate over values
for v in values:
	name = fname(v)
	path = os.path.join(cwd, name)

	print('Creating folder: %s'%name)

	# Copy folder
	shutil.copytree(mother, name)

	# Modify information

	# First, copy the information
	info_here = {}
	info_here.update(info)
	info_here['left cuff']['stimulation protocol']['0, 0']['currents'] = [-v]
	
	# Dump it into the file
	with open(name + '/EC/settings/electrodes.json', 'w') as f:
		json.dump(info_here, f, indent=4)
