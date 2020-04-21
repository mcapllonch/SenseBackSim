"""
Create folders for the new sets of results
Note: integer values are not created because I already have these results in previous sets
"""
import os
import numpy as np
import shutil
import json



# List of numbers
values = [x for x in np.around(np.arange(0.2, 4.2, 0.2), decimals=1) if not x.is_integer()]

# cwd
cwd = os.getcwd()

# Identify the folder that I want to replicate
mother = 'current_00500nA'

# Get the contents of the file we want
with open(mother + '/noEC/settings/electrodes.json', 'r') as f:
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
	with open(name + '/noEC/settings/electrodes.json', 'w') as f:
		json.dump(info_here, f, indent=4)
