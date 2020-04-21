"""
Read the results in csv from the simulation outputs and save the fields in xyz format
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import read_results



def save_vext_xyzt(data, name):
	""" Save all the data for vext in a xyzt file
	NOTE: This is a highly inefficient method. It also uses a lot of memory and creates large files """

	# Save results in xyzt format
	with open('results_vext_%s.xyzt'%name, 'w') as f:

		# Iterate over cables
		for i in all_cable_indices:
			print('Cable %i'%i)

			# Choose data
			try:
				data_ = data[i]['vext[0]']
			except KeyError:
				data_ = data[i]['vext[1]']

			# Iterate over z positions
			for iz, z in enumerate(zprofile[i].values()):

				# Iterate over time
				for t, v in zip(time, data_[iz]):

					f.write('%E, %E, %E, %E, %E\n'%(xx_[i], yy_[i], z, v, t))

def save_vext_xyz(data, name, t=None, it=None, tarray=None):
	""" Save all the data for vext for a certain time step in a xyz file """

	# Choose it from the input
	if it is not None:
		if t is not None:
			print('ERROR: In save_vext_xyz. Please choose either t or it, not both')
			sys.exit()
	elif t is not None:
		# A value of t has been chosen
		if tarray is None:
			print('ERROR: In save_vext_xyz. If a value of t is provided, please provide also the time array (tarray)')
			sys.exit()
		else:
			it = np.where(tarray == t)[0][0]

	# Save results in xyz format
	with open('results_vext_%s_time_step_%05i.xyz'%(name, it), 'w') as f:

		# Iterate over cables
		for i in all_cable_indices:
			print('Cable %i'%i)

			# Choose data
			try:
				data_ = data[i]['vext[0]']
			except KeyError:
				data_ = data[i]['vext[1]']

			# Iterate over z positions
			for z, v in zip(zprofile[i].values(), data_[:, it]):

				f.write('%E, %E, %E, %E\n'%(xx_[i], yy_[i], z, v))

def save_v_xyz(data, name, t=None, it=None, tarray=None):
	""" Save all the data for v for a certain time step in a xyz file """

	# Choose it from the input
	if it is not None:
		if t is not None:
			print('ERROR: In save_v_xyz. Please choose either t or it, not both')
			sys.exit()
	elif t is not None:
		# A value of t has been chosen
		if tarray is None:
			print('ERROR: In save_v_xyz. If a value of t is provided, please provide also the time array (tarray)')
			sys.exit()
		else:
			it = np.where(tarray == t)[0][0]

	# Save results in xyz format
	with open('results_v_%s_time_step_%05i.xyz'%(name, it), 'w') as f:

		# Iterate over cables
		for i in all_cable_indices:
			print('Cable %i'%i)

			# Choose data
			try:
				data_ = data[i]['v']
			except KeyError:
				pass
			else:

				# Iterate over z positions
				for z, v in zip(zprofile[i].values(), data_[:, it]):

					f.write('%E, %E, %E, %E\n'%(xx_[i], yy_[i], z, v))

# Results folder
folder = os.path.join(os.getcwd(), "code/data/results")

# Read the results
data, zprofile, xyr = read_results.read_results(folder)
print('Finished with read_results')
xx_, yy_, rr_ = xyr

# All cable indices
all_cable_indices = list(data.keys())

# Variable of interest
var_interest = 'v'
var_interest = 'vext[1]'

# Callables
callables = {
		'v': save_v_xyz, 
		'vext[0]': save_vext_xyz, 
		'vext[1]': save_vext_xyz, 
	}
# Time
nt = data[0][var_interest].shape[1]
dt = 0.005
tarray = np.arange(0, nt * dt, dt)

# Name of the simulation
name = 'dataset01'

# Save data in a xyz file for t = 1.
callables[var_interest](data, name, t=0.015, tarray=tarray)