"""
Miguel Capllonch Juan
29 September 2018
This computes the LFP from membrane currents in an anisotropic medium
"""

import numpy as np
import matplotlib.pyplot as plt
import csv



# Recording electrodes
rec_els = {
	'E': {'pos': (250, 0, 19000), 
		'lfp_indiv_ly1': {}, 'lfp_indiv_ly2': {}, 
		'lfp_total_ly1': None, 'lfp_total_ly2': None, 
		'color': 'r'},
	'S': {'pos': (0, -250, 19000), 
		'lfp_indiv_ly1': {}, 'lfp_indiv_ly2': {}, 
		'lfp_total_ly1': None, 'lfp_total_ly2': None, 
		'color': 'g'},
	'W': {'pos': (-250, 0, 19000), 
		'lfp_indiv_ly1': {}, 'lfp_indiv_ly2': {}, 
		'lfp_total_ly1': None, 'lfp_total_ly2': None, 
		'color': 'b'},
	'N': {'pos': (0, 250, 19000), 
		'lfp_indiv_ly1': {}, 'lfp_indiv_ly2': {}, 
		'lfp_total_ly1': None, 'lfp_total_ly2': None, 
		'color': 'k'}
}

# Stimulating electrode (amp in mA)
delay = 0.1
dur = 0.2
amp = -15.e-3
stimcurrent = None
stimelpos = (250, 0, 100)

# Medium conductivity tensor (1/(Ohm*um))
sigma_x = 1. / (1.211e3 * 1.e-6)
sigma_y = 1. / (1.211e3 * 1.e-6)
sigma_z = 1. / (0.175e3 * 1.e-6)

# Functions
def compute_lfp(currents, elpos):
	""" Compute the LFP from a time series of currents as recorded by 
	a point electrode situated at elpos
	The equation is taken from:
	Nicholson & Freeman (1975) 
	The current sources are all the segments in the axon """
	dx = x - elpos[0]
	dy = y - elpos[1]
	dz = z - elpos[2]
	denominator = 4. * np.pi * np.sqrt\
		(\
			sigma_y * sigma_z * dx ** 2 + \
			sigma_z * sigma_x * dy ** 2 + \
			sigma_x * sigma_y * dz ** 2\
		)
	# denominator = np.tile(denominator, nt).reshape(currents.shape)
	denominator = denominator.repeat(nt).reshape(currents.shape)
	# print dz.shape, (dz ** 2).shape, currents.shape
	return (currents / denominator).sum(axis=0)

def compute_lfp_fromtimeseries(currents, srcpos, elpos):
	""" Compute the LFP from a time series of currents as recorded by 
	a point electrode situated at elpos
	The equation is taken from:
	Nicholson & Freeman (1975) 
	This time, there is only one current point source """
	dx = srcpos[0] - elpos[0]
	dy = srcpos[1] - elpos[1]
	dz = srcpos[2] - elpos[2]
	denominator = 4. * np.pi * np.sqrt\
		(\
			sigma_y * sigma_z * dx ** 2 + \
			sigma_z * sigma_x * dy ** 2 + \
			sigma_x * sigma_y * dz ** 2\
		)
	# denominator = np.tile(denominator, nt).reshape(currents.shape)
	# denominator = denominator.repeat(nt).reshape(currents.shape)
	# denominator = denominator.repeat(nt)
	# print dz.shape, (dz ** 2).shape, currents.shape
	return currents / denominator

# Declare arrays
names = {}
ily1 = {}
ily2 = {}
balancely1 = {}
balancely2 = {}
lfp_indiv_ly1 = {}
lfp_indiv_ly2 = {}
x = []
y = []
z = []
ii_fibs = []

# Other parameters
dt = 0.005

# Get data
with open('./membranecurrents.csv', 'r') as f:
	frl = list(csv.reader(f))
	for i, row in enumerate(frl[1:]):
		ifib = int(row[0])
		ii_fibs.append(ifib)
		name = row[1]
		data = row[5:]
		ndata = len(data) / 2
		dataly1 = np.array([float(item) for item in data[:ndata]])
		dataly2 = np.array([float(item) for item in data[ndata:]])
		try:
			ily1[ifib].append(dataly1.copy())
			ily2[ifib].append(dataly2.copy())
			names[ifib].append(name)
		except KeyError:
			names[ifib] = [name]
			ily1[ifib] = [dataly1]
			ily2[ifib] = [dataly2]
		x.append(float(row[2]))
		y.append(float(row[3]))
		z.append(float(row[4]))

# Positions from lists to arrays
x = np.array(x)
y = np.array(y)
z = np.array(z)

# Finish setting parameters that depend on the data
tarray = np.arange(0, dt * ndata, dt)
nt = len(tarray)
nsegstotal = len(z)
stimcurrent = amp * np.ones_like(tarray)
# stimcurrent[np.where(tarray < delay + dur)] = amp
stimcurrent[np.where(tarray < delay)] = 0.
stimcurrent[np.where(tarray > delay + dur)] = 0.
# stimcurrent = stimcurrent()
# stimcurrent[np.where(delay < tarray < delay + dur)] = amp

# Positions of the nodes of Ranvier
zRN = {}
# and indices corresponding to them
indsRN = {}
for k, v in names.items():
	zRN[k] = []
	indsRN[k] = []
	for i, vv in enumerate(v):
		if 'node' in vv:
			zRN[k].append(z[i])
			indsRN[k].append(i)

for ifib in ii_fibs:

	# Current balances
	ily1[ifib] = np.array(ily1[ifib])
	ily2[ifib] = np.array(ily2[ifib])
	balancely1[ifib] = np.zeros(nt)
	balancely2[ifib] = np.zeros(nt)
	for i_t, t in enumerate(tarray):
		balancely1[ifib][i_t] = ily1[ifib][:, i_t].sum()
		balancely2[ifib][i_t] = ily2[ifib][:, i_t].sum()

	# Individual LFPs
	for k, re in rec_els.items():
		re['lfp_indiv_ly1'][ifib] = compute_lfp(ily1[ifib], re['pos'])
		re['lfp_indiv_ly2'][ifib] = compute_lfp(ily2[ifib], re['pos'])

# Finally, sum up the individual LFPs of the fibers into a total LFP
# for each electrode
for k, re in rec_els.items():
	re['lfp_total_ly1'] = np.zeros(nt)
	re['lfp_total_ly2'] = np.zeros(nt)
	for ifib in ii_fibs:
		re['lfp_total_ly1'] += re['lfp_indiv_ly1'][ifib]
		re['lfp_total_ly2'] += re['lfp_indiv_ly2'][ifib]
	# Add the contribution of the stimulating electrode
	re['lfp_total_ly1'] += compute_lfp_fromtimeseries(stimcurrent, 
		stimelpos, re['pos'])
	re['lfp_total_ly2'] += compute_lfp_fromtimeseries(stimcurrent, 
		stimelpos, re['pos'])

# What if I sum them?
resum = re['lfp_total_ly1'] + re['lfp_total_ly2']

###############################################################################
# Figures

# Time evolution at some point
fig, ax = plt.subplots()
for k, re in rec_els.items():
	ax.plot(tarray, re['lfp_total_ly1'], c=re['color'] , ls='-', 
		label=k + '. Layer 1')
	ax.plot(tarray, re['lfp_total_ly2'], c=re['color'] , ls='--', 
		label=k + '. Layer 2')
ax.plot(tarray, resum, 'r', lw=3, label='Sum')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Extracellular recordings (mV)')
ax.set_title('Extracellular recordings')
ax.legend()
fig.tight_layout()
# plt.show()
fig.savefig('recordings_VC.png')