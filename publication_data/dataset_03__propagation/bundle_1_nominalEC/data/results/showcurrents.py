"""
Miguel Capllonch Juan
28 September 2018
"""

import numpy as np
import matplotlib.pyplot as plt
import csv




# Declare arrays
names = {}
ily1 = {}
ily2 = {}
balancely1 = {}
balancely2 = {}
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
		print 'A mem,', len(dataly1)
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

# Finish setting parameters that depend on the data
tarray = np.arange(0, dt * ndata, dt)
nt = len(tarray)

nsegstotal = len(z)

# Positions of the nodes of Ranvier
zRN = {}
# and indices corresponding to them
indsRN = {}
for k, v in names.items():
	zRN[k] = []
	indsRN[k] = []
	for i, vv in enumerate(v):
		if 'NODE' in vv:
			zRN[k].append(z[i])
			indsRN[k].append(i)

# Current balances
for ifib in ii_fibs:
	ily1[ifib] = np.array(ily1[ifib])
	ily2[ifib] = np.array(ily2[ifib])
	balancely1[ifib] = np.zeros(nt)
	balancely2[ifib] = np.zeros(nt)
	for i_t, t in enumerate(tarray):
		balancely1[ifib][i_t] = ily1[ifib][:, i_t].sum()
		balancely2[ifib][i_t] = ily2[ifib][:, i_t].sum()

###############################################################################
# Figures

# Time evolution at some point
fig, ax = plt.subplots(2, 1, sharex=True)
ax1, ax2 = ax.flatten()
i_t = nt / 2
ifib = ii_fibs[0]
iseg = indsRN[ifib][len(indsRN[ifib]) / 2 - 2]
isegRN = indsRN[ifib].index(iseg)
ax1.plot(tarray, ily1[ifib][iseg], label='Layer 1')
ax1.plot(tarray, ily2[ifib][iseg], label='Layer 2')
ax1.scatter(tarray[i_t], ily1[ifib][iseg][i_t], c='g', marker='*', s=100)
ax1.scatter(tarray[i_t], ily2[ifib][iseg][i_t], c='g', marker='*', s=100)
ax1.set_ylabel('i_mem (mA)')
ax1.set_title('Time evolution of i_mem at segment ' + names[ifib][iseg])
ax1.legend()
# Balance
ax2.plot(tarray, balancely1[ifib], label='Layer 1')
ax2.plot(tarray, balancely2[ifib], label='Layer 2')
ax2.scatter(tarray[i_t], balancely1[ifib][i_t], c='g', marker='*', s=100)
ax2.scatter(tarray[i_t], balancely2[ifib][i_t], c='g', marker='*', s=100)
ax2.set_xlim(tarray.min(), tarray.max())
ax2.set_xlabel('t (ms)')
ax2.set_ylabel('Current balance (mA)')
ax2.set_title('Time evolution of the current balance over the whole fiber %i'%ifib)
fig.tight_layout()
# plt.show()
fig.savefig('imem_vs_t.png')

# Along the axon at some point in time
fig, ax = plt.subplots()
y1 = ily1[ifib][i_t]
y2 = ily2[ifib][i_t]
dy = y1 - y2
# Vertical lines indicating the nodes or Ranvier
for zz in zRN[ifib]:
	ax.axvline(zz, c=(0.7, 0.7, 0.7))
# Data
print 'Shape check:', nsegstotal, len(z), len(y1)
ax.plot(z, [y1[i_z] for i_z in xrange(nsegstotal)], 
	label='Layer 1')
ax.plot(z, [y2[i_z] for i_z in xrange(nsegstotal)], 
	label='Layer 2')
ax.plot(z, [dy[i_z] for i_z in xrange(nsegstotal)], 
	label='Layer 1 - Layer 2')
ax.scatter(zRN[ifib], [y1[i_z] for i_z in indsRN[ifib]], c='r')
ax.scatter(zRN[ifib], [y2[i_z] for i_z in indsRN[ifib]], c='r')
ax.scatter(zRN[ifib], [dy[i_z] for i_z in indsRN[ifib]], c='r')
print len(indsRN), iseg
print zRN[ifib][isegRN], y1[iseg]
ax.scatter(zRN[ifib][isegRN], y1[iseg], c='g', marker='*', s=100)
ax.scatter(zRN[ifib][isegRN], y2[iseg], c='g', marker='*', s=100)
ax.scatter(zRN[ifib][isegRN], dy[iseg], c='g', marker='*', s=100)
# Arrange
ax.set_xlabel('z (um)')
ax.set_ylabel('i_mem (mA)')
ax.set_title('Profile of i_mem along the axon\nt = %0.03f ms'%tarray[i_t])
ax.set_xlim(min(z), max(z))
ax.legend()
# plt.show()
fig.savefig('imem_vs_z.png')