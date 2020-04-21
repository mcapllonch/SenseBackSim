"""
Draw the stimulation results
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict


keys = [
	'nerve', 
	'Fascicle_00', 
	'Fascicle_01', 
	'Fascicle_02', 
	'Fascicle_03', 
	'Fascicle_04', 
	'Fascicle_05', 
	'Fascicle_06', 
]

sorted_keys = [
	'nerve', 
	'Fascicle_01', 
	'Fascicle_02', 
	'Fascicle_03', 
	'Fascicle_04', 
	'Fascicle_05', 
	'Fascicle_06', 
	'Fascicle_00', 
]

nkeys = len(keys)

ec_keys = [
	'EC', 
	'noEC', 
]

names = {
	'nerve': 'Nerve',
	'Fascicle_00': 'Fascicle 7',
	'Fascicle_01': 'Fascicle 1',
	'Fascicle_02': 'Fascicle 2',
	'Fascicle_03': 'Fascicle 3',
	'Fascicle_04': 'Fascicle 4',
	'Fascicle_05': 'Fascicle 5',
	'Fascicle_06': 'Fascicle 6',
}

naxons = {
	'nerve': 658,
	'Fascicle_00': 91,
	'Fascicle_01': 82,
	'Fascicle_02': 118,
	'Fascicle_03': 99,
	'Fascicle_04': 87,
	'Fascicle_05': 98,
	'Fascicle_06': 83,
}

colors = {
	'nerve': 'k', 
	'Fascicle_00': 'r', 
	'Fascicle_01': 'g', 
	'Fascicle_02': 'b', 
	'Fascicle_03': 'yellow', 
	'Fascicle_04': 'orange', 
	'Fascicle_05': 'purple', 
	'Fascicle_06': 'cyan'
}

markers = {
	'nerve': 'o', 
	'Fascicle_00': 'D', 
	'Fascicle_01': 's', 
	'Fascicle_02': '^', 
	'Fascicle_03': 'v', 
	'Fascicle_04': '<', 
	'Fascicle_05': '>', 
	'Fascicle_06': 'x'
}

styles = {
	'noEC': '--', 
	'EC': '-', 
}

simnames = {
	'noEC': 'SNOEC', 
	'EC': 'SEC', 
}

currents = []
recruitment = {}

# Initialize recruitment data
for k in keys:
	recruitment[k] = {
		'noEC': {
			'absolute': [], 
			'fractional': [], 
			}, 
		'EC': {
			'absolute': [], 
			'fractional': [], 
			}, 
	}

# Read data
with open('recruitment_data.csv', 'r') as f:
	frl = list(csv.reader(f))
	for row in frl[2:]:
		current = row[0]
		currents.append(np.abs(float(current)))
		for i, x in enumerate(row[1:]):
			thiskey = keys[i % nkeys]
			ec_key = ec_keys[::-1][i // nkeys]
			try:
				recruitment[thiskey][ec_key]['absolute'].append(int(x))
			except ValueError:
				recruitment[thiskey][ec_key]['absolute'].append(np.nan)
				recruitment[thiskey][ec_key]['fractional'].append(np.nan)
			else:
				recruitment[thiskey][ec_key]['fractional'].append(int(x) / naxons[thiskey])

# Once finished reading, turn every list to array and plot it against the currents
currents = np.array(currents)


##########################################################################
# Figures

# Absolute recruitment to arrays
for i, k in enumerate(keys):
	for ec_key in ec_keys:
		recruitment[k][ec_key]['absolute'] = np.array(recruitment[k][ec_key]['absolute'])

# Fractional recruitment
plt.style.use('./PaperFig_0.8_textwidth_v4.mplstyle')
fig2, ax2 = plt.subplots(4, 2)
ax2 = ax2.flatten()
linecolors = {'noEC': 'b', 'EC': 'k'}
diffs = OrderedDict()
for i, k in enumerate(sorted_keys):
	for ec_key in ec_keys:
		recruitment[k][ec_key]['fractional'] = np.array(recruitment[k][ec_key]['fractional'])
		try:
			ax2[i].plot(currents, recruitment[k][ec_key]['fractional'], c=linecolors[ec_key], ls='-', label=simnames[ec_key])
		except KeyError:
			ax2[i].plot(currents, recruitment[k][ec_key]['fractional'], c=linecolors[ec_key], ls='-')
		ax2[i].text(2.5, 0.1, names[k])

	# Differences
	diff = recruitment[k]['EC']['fractional'] - recruitment[k]['noEC']['fractional']
	diffs[k] = diff
	# ax2[i].plot(currents, diff, c='k', ls=':')
	# ax2[i].plot(currents, diff, c='k', ls='--')
	ax2[i].plot(currents, diff, c='r', ls='-', label='Difference')

	# Arrange figure
	ax2[i].set_yticks([])
	if i < 6:
		ax2[i].set_xticklabels([])
	if i % 2 != 0:
		ax2[i].set_yticklabels([])
	else:
		# ax2[i].set_yticks([0, 0.5, 1])
		# ax2[i].set_yticklabels([0, 0.5, 1])
		ax2[i].set_yticks([0, 1])
		ax2[i].set_yticklabels([0, 1])
		# ax2[i].set_title(k)
	ax2[i].set_xlim(0, 4)
	ax2[i].set_ylim(0, 1)
fig2.text(0.5, 0.0, 'Pulse Amplitude Absolute Value (' + r'$\rm \mu A$' + ')', ha='center')
fig2.text(0.0, 0.5, 'Scaled Recruitment', va='center', rotation='vertical')
ax2[1].legend(edgecolor='inherit', fancybox=False)
fig2.tight_layout()

# Inter-fascicular selectivity
# First, compute it
keys_nonerve = [x for x in keys if 'nerve' not in x]
nfas = len(keys_nonerve)
selectivity = {}
for k in keys_nonerve:
	selectivity[k] = {}
	for ec_key in ec_keys:
		selectivity[k][ec_key] = recruitment[k][ec_key]['fractional'] - sum([recruitment[l][ec_key]['fractional'] for l in keys_nonerve if l != k]) / (nfas - 1)

# Then show it
plt.style.use('./PaperFig_0.8_textwidth.mplstyle')

fig3, ax3 = plt.subplots()
sel = selectivity['Fascicle_01']
for ec_key in ec_keys:
	sel[ec_key] = np.array(sel[ec_key])
	ax3.plot(currents, sel[ec_key], c=linecolors[ec_key], ls='-', label='S' + ec_key.upper())
	print('Maximum selectivity for %s: %f'%(ec_key, sel[ec_key].max()))
ax3.set_xlabel('Pulse Amplitude Absolute Value (' + r'$\rm \mu A$' + ')')
ax3.set_ylabel('Selectivity')
ax3.legend(edgecolor='inherit', fancybox=False)
ax3.set_xlim(0, 4)
ax3.set_ylim(0, 1)
fig3.tight_layout()

# Recruitment differences
fig4, ax4 = plt.subplots()
for i, k in enumerate(sorted_keys):
	# Point where the recruitment is maximum
	limit = np.where(recruitment[k]['EC']['fractional'] == recruitment[k]['EC']['fractional'].max())[0].min()
	limit = np.s_[:limit]
	# Difference
	diff = recruitment[k]['EC']['fractional'] - recruitment[k]['noEC']['fractional']
	if 'nerve' in k:
		lw = 2.
	else:
		lw = 0.5
	ax4.plot(currents[limit], diff[limit], c='k', marker=markers[k], markersize=5, lw=lw, label=names[k])

ax4.set_xlabel('Pulse Amplitude Absolute Value (' + r'$\rm \mu A$' + ')')
ax4.set_ylabel('Scaled Recruitment Difference')
ax4.legend(scatterpoints=1)
fig4.tight_layout()

########################################################################
# Finish up
print('Maximum scaled recruitment differences:')
for k, v in diffs.items():
	print('\t%s, %0.3f'%(names[k], v.max()))
print('')

plt.show()
