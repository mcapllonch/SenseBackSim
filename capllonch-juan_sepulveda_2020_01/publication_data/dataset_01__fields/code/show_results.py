"""
Open a xyz file and visualise a slice of it
Comment or uncomment lines at will for different visualizations
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm

import visualisation as vis
import workspace as ws
import anatomy


filename = 'results_vext_dataset01_time_step_00003.xyz'

number = 1

# Values of z and t
z_slice = 4000.
z_slice = 6000.
z_slice = 7000.
z_slice = 3875.
z_slice = 5000.
t = 0.015

# Value of 'where'
where = z_slice / 1.e4

fig, ax = plt.subplots()
zorder = 0

# Visualise data in colour maps
vis.show_slice(
		ax, 
		filename, 'z', 
		where=where, 
		method='linear', 
		rsx=60, 
		rsy=60, 
		rsz=150, 
		window=(where - 0.05, where + 0.05), 
		nlevels=100, 
		# vmax=-20, 
		# vmax=-5, 
		vmax=-20, 
		# vmin=-1000, 
		# vmin=-2500, 
		vmin=-2400, 
		# vmax=-1, 
		# vmin=-20, 
		# extend='both', 
		save=True, 
		xyzunits=r'$\rm \mu m$', 
		vunits='mV', 
		cmap=plt.cm.jet, 
		show_intp_points=False, 
		logscale=True, 
		zorder=zorder, 
	)
zorder += 1

# Include the contours
# Prepare the necessary stuff for the simulation
import simcontrol
# simcontrol.prepare_simulation()
simcontrol.prepare_workspace(remove_previous=False)
anatomy.create_nerve()
# Draw contours
if False:
	ws.nvt.draw_axons_and_points(
			ax, 
			facecolor='none', 
			edgecolor='k', 
			zorder=zorder
		)
if False:
	ws.nvt.draw_triangulation(
			ax, 
			c='white', 
			lw=0.5, 
			alpha=0.4, 
			zorder=zorder, 
		)
	zorder += 1
	ws.nvt.draw_contours(ax, c='k', zorder=zorder)
	zorder += 1

if False:
	vis.show_slice(
			ax, 
			filename, 'x', 
			0.95, 
			method='linear', 
			rs=100, 
			window=(0.91, 0.99), 
			nlevels=300, 
			save=True, 
			xyzunits=r'$\rm \mu m$', 
			vunits='mV', 

		)

# Visualise the longitudinal profile of the potential along the NAELC that gets stimulated, which is cable # 658 (NAELC number zero)
# Or basically, it's the cable situated at x = 250, y = 0

# fig.savefig('interp_data_cut_%s_%i.png'%(filename.split('.')[0], number))

# Load data
data = np.loadtxt(filename, delimiter=',')
x = data[:,0]
y = data[:,1]
z = data[:,2]
v = data[:,3]

plt.style.use('./PaperFig_0.8_textwidth.mplstyle')
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

for xx, yy in [[250, 0], [0, 0], [-250, 0]]:
	distances = np.hypot(x - xx, y - yy)
	axon_i = np.where(distances == distances.min())

	# Modify xx and yy so that they match the real ones
	xx = x[axon_i][0]
	yy = y[axon_i][0]
	
	ax2.plot(z[axon_i], v[axon_i], label=r'$x = %i$ '%xx + r'$\rm \mu m$')
	ax3.plot(z[axon_i], np.abs(v[axon_i]), label=r'$x = %i$ '%xx + r'$\rm \mu m$')

ax2.set_xlim(0.2875e4, 0.7125e4)
ax2.set_xlabel(r'$\rm z$' + ' (' + r'$\rm \mu m$' + ')')
ax2.set_ylabel(r'$\rm v_{E}$' + ' (' + r'$\rm mV$' + ')')
# ax2.legend(loc='best', fontsize=8)
ax2.legend(loc='best')
# ax2.set_xlabel('z')
# ax2.set_ylabel('v')

ax3.set_xlim(0.2875e4, 0.7125e4)
ax3.set_ylim(1.e0, 2.e3)
ax3.set_xlabel(r'$\rm z$' + ' (' + r'$\rm \mu m$' + ')')
ax3.set_yscale('log')
ax3.set_ylabel(r'$|v_{E}|$' + ' (' + r'$\rm mV$' + ')')
# ax3.legend(loc='best', fontsize=8)
ax3.legend(loc='best')

# fig2.savefig('vext_zprofile.png', bbox_inches='tight')
# fig3.savefig('vext_zprofile_logscale.png', bbox_inches='tight')

# Results folder
folder = os.path.join(os.getcwd(), "data/results")

# Time
nt = 251
dt = 0.005
tarray = np.arange(0, nt * dt, dt)

# Create figure
plt.style.use('../PaperFig_1textwidth_v4.mplstyle')
fig4, ax4 = plt.subplots()
xy, data_intp = vis.z_slice(ax4, filename, folder, z=z_slice, t=0.015, tarray=tarray)

# Identify the axons with their positions
pd = ws.nvt.pd
xy_indices = []

for xx, yy in zip(pd.x, pd.y):
	distances = np.hypot(xy[:, 0] - xx, xy[:, 1] - yy)
	i = np.where(distances == distances.min())[0][0]
	xy_indices.append(i)

# Sort data_intp
data_intp = data_intp[xy_indices]

# Draw the power diagram
ws.nvt.pd.draw(
		ax4, 
		draw_polygons=False, 
		linewidth=0.5, 
		line_alpha=1., 
		colour='same_as_facecolors', 
		# colour='k', 
		values='precomputed', 
		precompv=data_intp, 
		# cmap=plt.cm.jet, 
		# cmap=plt.cm.gist_stern, 
		# cmap=plt.cm.afmhot, 
		# cmap=plt.cm.afmhot_r, 
		cmap=plt.cm.gnuplot, 
		vmin=25., 
		# vmin=5., 
		# vmin=20., 
		# # vmax=1424.045, 
		vmax=1000., 
		# vmax=2400., 
		# vmin=1, 
		# vmax=20, 
		extend='max', 
		title=r'$|v_{E}|$ at ' + r'$z =$ $%i$ $\rm \mu m$ and $t =$ $%0.3f$ $\rm ms$'%(z_slice, t), 
		cblabel=r'$|v_{E}|$ ($\rm mV$)', 
		logscale=True, 
		allowed_cbticklabels=(10, 30, 100, 300, 1000), 
		zorder=zorder, 
	)
zorder += 1
ws.nvt.draw_contours(ax4, c='k', zorder=zorder)
zorder += 1

# Show diamond at the position of the active pad
ax4.scatter(250, 0, marker='D', c='lime', s=30, zorder=zorder)
# Aspect ratio
ax4.set_aspect('equal')

# Print the limits of v
vabs = np.abs(v)
v_masked = np.ma.masked_where(vabs < 1.e-6, vabs)

print('v_min:', v_masked.min())
print('v_mean:', v_masked.mean())
print('v_max:', v_masked.max())

vabs = np.abs(v_masked)
print('|v|_min:', vabs.min())
print('|v|_mean:', vabs.mean())
print('|v|_max:', vabs.max())


plt.show()
plt.close('all')