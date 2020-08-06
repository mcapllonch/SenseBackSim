import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from datetime import datetime
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter
from matplotlib import ticker, cm

import workspace as ws
import read_results



def show_slice(
		ax, 
		filename, 
		axis, 
		where, 
		method='nearest', 
		rsx=100, 
		rsy=100, 
		rsz=100, 
		window=(0., 1.), 
		nlevels=100, 
		vmin=None, 
		vmax=None, 
		extend='neither', 
		logscale=False, 
		save=False, 
		zorder=0, 
		**kwargs
	):
	""" Opens a file 'filename' and shows the data at a slice perpendicular
	to the 'axis' axis at a relative position 'where', where 'where' ranges
	from 0. to 1. """

	# Keyword arguments
	if len(kwargs) >= 2:
		xyzunits = kwargs['xyzunits']
		vunits = kwargs['vunits']
	else:
		xyzunits = ''
		vunits = ''
	try:
		cmap = kwargs['cmap']
	except KeyError:
		cmap = plt.cm.jet
	try:
		show_data_points = kwargs['show_data_points']
	except KeyError:
		show_data_points = False
	try:
		show_intp_points = kwargs['show_intp_points']
	except KeyError:
		show_intp_points = False

	# Set colors for bad values
	if False:
		cmap.set_under(color='black')
		cmap.set_over(color='black')
	if False:
		cmap.set_bad(color='black')

	# Load data
	data = np.loadtxt(filename, delimiter=',')
	x = data[:,0]
	y = data[:,1]
	z = data[:,2]
	v = data[:,3]

	# Sides
	xspan = x.max() - x.min()
	yspan = y.max() - y.min()
	zspan = z.max() - z.min()

	# Window
	wl, wr = window
	# Crop the window limits just in case
	wl = max(wl, 0.)
	wl = min(wl, 1.)
	wr = max(wr, 0.)
	wr = min(wr, 1.)

	# Regular mesh
	xi_ = np.linspace(x.min(), x.max(), rsx)
	yi_ = np.linspace(y.min(), y.max(), rsy)
	zi_ = np.linspace(z.min(), z.max(), rsz)
	
	if axis == 'x':
		# Indices for the windows
		wl_i = int(rsx * wl)
		wr_i = int(rsx * wr)
		# Create the mask for the data
		xleft = x.min() + xspan * wl
		xrght = x.min() + xspan * wr
		mask = np.ma.masked_where((x >= xleft) & (x <= xrght), x).mask
		# Crop the corresponding array for the regular mesh
		xi_ = xi_[wl_i:wr_i]
	elif axis == 'y':
		# Indices for the windows
		wl_i = int(rsy * wl)
		wr_i = int(rsy * wr)
		# Create the mask for the data
		yleft = y.min() + yspan * wl
		yrght = y.min() + yspan * wr
		mask = np.ma.masked_where((y >= yleft) & (y <= yrght), y).mask
		# Crop the corresponding array for the regular mesh
		yi_ = yi_[wl_i:wr_i]
	elif axis == 'z':
		# Indices for the windows
		wl_i = int(rsz * wl)
		wr_i = int(rsz * wr)
		# Create the mask for the data
		zleft = z.min() + zspan * wl
		zrght = z.min() + zspan * wr
		mask = np.ma.masked_where((z >= zleft) & (z <= zrght), z).mask
		# Crop the corresponding array for the regular mesh
		zi_ = zi_[wl_i:wr_i]

	# 'mesh' the regular mesh
	xi, yi, zi = np.meshgrid(xi_, yi_, zi_)

	# Mask data
	x = x[mask]
	y = y[mask]
	z = z[mask]
	v = v[mask]

	# Interpolation
	t0 = datetime.now()
	vi = griddata((x,y,z), v, (xi,yi,zi), method=method)
	t1 = datetime.now()
	print('Time needed: %f s'%((t1 - t0).total_seconds()))

	# Figure

	# Levels for colouring
	vi_mkd_cpd = np.ma.masked_where(np.isnan(vi), vi).compressed()
	vmin_, vmax_ = vi_mkd_cpd.min(), vi_mkd_cpd.max()
	if vmax is None:
		vmax = vmax_
	if vmin is None:
		vmin = vmin_
	# levels_i = np.linspace(vmin, vmax, nlevels)
	# Avoid an error due to a flat level list
	if vmin == vmax:
		vmin -= 0.1
		vmax += 0.1
	# Re-build the levels
	print('Limits:', vmin, vmax)
	levels_i = np.linspace(vmin, vmax, nlevels)
	loglevels = np.logspace(
				np.log10(min(np.abs(vmin), np.abs(vmax))), 
				np.log10(max(np.abs(vmin), np.abs(vmax))), 
				nlevels
			)
	# Minimum and maximum exponents
	min_exp = int(np.log10(min(np.abs(vmin), np.abs(vmax))))
	max_exp = int(np.log10(max(np.abs(vmin), np.abs(vmax)))) + 1
	# Ticks for the colorbar in case it's logscale
	cbticks = 10. ** np.arange(min_exp, max_exp + 1, 1)
	cbticks_ = cbticks.tolist()
	for i in range(2, 10, 1):
		cbticks_ = cbticks_ + (float(i) * cbticks).tolist()
	cbticks_.sort()
	# print(cbticks_, min_exp, max_exp)
	cbticklabels = []
	for c in cbticks_:
		if c in (20, 50, 100, 200, 500, 1000):
			cbticklabels.append('%i'%c)
		else:
			cbticklabels.append('')

	if axis == 'x':
		s = int(where * rsx) - wl_i
		if not logscale:
			cf = ax.contourf(
					yi[:, s, :], 
					zi[:, s, :], 
					vi[:, s, :], 
					levels_i, 
					cmap=cmap, 
					extend=extend, 
					zorder=zorder
				)
		else:
			cf = ax.contourf(
					yi[:, s, :], 
					zi[:, s, :], 
					np.abs(vi[:, s, :]), 
					loglevels, 
					cmap=cmap, 
					norm=LogNorm(
						vmin=min(np.abs(vmin), np.abs(vmax)), 
						vmax=max(np.abs(vmin), np.abs(vmax))
						), 
					# locator=ticker.LogLocator(), 
					zorder=zorder
				)
		if show_data_points:
			ax.plot(y, z, 'k.', markersize=0.6)
		if show_intp_points:
			ax.scatter(yi, zi, c='k', s=1)
		ax.set_title('x = %f'%xi[s, s, s] + ' ' + xyzunits)
		ax.set_xlabel('y (' + xyzunits + ')')
		ax.set_ylabel('z (' + xyzunits + ')')
	elif axis == 'y':
		s = int(where * rsy) - wl_i
		if not logscale:
			cf = ax.contourf(
					xi[s, :, :], 
					zi[s, :, :], 
					vi[s, :, :], 
					levels_i, 
					cmap=cmap, 
					extend=extend, 
					zorder=zorder
				)
		else:
			cf = ax.contourf(
					xi[s, :, :], 
					zi[s, :, :], 
					np.abs(vi[s, :, :]), 
					loglevels, 
					cmap=cmap, 
					norm=LogNorm(
						vmin=min(np.abs(vmin), np.abs(vmax)), 
						vmax=max(np.abs(vmin), np.abs(vmax))
						), 
					# locator=ticker.LogLocator(), 
					zorder=zorder
				)
		if show_data_points:
			ax.plot(x, z, 'k.', markersize=0.6)
		if show_intp_points:
			ax.scatter(xi, zi, c='k', s=1)
		ax.set_title('y = %f'%yi[s, s, s] + ' ' + xyzunits)
		ax.set_xlabel('x (' + xyzunits + ')')
		ax.set_ylabel('z (' + xyzunits + ')')
	elif axis == 'z':
		s = int(where * rsz) - wl_i
		if not logscale:
			cf = ax.contourf(
					xi[:, :, s], 
					yi[:, :, s], 
					vi[:, :, s], 
					levels_i, 
					cmap=cmap, 
					extend=extend, 
					zorder=zorder
				)
		else:
			cf = ax.contourf(
					xi[:, :, s], 
					yi[:, :, s], 
					np.abs(vi[:, :, s]), 
					loglevels, 
					cmap=cmap, 
					norm=LogNorm(
						vmin=min(np.abs(vmin), np.abs(vmax)), 
						vmax=max(np.abs(vmin), np.abs(vmax))
						), 
					# locator=ticker.LogLocator(), 
					zorder=zorder
				)
		if show_data_points:
			ax.plot(x, y, 'k.', markersize=0.6)
		if show_intp_points:
			ax.scatter(xi, yi, c='k', s=1)
		ax.set_title('z = %f'%zi[s, s, s] + ' ' + xyzunits)
		ax.set_xlabel('x (' + xyzunits + ')')
		ax.set_ylabel('y (' + xyzunits + ')')
		ax.set_aspect('equal')

	if logscale:
		l_f = LogFormatter(10, labelOnlyBase=False)
		# cb = plt.colorbar(cf, ticks=cbticks_, format=l_f)
		cb = plt.colorbar(cf, format=l_f)
		# cb = plt.colorbar(cf, ticks=cbticks)
		# cb = plt.colorbar(cf)
		cb.set_ticks(cbticks_)
		# cb.set_ticklabels(cbticklabels)
		# cb.set_norm(LogNorm())
		# print(cbticklabels)
	else:
		cb = plt.colorbar(cf, ax=ax)

	cb.set_label(vunits)

def z_slice(ax, filename, folder, z, t, tarray):
	""" This slicing function uses a much more efficient 
	method to interpolate along the z-axiz """

	# This is from a xyz file

	# Load data
	data = np.loadtxt(filename, delimiter=',')
	x = data[:,0]
	y = data[:,1]
	z_ = data[:,2]
	v = data[:,3]

	# Set the cables
	xy_all = np.array(list(zip(x, y)))
	xy_unique = np.unique(xy_all, axis=0)
	nc = len(xy_unique)

	data_intp = np.zeros(nc)

	# Construct zprofile
	for i, xy in enumerate(xy_unique):

		# Select the corresponding rows
		rows = np.where((xy == xy_all).all(axis=1))[0]

		# Sort z_ just in case
		# argsort = np.argsort(z_[rows])
		# z_.sort()
		# Sort v accordingly
		# Not yet
		# v.sort(argsort)

		# Select z_
		z__ = z_[rows]
		z_lower = z__[np.where(z__ < z)]
		z__minus_z = np.abs(z__ - z)
		z_lower_minus_z = np.abs(z_lower - z).min()
		here = np.where(z__minus_z == z_lower_minus_z)[0][0]


		# Interpolate
		v__ = v[rows]
		data_intp[i] = (z__minus_z[here + 1] * v__[here] + z__minus_z[here] * v__[here + 1]) / (z__[here + 1] - z__[here])
		# if v[here] < -100.:

	# Show result
	ax.scatter(xy_unique[:, 0], xy_unique[:, 1], c=data_intp, s=10, cmap=plt.cm.jet_r)
	ax.set_aspect('equal')

	# KEEP THE CODE BELOW, DON'T DELETE IT.
	# It works only when there are results
	if False:
		# Read results
		data, zprofile, xyr = read_results.read_results(folder)
		xx_, yy_, rr_ = xyr

		nc = len(data)
		data_intp = np.zeros(nc)

		for i in range(nc):

			# Choose data
			try:
				data_ = data[i]['vext[0]']
			except KeyError:
				data_ = data[i]['vext[1]']

			# Find the zprofile[i] value lying just below z
			zpr = np.array(zprofile[i].values())
			zpr_minus_z = np.abs(zpr - z)
			here = np.where(zpr_minus_z == zpr_minus_z.min())

			# Interpolate v
			it = np.where(tarray == t)[0][0]
			v = data_[:, it]
			# data_intp[i] = (zpr_minus_z[here] * v[here] + zpr_minus_z[here + 1] * v[here + 1]) / (zpr[here + 1] - zpr[here])
			data_intp[i] = (zpr_minus_z[here + 1] * v[here] + zpr_minus_z[here] * v[here + 1]) / (zpr[here + 1] - zpr[here])

		# Show result
		ax.scatter(xx_, yy_, c=data_intp)

	return xy_unique, data_intp

def axons_activity(time, data, tlims):
	""" Show axons' activity """		
	labels = {'Sti': 'Stimulation point', 'Rec': 'Recording point'}
	colors = {'Sti': 'k', 'Rec': 'b'}
	fig, ax = plt.subplots()
	# Potential
	for vecs in data:
		for k, vv in vecs.items():
			ax.plot(time, vv, colors[k])
	ax.set_xlim(tlims[0], tlims[1])
	ax.set_xlabel('t (ms)')
	ax.set_ylabel('Vm (mV)')
	ax.set_title('Axons Vm (mV)')
	plt.savefig(os.path.join(ws.folders["data/results"], 'images/Axons.png'))

def cross_sec_results(data):
	""" Cross-section of the nerve and results """
	# Map of axons
	act_axs = []
	k = 0
	# Count it as having an AP if Vm is above 20 mV
	for i, v in enumerate(data):
		# Don't use the stimulation point to check if there is an AP, because
		# the stimulation artefact can give a false positive.
		# v = v['Sti']
		# Use instead the recording position. This is in better accordance to 
		# Raspopovic & al. (2011), who say that the AP needs to travel the whole
		# length for the axon to be tagged as activated
		v = v['Rec']
		if v.max() > -30.:
			act_axs.append(True)
		else:
			act_axs.append(False)

	# Figure
	fig, ax = plt.subplots()

	zorder = 0

	# Triangulation
	ws.nvt.draw_contours(ax, c='k', zorder=zorder)
	ws.nvt.draw_axons_and_points(ax, c='k', act_axons=act_axs, zorder=zorder)
	ws.nvt.draw_triangulation(ax, c='darkgrey', lw=0.4, zorder=zorder)
	zorder += 1

	# Scatter all points in black
	x, y = ws.nvt.x, ws.nvt.y
	ax.scatter(x, y, c='k', s=1, zorder=zorder)

	# Arrange figure
	xleft = x.min()
	xrght = x.max()
	ybot = y.min()
	ytop = y.max()
	margin = 100

	ax.set_xlim(xleft - margin, xrght + margin)
	ax.set_ylim(ybot - margin, ytop + margin)
	ax.set_aspect('equal')
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	title_ups = 'Tessellated nerve and current injections'
	title_len = 'Length (z-axis): %0.00f um'%ws.length
	title_inj = 'Injected current(s): %0.00f uA'%(ws.icamp * 1.e-3)
	title = title_ups + '\n' + title_len + '\n' + title_inj
	ax.set_title(title)

	# Field
	ws.nvt.pd.draw(ax, colour='dimgrey', alpha=1., linestyle='-', linewidth=.5,
			markers=None, title=title, zorder=zorder)
	fig.tight_layout()
	fig.savefig(os.path.join(ws.folders["data/results"], 'images/Nerve_and_tessellation.png'))