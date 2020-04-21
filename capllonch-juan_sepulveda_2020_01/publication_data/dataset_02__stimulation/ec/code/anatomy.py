"""
Handle all the anatomical aspects of the nerve
z-profile:

Cross-section:
Fill the fascicles and discretise the whole cross-section with a 
tessellation
"""
import os
import json
import random
import numpy as np
import matplotlib.pyplot as plt
import triangle
from collections import OrderedDict
import planar as pl
import copy
import warnings
import csv

import workspace as ws
import contourhandler as cth
import tessellations as tess
import circlepacker as cp
import geometry as geo
import tools

# Ignore warnings
warnings.filterwarnings('ignore')



###############################################################################
# CROSS-SECTION

class NerveTess():
	""" 
	Tessellation of a nerve's cross-section, including:
	 - Power diagram
	 - Its dual triangulation
	"""
	def __init__(self):
		self.cpath = ws.anatomy_settings["cross-section"]["contours file"]
		
		# Dictionaries for:
		# The fascicles to which each cable belongs (None if necessary)
		self.cables_fascicles = OrderedDict()
		# Tissues surrounding each cable
		self.cables_tissues = OrderedDict()

	def build_contours(self):
		"""
		Build the contours of the nerve
		Contour processing
		We obtain five contour variables:
		1. contour: The actual contours dictionary we are going to use. Its number of
						points may be lower than it contains in the file.
		2. contour_hd: The contours with all their points, no reduction.
		3. contour_pslg: The contours in PSLG format
		4. contour_nerve: Contour for the nerve only, in order to triangulate it and
						fill it with points.
		5. contour_pslg_nerve: Contour for the nerve only in PSLG format.
		"""

		c_reduction = self.c_reduction
		# Open contour
		contour = cth.load_contours(self.cpath)
		# Save the original contours, I will need them to prevent axons from protruding
		# out of the fascicles
		contour_hd = copy.deepcopy(contour)
		# Take just one fraction of the points in case they are too many
		contour = cth.reduce_points(contour, c_reduction)
		# Convert to PSLG
		contour_pslg = cth.c2pslg(contour)
		# Ony the nerve
		contour_nerve = {'Nerve': contour['Nerve']}
		contour_pslg_nerve = cth.c2pslg(contour_nerve)

		# Save all the types of contours
		self.contour = contour
		self.contour_nerve = contour_nerve
		self.contour_hd = contour_hd
		self.contour_pslg = contour_pslg
		self.contour_pslg_nerve = contour_pslg_nerve

		# Save other properties
		self.polygon = geo.Polygon(np.array(contour_nerve['Nerve']))
		# self.centroid = self.polygon.centroid
		self.crss_area = self.polygon.area
		# self.polygon.get_circumcircle()
		self.circumcircle = self.polygon.circumcircle
		self.circumcenter = self.circumcircle.c
		self.circumradius = self.circumcircle.r
		self.circumdiameter = self.polygon.circumdiameter

		# Instantiate fascicles
		self.build_fascicles()
		
	def build(self, params):
		"""
		Build the tessellation for the nerve
		"""

		# Get parameters
		# self.get_params(params)
		self.__dict__.update(params)

		mnd = self.mnd
		min_sep = self.min_sep
		rmin = self.rmin
		rmax = self.rmax
		circp_tol = self.circp_tol
		max_axs_pf = self.max_axs_pf
		numberofaxons = self.numberofaxons
		locations = self.locations
		radii = self.radii
		# models = self.models

		# Build contours
		self.build_contours()

		contour = self.contour
		contour_hd = self.contour_hd
		contour_pslg = self.contour_pslg
		contour_pslg_nerve = self.contour_pslg_nerve

		##################################################################
		# Triangulate

		# Max. area for the triangles, 
		# inverse to the minimum NAELC density
		maxarea = 1. / mnd
		# Triangulate
		# Instead, triangulate without taking the fascicles' contours
		# into account
		tri = triangle.triangulate(contour_pslg_nerve, 'a%f'%maxarea)
		tri = triangle.triangulate(contour_pslg, 'a%f'%maxarea)
		# Vertices
		tv = tri['vertices']

		self.original_points = tv.T

		##################################################################
		# Fill fascicles

		# If the axons have fixed locations, get their locations and
		# radii first, and then they will be added to the different 
		# fascicles when needed
		if self.packing_type == "fixed locations":
			try:
				xx_, yy_ = np.array(locations).T
			except ValueError:
				# Something went wrong or simply there are no axons
				xx_ = yy_ = rr_ = np.array([])
			else:
				rr_ = np.array(radii)

			# Axon models
			self.models = self.models['fixed']

		# else:
		# 	# Axon models. Create them according to the proportions


		# Remove points inside the fascicles
		remove_these = []
		axons = {}
		naxons = {}
		naxons_total = 0
		print('about to fill the contours')
		for k, v in contour.items():
			varr = np.array(v)
			# Remove points outside the nerve
			if 'Nerve' in k:
				plpol = pl.Polygon(v)
				for i, p in enumerate(tv):
					# Remove any points in tv falling outside the nerve
					# and outside its boundaries
					# Note: that means that I don't remove the points ON the
					# boundaries
					if not (plpol.contains_point(p) or np.isin(p, varr).all()):
						remove_these.append(i)
			# Remove points from the fascicles
			if 'Fascicle' in k:
				print(k)
				plpol = pl.Polygon(v)
				for i, p in enumerate(tv):
					# Remove the points of the fascicle's contours from tv
					inclc = plpol.contains_point(p) and (not np.isin(p, varr).all())
					# Actually, don't remove the fascicle's contours
					# Remove any points in tv contained in a fascicle
					notic = plpol.contains_point(p) or np.isin(p, varr).all()
					if notic:
						remove_these.append(i)
				# Fill the fascicle
				# Dictionary of axons
				axons[k] = {}
				# Create the circles
				# Different packing strategies yield different results
				if self.packing_type == "uniform":
					xx, yy, rr = cp.fill(contour_hd[k], rmin, rmax, 
						min_sep, nmaxpf=max_axs_pf, tolerance=circp_tol)
				elif self.packing_type == "gamma":
					distr_params = {
						"mean": self.avg_r, 
						"shape": self.gamma_shape
					}
					xx, yy, rr = cp.fill(contour_hd[k], rmin, rmax, 
						min_sep, nmaxpf=max_axs_pf, tolerance=circp_tol, 
						distribution="gamma", distr_params=distr_params)
					print('Filled %s'%k)
				elif self.packing_type == "fixed locations":
					# Iterate over axons and get those inside 
					# the fascicle
					xx, yy, rr = [], [], []
					for x_, y_, r_ in zip(xx_, yy_, rr_):
						if plpol.contains_point((x_, y_)):
							xx.append(x_)
							yy.append(y_)
							rr.append(r_)
					xx = np.array(xx)
					yy = np.array(yy)
					rr = np.array(rr)

				# Store information in a clean way
				axons[k]['x'] = xx
				axons[k]['y'] = yy
				axons[k]['r'] = rr
				# axons[k]['models'] = models[:]
				naxons[k] = len(xx)
				naxons_total += naxons[k]
				
		rps = tv[remove_these]

		self.removed_points = rps

		keep_these = np.array(list(set(range(tv.shape[0])) - set(remove_these)))
		# print("keep_these:", keep_these)
		# print("remove_these:", remove_these)
		tv = tv[keep_these]

		# List x, y and r once some points have been removed
		tv = tv.T
		x, y = tv
		nc = x.size
		r = np.zeros_like(x)

		# Dictionaries for:
		# cables (their type)

		cables = OrderedDict()
		models = {}
		
		# Points corresponding to epineurium
		for ic in range(nc):
			cables[ic] = 'NAELC'
			# NAELC model indexing
			models[ic] = 'NAELC'
			print(ic, models, self.models)
		
		# Add axons to the existing points for the nerve
		for k in axons:
			x = np.array(x.tolist() + axons[k]['x'].tolist())
			y = np.array(y.tolist() + axons[k]['y'].tolist())
			r = np.array(r.tolist() + axons[k]['r'].tolist())

		nc = x.size
		nNAELC = nc - naxons_total

		# Axon models
		if self.packing_type != 'fixed locations':
			# Now that the axons have been placed, determine their models according to the proportions
			# ninst: 'number of instances' (of each model)
			ninst = {}
			proportions = self.models['proportions']
			keys = list(proportions.keys())
			for m, p in proportions.items():
				ninst[m] = int(p * naxons_total)
			# Remainder
			rem = naxons_total - sum(list(ninst.values()))
			
			# Just add the remaining in an arbitrary (not random) way
			for i in range(rem):
				ninst[keys[i]] += 1

			# Now select the indices of the axons for each model
			axon_indices = (np.arange(naxons_total) + nNAELC).tolist()
			remaining_axon_indices = axon_indices[:]
			# Dictionary for the axon indices for each model
			inds = {}
			# Dictionary for the model names for each index
			model_by_index = {}
			for m, p in proportions.items():
				sample = random.sample(remaining_axon_indices, ninst[m])
				remaining_axon_indices = list(set(remaining_axon_indices) - set(sample))
				inds[m] = sample[:]
				for i in sample:
					model_by_index[i] = m
				
		# Points corresponding to axons
		for ik, i in enumerate(np.arange(nNAELC, nc, 1)):
			cables[i] = 'Axon'
			# Axon model indexing
			if self.packing_type == 'fixed locations':
				models[i] = self.models[ik]
			else:
				models[i] = model_by_index[i]
			print(ik, i, models, self.models)

		##################################################################
		# Power diagram
		# Zero-valued radii: Voronoi diagram

		# Build power diagram
		pd = tess.PowerDiagram(x, y, r, contour_pslg_nerve)
		pd.build()

		# Lengths of the segments and the connections
		pdpairs = pd.pairs
		pdsegments = pd.segments
		pairs = []
		segments = {}
		len_seg = {}
		len_con = {}
		for pair in pdpairs:
			# if len(pdsegments[pair]) > 1:
			if pdsegments[pair] is not None:
				seg = pdsegments[pair]
				pairs.append(pair)
				segments[pair] = seg
				a, b = seg.a, seg.b
				len_seg[pair] = geo.dist(a, b).mean()
				i, j = pair
				len_con[pair] = geo.dist((x[i], y[i]), (x[j], y[j]))

		# Store relevant stuff in the attributes
		self.pd = pd
		self.x = x
		self.y = y
		self.r = r
		self.nc = nc
		self.axons_dict = axons
		self.pairs = pairs
		self.cables = cables
		self.models = models
		# Unique list of models
		self.models_set = set(self.models.values())
		self.segments = segments
		self.trios = pd.trios
		self.len_con = len_con
		self.len_seg = len_seg
		self.free_areas = pd.free_areas
		self.circ_areas = pd.circ_areas
		# Total endoneurial cross-sectional free area
		self.endo_free_cs_area = self.fas_total_area - self.circ_areas.sum()
		self.naxons = naxons
		self.naxons_total = naxons_total
		self.nNAELC = nNAELC

	def save_to_file(self, spath):
		"""
		Save the nerve with its properties to a file
		"""
		with open(spath, 'w') as f:
			fw = csv.writer(f, delimiter=';')
			# Headers
			fw.writerow(['Cable number, x, y, r, free extracellular area, start position, endoneurium, epineurium'])
			# Cables
			for i in range(self.nc):
				fw.writerow([self.cables[i], self.x[i], self.y[i], self.r[i], self.free_areas[i], self.start_positions[i], self.cables_tissues[i]['epineurium'], str(self.cables_tissues[i]['endoneurium'])])
			# Pairs
			for pair in self.pairs:
				i, j = pair
				fw.writerow(['Pair', i, j] + \
					[tools.arrtonum(item) for sublist in self.segments[pair].points_list for item in sublist] + \
					[self.len_seg[pair], self.len_con[pair]])

	def save_to_json(self, path):
		"""
		New function. 7 October 2019.
		Save the nerve with its properties to a file
		NOT FINISHED
		"""
		print('saving json in: %s'%path)
		# Create dictionary to be dumped
		topology = OrderedDict()
		topology['cables'] = OrderedDict()
		topology['pairs'] = OrderedDict()

		for i in range(self.nc):
			if isinstance(i, np.int64):
				print('yes', i)
				# This is necessary since json can't serialise np.int64
				i = int(i)
			topology['cables'][i] = OrderedDict()
			topology['cables'][i]['type'] = self.cables[i]
			topology['cables'][i]['x'] = self.x[i]
			topology['cables'][i]['y'] = self.y[i]
			topology['cables'][i]['r'] = self.r[i]
			topology['cables'][i]['model'] = self.models[i]
			topology['cables'][i]['free extracellular area'] = self.free_areas[i]
			topology['cables'][i]['start position'] = self.start_positions[i]
			for s in ('endoneurium', 'epineurium'):
				topology['cables'][i][s] = self.cables_tissues[i][s]

		for pair in self.pairs:
			i, j = pair
			if isinstance(i, np.int64):
				i = int(i)
			if isinstance(i, np.integer):
				i = int(i)
			if isinstance(j, np.int64):
				j = int(j)
			if isinstance(j, np.integer):
				j = int(j)
			seg = self.segments[pair]
			str_pair = str(pair)
			topology['pairs'][str_pair] = {
					'pair': (i, j), 
					'separator segment': OrderedDict()
				}
			topology['pairs'][str_pair]['separator segment']['a'] = OrderedDict()
			topology['pairs'][str_pair]['separator segment']['a']['x'] = seg.a[0]
			topology['pairs'][str_pair]['separator segment']['a']['y'] = seg.a[1]
			topology['pairs'][str_pair]['separator segment']['b'] = OrderedDict()
			topology['pairs'][str_pair]['separator segment']['b']['x'] = seg.b[0]
			topology['pairs'][str_pair]['separator segment']['b']['y'] = seg.b[1]
			
		# Dump
		with open(path, 'w') as f:
			json.dump(topology, f, indent=4)

	def build_preexisting(self):
		"""
		Build the nerve using the parameters stored in files
		"""

		# Build contours
		self.c_reduction = ws.anatomy_settings["cross-section"]["contours point reduction"]
		self.build_contours()

		contour = self.contour
		contour_hd = self.contour_hd
		contour_pslg = self.contour_pslg
		contour_pslg_nerve = self.contour_pslg_nerve

		# Build internal elements
		x = []
		y = []
		r = []
		cables = []
		free_areas = []
		start_positions = []
		cables_tissues = OrderedDict()
		segments = {}
		len_seg = {}
		len_con = {}
		numberof = {
			'Axon': 0,
			'NAELC': 0
		}
		itpath = ws.anatomy_settings["cross-section"]["internal topology file"]
		with open(itpath, 'r') as f:
			# Skip header
			frl = list(csv.reader(f, delimiter=';'))[1:]
			# k is a cable counter
			k = 0
			for row in frl:
				key = row[0]
				if key != 'Pair':
					cables.append(key)
					x.append(float(row[1]))
					y.append(float(row[2]))
					r.append(float(row[3]))
					free_areas.append(float(row[4]))
					start_positions.append(float(row[5]))
					try:
						cables_tissues[k] = {
							'epineurium': float(row[7]), 
							'endoneurium': float(row[6])
						}
					except IndexError:
						# There's no such information
						cables_tissues[k] = {
							'epineurium': 0., 
							'endoneurium': 0.
						}
					numberof[key] += 1
					k += 1
				else:
					i = int(row[1])
					j = int(row[2])
					segments[(i, j)] = geo.Segment([(float(row[3]), float(row[4])), (float(row[5]), float(row[6]))])
					len_seg[(i, j)] = float(row[7])
					len_con[(i, j)] = float(row[8])

		# Save things as attributes
		
		self.x = np.array(x)
		self.y = np.array(y)
		self.r = np.array(r)
		self.free_areas = np.array(free_areas)
		self.start_positions = np.array(start_positions)
		self.cables_tissues = cables_tissues
		self.segments = segments
		self.len_seg = len_seg
		self.len_con = len_con
		self.pairs = len_con.keys()
		self.cables = cables
		self.nc = len(cables)
		self.naxons_total = numberof['Axon']
		self.nNAELC = numberof['NAELC']
		# Build power diagram
		pd = tess.PowerDiagram(self.x, self.y, self.r, contour_pslg_nerve)
		for p, s in zip(self.pairs, self.segments.values()):
			print(p, str(s))
		pd.build_preexisting(self.pairs, self.segments)
		self.trios = pd.trios
		self.pd = pd
		self.circ_areas = pd.circ_areas
		# Total endoneurial cross-sectional free area
		self.endo_free_cs_area = self.fas_total_area - self.circ_areas.sum()

	def build_from_json(self):
		"""
		Build the nerve using the parameters stored in json files
		"""

		# Build contours
		self.c_reduction = ws.anatomy_settings["cross-section"]["contours point reduction"]
		self.build_contours()

		contour = self.contour
		contour_hd = self.contour_hd
		contour_pslg = self.contour_pslg
		contour_pslg_nerve = self.contour_pslg_nerve

		# Build internal elements
		x = []
		y = []
		r = []
		cables = OrderedDict()
		free_areas = []
		start_positions = []
		cables_tissues = OrderedDict()
		segments = {}
		len_seg = {}
		len_con = {}
		numberof = {
			'Axon': 0,
			'NAELC': 0
		}
		models = {}

		itpath = ws.anatomy_settings["cross-section"]["internal topology file"]
		topology = read_from_json(itpath, object_pairs_hook=OrderedDict)

		# Read the dictionary and crete the necessary variables from it

		for i, c in topology['cables'].items():
			i = int(i)
			cables[i] = c['type']
			x.append(c['x'])
			y.append(c['y'])
			r.append(c['r'])
			free_areas.append(c['free extracellular area'])
			start_positions.append(c['start position'])
			cables_tissues[i] = OrderedDict()
			cables_tissues[i]['endoneurium'] = c['endoneurium']
			cables_tissues[i]['epineurium'] = c['epineurium']
			numberof[c['type']] += 1
			# Axon model
			try:
				models[i] = c['model']
			except KeyError:
				# It's not an axon
				models[i] = cables[i]

		# Turn the cables dictionary into a sorted list
		# And also sort everything else
		sortorder = np.argsort(np.array(list(cables.keys())))
		cables = np.array(list(cables.values()))[sortorder]
		x = np.array(x)[sortorder]
		y = np.array(y)[sortorder]
		r = np.array(r)[sortorder]
		free_areas = np.array(free_areas)[sortorder]
		start_positions = np.array(start_positions)[sortorder]

		# Now pairs...
		for p in topology['pairs'].values():
			i, j = p['pair']
			pair = (i, j)
			s = p['separator segment']
			a = s['a']
			b = s['b']
			seg = geo.Segment((
					(a['x'], a['y']), 
					(b['x'], b['y'])
				))
			segments[pair] = seg
			len_seg[pair] = seg.length
			len_con[pair] = geo.dist(
					(x[i], y[i]), 
					(x[j], y[j])
				)
		
		# Save things as attributes
		# self.x = np.array(x)
		# self.y = np.array(y)
		# self.r = np.array(r)
		# self.free_areas = np.array(free_areas)
		# self.start_positions = np.array(start_positions)
		self.x = x
		self.y = y
		self.r = r
		self.free_areas = free_areas
		self.start_positions = start_positions
		self.cables_tissues = cables_tissues
		self.segments = segments
		self.len_seg = len_seg
		self.len_con = len_con
		self.pairs = len_con.keys()
		self.cables = cables
		self.models = models
		# Unique list of models
		self.models_set = set(self.models.values())
		# print('cables:', cables)
		# print('r:', r)
		# for c_, r_ in zip(cables, r):
		# 	print(c_, r_)
		self.nc = len(cables)
		self.naxons_total = numberof['Axon']
		self.nNAELC = numberof['NAELC']
		# Build power diagram
		pd = tess.PowerDiagram(self.x, self.y, self.r, contour_pslg_nerve)
		# for p, s in zip(self.pairs, self.segments.values()):
		# 	print(p, str(s))
		pd.build_preexisting(self.pairs, self.segments)
		self.trios = pd.trios
		self.pd = pd
		self.circ_areas = pd.circ_areas
		# Total endoneurial cross-sectional free area
		self.endo_free_cs_area = self.fas_total_area - self.circ_areas.sum()

	def build_fascicles(self):
		""" Create instances of the class Fascicle for this nerve 
		This method was created on 13 November 2018 """
		fascicles = OrderedDict()
		fas_total_area = 0.
		for key, value in self.contour.items():
			if 'Fascicle' in key:
				fas = Fascicle(key, value)
				fascicles[key] = fas
				fas_total_area += fas.area
		# Store attributes
		# List of Fascicle instances
		self.fascicles = fascicles
		# Sum of fascicular areas
		self.fas_total_area = fas_total_area
		# All the area outside the fascicles
		self.interfas_area = self.crss_area - fas_total_area

	def allocate_cables_in_fascicles(self):
		""" Allocate cables inside their corresponding fascicles.
		Also find the weighted area lying in each tissue for each cable """
		
		# Iterate over cables
		for i in range(self.nc):
			# By default, no fascicle (it will remain so if that's the case)
			self.cables_fascicles[i] = None


			# if True:
			if ws.anatomy_settings["cross-section"]["use internal topology"] == "generation":

				# Tissues surrounding the cable (as a fraction of the free area)
				self.cables_tissues[i] = OrderedDict()
				self.cables_tissues[i]['epineurium'] = 0.
				self.cables_tissues[i]['endoneurium'] = OrderedDict()

			# Find the fascicle to which this cable belongs
			for k, fas in self.fascicles.items():
				if fas.polygon.plpol.contains_point((self.x[i], self.y[i])):
					self.cables_fascicles[i] = k
					self.fascicles[k].add_cable(i)
					# break

				# Compute the (weighted) area overlapping with this fascicle if we haven't this information yet; this should only be the case when the nerve is first generated

				# Find the intersection
				# if True:
				if ws.anatomy_settings["cross-section"]["use internal topology"] == "generation":

					# Initialise it
					self.cables_tissues[i]['endoneurium'][k] = 0.

					# The endoneurium area is found from the overlap with all fascicles

					# Polygons intersection
					intersection = geo.intersection_polygons(
							self.pd.polygons[i], 
							fas.polygon
						)
					# Add this area (if there is any intersection)
					# Note that the cable's circular area needs to be subtracted
					# If it's a NAELC, this is zero. If it's an axon, it needs to 
					# be subtracted in full since, for sure, the axon is fully 
					# inside the fascicle
					if intersection is not None:
						self.cables_tissues[i]['endoneurium'][k] = (intersection.area - self.pd.circ_areas[i]) / self.pd.free_areas[i]

			# Finally, the epineurium (weighted) area
			# if True:
			if ws.anatomy_settings["cross-section"]["use internal topology"] == "generation":
				
				self.cables_tissues[i]['epineurium'] = 1. - sum(self.cables_tissues[i]['endoneurium'].values())

			# Print results
			# if True:
			if ws.anatomy_settings["cross-section"]["use internal topology"] == "generation":
				print('tissues for cable %i (%s): epi: %f, endo: %s'%(i, self.cables[i], self.cables_tissues[i]['epineurium'], str(self.cables_tissues[i]['endoneurium'])))

	def draw_contours(self, ax, c='k', lw=1., zorder=0):
		"""
		Draw the contours of the nerve
		"""
		for v in self.contour.values():
			varr = np.array(v).T.tolist()
			xx, yy = varr
			xx.append(varr[0][0])
			yy.append(varr[1][0])
			ax.plot(xx, yy, c, lw=lw, zorder=zorder)

	def draw_axons_and_points(self, ax, act_axons=None, facecolor='k', 
		edgecolor='k', lw=1., s=1, text=False, zorder=0):
		"""
		Draw the axons
		"""
		nc, x, y, r = self.nc, self.x, self.y, self.r
		jaxon = 0
		for i in range(nc):
			if r[i] == 0:
				ax.scatter(x[i], y[i], c='k', s=s, zorder=zorder)
			# Reset color to default
			col = facecolor
			ec = edgecolor
			# If this is an axon...
			if r[i] > 0:
				# If I measured activation...
				if act_axons is not None:
					# If this axon is activated,
					if act_axons[jaxon]:
						col = 'lightgreen'
						col = 'lime'
						col = 'red'
						col = 'magenta'
						ec = 'red'
				jaxon += 1
			# Circle
			circle = plt.Circle((x[i], y[i]), r[i], facecolor=col, alpha=1., 
				zorder=zorder + 1, linewidth=lw, edgecolor=ec)
			ax.add_artist(circle)
			# Text
			if text:
				ax.text(x[i], y[i], i)

	def draw_triangulation_deprecated(self, ax, lw=1., c='k', zorder=0):
		"""
		Draw the triangulation of all the points
		"""
		x, y, trios = self.x, self.y, self.trios
		for trio in trios:
			print("Trio:", trio)
		ax.triplot(x, y, trios, color=c, linewidth=lw, zorder=zorder)

	def draw_triangulation(self, ax, lw=1., ls='-', c='k', alpha=1, 
		dashes=(1, 1), zorder=0):
		"""
		Draw the triangulation of all the points
		"""
		for pair in self.pairs:
			i, j = pair
			if ls == '--':
				ax.plot((self.x[i], self.x[j]), (self.y[i], self.y[j]), 
					color=c, lw=lw, ls=ls, dashes=dashes, alpha=alpha, zorder=zorder)
			else:
				ax.plot((self.x[i], self.x[j]), (self.y[i], self.y[j]), 
					color=c, lw=lw, ls=ls, alpha=alpha, zorder=zorder)

	def print_properties(self):
		""" Print properties of the nerve, same as with the fascicles"""
		ws.log('--------------------------------------')
		ws.log('Nerve')
		ws.log('Circumcircle\'s Diameter: %0.2f'%self.circumdiameter)
		ws.log('Circumcircle\'s Center: (%f, %f)'%(self.circumcenter[0], self.circumcenter[1]))
		ws.log('Area: %0.3f um2'%self.crss_area)
		ws.log('Number of Axons: %i'%self.naxons_total)
		ws.log('Number of NAELC: %i'%self.nNAELC)
		ws.log('Total number of cables: %i'%self.nc)
		ws.log('--------------------------------------')


class Fascicle():
	"""
	This is the class for a fascicle. This class has been created 
	in order to easily access some properties of it.
	This was created on 13 November 2018
	"""
	def __init__(self, key, contour):
		self.id = key
		self.contour = contour
		self.polygon = geo.Polygon(contour)
		self.polygon.get_circumcircle()
		self.circumcircle = self.polygon.circumcircle
		self.circumcenter = self.circumcircle.c
		self.circumdiameter = self.polygon.circumdiameter
		self.area = self.polygon.area
		self.cables = []
		self.ncables = 0
		self.naxons = 0

	def add_cable(self, i):
		""" Add one cable to this fascicle """
		self.cables.append(i)
		self.ncables += 1
		if ws.nvt.cables[i] == "Axon":
			self.naxons += 1

	def get_free_area(self):
		""" Compute the fiber area (area occupated by the fibers) and the free (free from axons) endoneurial area 
		of this fascicle. This can only be done when all the 
		cables have been allocated """

		# First, fiber area
		self.fiber_area = 0.

		# for i in self.cables:
		# 	self.fiber_area += ws.nvt.r[i] ** 2
		# self.fiber_area *= np.pi
		self.fiber_area = np.pi * (ws.nvt.r[self.cables] ** 2).sum()
		# print(self.fiber_area)

		# Second, free extracellular area
		self.free_area = self.area - self.fiber_area

		# Packing ratio
		self.packing_ratio = self.fiber_area / self.area

	def get_intracellular_area(self):
		""" Add the intracellular area. Skip this if it already exists """
		try:
			self.intracellular_area
		except AttributeError:
			# Create it if it does not exist
			self.intracellular_area = 0.
			for i in self.cables:
				if ws.nvt.cables[i] == "Axon":
					self.intracellular_area += ws.cables[i].intracellular_area
					
		# Intracellular to extracellular areas ratio
		self.intra_extra_ratio = self.intracellular_area / self.free_area

	def print_properties(self):
		""" Print the main properties of the fascicle: 
		size, packing ratios, number of axons, etc."""
		try:
			ws.log('--------------------------------------')
			ws.log('Fascicle %s'%self.id)
			ws.log('Circumcircle\'s Diameter: %0.2f'%self.circumdiameter)
			ws.log('Circumcircle\'s Center: (%f, %f)'%(self.circumcenter[0], self.circumcenter[1]))
			# ws.log('Centroid: (%f, %f)'%(self.polygon.centroid[0], self.polygon.centroid[1]))
			ws.log('Area: %0.3f um2'%self.area)
			ws.log('Number of Axons: %i'%self.naxons)
			ws.log('Fiber Area: %0.3f um2'%self.fiber_area)
			ws.log('Extracellular Area: %0.3f um2'%self.free_area)
			ws.log('Fiber Packing Ratio: %0.3f'%self.packing_ratio)
			ws.log('Intracellular Area: %0.3f um2'%self.intracellular_area)
			ws.log('Intracellular to Extracellular Areas Ratio: %0.3f'%self.intra_extra_ratio)
			ws.log('--------------------------------------')
		except AttributeError:
			ws.log('ERROR: Fascicle %s is missing certain attributes that were intended to print. Skipping the rest.'%self.id)

###############################################################################
# Z-AXIS

class Section():
	"""Section"""
	def __init__(self, sectype, length, xstart):
		self.sectype = sectype
		self.length = length
		self.xstart = xstart

class Axon():
	""" Axon, which contains a chain of sections """
	def __init__(self, sections):
		self.sections = sections
		
########################################################################
# z-profile
########################################################################
# Just to select the section types with numbers
# This is more convenient for hoc, I think

def read_from_json(path, object_pairs_hook=None):
	""" Read topology file in a json format """
	with open(path, 'r') as f:
		return json.load(f, object_pairs_hook=object_pairs_hook)

def create_nerve():
	"""
	Create the nerve using the available information
	Also create the necessary variables for the model to work and to be 
	able to perform further actions
	"""

	ws.log("Creating Nerve")

	# Build a nerve from the specifications in the anatomy settings

	# Contours

	# Generation
	if ws.anatomy_settings["cross-section"]["use contours"] == "generation":

		generate_contours()
		# Change the settings so we load the just generated file
		ws.anatomy_settings["cross-section"]["contours file"] = os.path.join(
			ws.folders["data/saved"], 
			ws.anatomy_settings["cross-section"]["contours generation"]["filename"])

	# Load from a file
	elif ws.anatomy_settings["cross-section"]["use contours"] == "load file":
		pass

	# Create the tessellated nerve
	ws.nvt = NerveTess()
	nvt = ws.nvt
	
	# Filling
	
	# Generation
	if ws.anatomy_settings["cross-section"]["use internal topology"] == "generation":
		# Create the cables and tessellation
		print('generation of internal topology')
		fill_contours()
		ws.anatomy_settings["cross-section"]["internal topology file"] = \
			os.path.join(ws.folders["data/saved"], 
				ws.anatomy_settings["cross-section"]["fibers distribution"]["filename"])

	# Load from file and build
	if ws.anatomy_settings["cross-section"]["use internal topology"] == "load file":
		itpath = ws.anatomy_settings["cross-section"]["internal topology file"]
		if '.csv' in itpath:
			ws.nvt.build_preexisting()
		elif '.json' in itpath:
			ws.nvt.build_from_json()

	# Allocate cables to fascicles
	ws.nvt.allocate_cables_in_fascicles()

	# Save the generated topology
	if ws.anatomy_settings["cross-section"]["use internal topology"] == "generation":

		savepath = ws.anatomy_settings["cross-section"]["fibers distribution"]["filename"]

		# To a csv file
		spath = os.path.join(ws.folders["data/saved"], savepath.replace('.json', '.csv'))
		ws.nvt.save_to_file(spath)

		# To a json file
		if ws.anatomy_settings["cross-section"]["use internal topology"] == "generation":
			jsonpath = os.path.join(ws.folders["data/saved"], savepath.replace('.csv', '.json'))
			ws.nvt.save_to_json(jsonpath)

	# Get fascicles' free areas
	for k in ws.nvt.fascicles:
		ws.nvt.fascicles[k].get_free_area()

	# Get the necessary stuff from the nerve
	ws.nNAELC = nvt.nNAELC
	ws.naxons_total = nvt.naxons_total
	ws.nc = nvt.nc
	r = nvt.r
	ws.pairs = nvt.pairs
	contour_nerve = nvt.contour_nerve

	if ws.settings["graphics"]:
		fig, ax = plt.subplots()
		diams = 2. * r[np.where(r > 0)]
		count, bins, ignored = ax.hist(diams, 50, normed=True)
		binswidth = bins[1:] - bins[:-1]
		import scipy.special as sps
		if ws.anatomy_settings["cross-section"]["use internal topology"] == "generation" \
		and ws.anatomy_settings["cross-section"]["fibers distribution"]["axon placements"]["packing"]["type"] == "gamma":
			ws.log("Axon diameters:")
			ws.log("\tThis should be close to 1: %f:"%(count * binswidth).sum())
			ws.log("\tThis should be close to %f: %f"%(2 * ws.nvt.avg_r, diams.mean()))
			shape = ws.nvt.gamma_shape
			mean = 2. * ws.nvt.avg_r
			scale = mean / shape
		else:
			shape = 2.5
			mean = 2. * 3.65
			scale = mean / shape
		y = bins**(shape-1)*(np.exp(-bins/scale) / (sps.gamma(shape)*scale**shape))
		ax.plot(bins, y, linewidth=2, color='r')
		ax.set_xlabel("Axon Diameters (um)")
		plt.show()

	# Find the axon with the lowest radius
	try:
		ws.raxmin = r[np.where(r != 0)].min()
	except ValueError:
		# No axons
		ws.raxmin = 0.

	# Number of points in the nerve's contour
	ws.npc = len(contour_nerve['Nerve'])

	# An array with the angular positions of the points on the nerve's contour
	# I will need this for the stimulation with cuff electrodes
	ws.contour_angles = np.zeros(ws.npc)
	for i, (xc, yc) in enumerate(contour_nerve['Nerve']):
		# ws.contour_angles[i] = geo.pts_angle(nvt.centroid, (xc, yc))
		ws.contour_angles[i] = geo.pts_angle(nvt.circumcenter, (xc, yc))

	# Save in the workspace
	ws.r = nvt.r
	ws.contour_nerve = contour_nerve

	# z-axis
	ws.z = np.linspace(0., ws.length, ws.nseg)

	
	# Create figure to visualise the tessellation process
	if ws.settings["graphics"]:
		fig, ax = plt.subplots()
		ws.ax = ax
		ws.nvt.pd.draw(ax)
		# ws.nvt.draw_axons_and_points(ax, facecolor='grey')
		ws.nvt.draw_axons_and_points(ax, facecolor='grey', lw=0)
		ws.nvt.draw_contours(ax, c='r')
		if False:
			for i, (x, y) in enumerate(zip(ws.nvt.x, ws.nvt.y)):
				ax.text(x, y, i)
				ax.text(x, y, ws.nvt.models[i][0], color='r')
				if ws.nvt.models[i] == 'gaines_sensory':
					ax.scatter(x, y, c='b', s=np.pi * ws.nvt.r[i]**2, zorder=1e9)
				elif ws.nvt.models[i] == 'MRG':
					ax.scatter(x, y, c='r', s=np.pi * ws.nvt.r[i]**2, zorder=1e9)
			# Draw the contour from the PD
			# ws.nvt.pd.draw(ax, colour='g', linestyle='-', linewidth=3)
		if False:
			ws.nvt.pd.draw(ax, colour='k', linestyle='-', linewidth=0, values='precomputed', precompv=ws.nvt.free_areas, alpha=0.2)
		ax.set_aspect('equal')
		ax.legend()
		plt.show()

	ws.log("Nerve has been created")

def add_sec(sectype, l, sections_, xstart):
	""" Add a section 'sectype' of length 'l' to the list """
	# global sections_, xstart
	sections_.append(sections[sectype](l, xstart))
	xstart += l
	return sections_, xstart

def build_sequence(axon):
	""" Builds the sequence of sections for an axon.
	Completely modified and adapted to any type of basic sequence 
	on 6 January 2019 """

	# Variables shortnames
	start = axon.properties['start']
	L = ws.length
	basic_sequence = axon.properties['basic sequence']
	lengths = {}
	for key, subkey in axon.properties['section lengths'].items():
		lengths[key] = axon.properties[subkey]

	# Actual sequence
	sequence = []
	# Positions where the sections start
	zz = []

	# Iterate while not finished
	z = 0
	i = 0
	finished = False
	while not finished:
		# Select the type of section
		section = basic_sequence[i % len(basic_sequence)]
		# Append it to the sequence
		sequence.append(section)
		zz.append(z)
		# Update the position
		z += lengths[section]
		# Update counter and check if we finished
		i += 1
		finished = (z >= L + start)

	# zz as an array
	zz = np.array(zz)

	# Identify the first section according to the start point
	first_section_index = np.where(zz > start)[0][0] - 1
	# Update the sequence and zz
	sequence_v2 = sequence[first_section_index:]
	zz_v2 = zz[first_section_index:]

	# Finally, update zz by cutting the first section from the start point
	zz_v3 = zz_v2.copy()
	zz_v3[0] = start

	# Once all this is done, create an array with the actual lengths for 
	# all sections
	actual_lengths = zz_v3[1:] - zz_v3[:-1]
	actual_lengths = np.append(actual_lengths, L + start - zz_v3[-1])

	# Do something: Make the sequence be objects
	sections = []
	for sec, l, zvalue in zip(sequence_v2, actual_lengths, zz_v3):
		sections.append(Section(sec, l, zvalue))

	return sections

def all_vars(axon):
	""" From the sequence, get all the wanted variables """

	# Build the actual sequence of sections and section lengths
	sections = build_sequence(axon)
	
	# Section counters, classified by section types
	section_counter = {}
	for key in ws.axonmodel_settings[axon.model]['section types']:
		section_counter[key] = 0

	# List of lenghts for each section of each section type
	# Just create the lists for now
	lengths = {}
	for key in ws.axonmodel_settings[axon.model]['section types']:
		lengths[key] = []
		
	# Add values to the counters, length lists 
	for section in sections:
		section_counter[section.sectype] += 1
		lengths[section.sectype].append(section.length)

	# Lengths: list to array
	for key, values in lengths.items():
		lengths[key] = np.array(values)

	# Largest n
	n_max = max(section_counter.values())

	return sections, section_counter, n_max, lengths

def locate(sclns, z):
	""" Locate the section and segment (pos) on cable where z 
	lies on """
	pos = 0.
	for i, l in enumerate(sclns):
		if pos <= z <= pos + l:
			return i, (z - pos) / l
		pos += l

def zprofile(cell):
	""" Get the z-positions of all the segments of a cell """
	z = []
	z_accum = 0
	for sec in list(cell.all):
		l = sec.L
		for seg in list(sec.allseg())[1:-1]:
			z.append(z_accum + seg.x * l)
		z_accum += l
	return np.array(z)

def lambda_f(f, d, ra, cm):
	return 1e5 * np.sqrt(d / (4. * np.pi * f * ra * cm))

def nseg_dlambda_rule(l, d, ra, cm):
	return int((l / (0.1 * lambda_f(100., d, ra, cm)) + 0.9) / 2) * 2 + 1

def generate_contours():
	""" Generate a nerve with fascicles inside, all given by the  
	settings in anatomy.json """

	# From the settings, select only the part regarding the generation 
	# of the cross-section
	settings = ws.anatomy_settings["cross-section"]["contours generation"]

	##########################
	# Nerve

	# dimensions and centre
	rn = 0.5 * settings["nerve"]["diameter"]
	c = tuple(settings["nerve"]["center"])

	# Angles
	n = 360
	angles = 2. * np.pi * np.arange(n) / n

	# Points
	x = c[0] + rn * np.cos(angles)
	y = c[1] + rn * np.sin(angles)

	# Store that into a larger array
	contours = np.array([x, y]).T

	# File name
	fname = os.path.join(ws.folders["data/saved"], settings["filename"])

	# Write into csv file
	with open(fname, "w") as f:
		fw = csv.writer(f)
		fw.writerow(["Nerve"])
		for (x, y) in contours:
			fw.writerow([x, y])

	cnerve = contours.copy()

	##########################
	# Fasicles

	# Number of fascicles
	nfas = settings["fascicles"]["number"]
	fdiams = settings["fascicles"]["diameters"]
	fcenters = settings["fascicles"]["centers"]

	# Points
	x = np.zeros((nfas, n))
	y = np.zeros((nfas, n))
	contours = np.zeros((nfas, n, 2))
	for i in range(nfas):
		x[i] = fcenters[i][0] + 0.5 * fdiams[i] * np.cos(angles)
		y[i] = fcenters[i][1] + 0.5 * fdiams[i] * np.sin(angles)
		# Store that into a larger array
		contours[i] = np.array([x[i], y[i]]).T

	##########################
	# Write into csv file
	with open(fname, "a") as f:
		fw = csv.writer(f)
		for i in range(nfas):
			fw.writerow(["Fascicle_%02i"%i])
			for (xx, yy) in contours[i]:
				fw.writerow([xx, yy])

def fill_contours():
	""" Fill the contours with axons and NAELC as specified by 
	anatomy.json"""

	# Load settings
	settings = ws.anatomy_settings["cross-section"]["fibers distribution"]
	apsettings = settings["axon placements"]
	packsettings = apsettings["packing"]

	# AXON PACKING AND TESSELLATION
	# Parameters
	params = {}

	# Fascicle filling

	# Minimum separation
	params['min_sep'] = apsettings["minimum separation"]
	# Triangulation
	# Apprx. max. no. of points or triangles (not counting axons)
	params['mnd'] = settings["min NAELC density"]
	# Reduction of the number of points in the contours
	params['c_reduction'] = ws.anatomy_settings["cross-section"]["contours point reduction"]
	# Type of axon packing
	params['packing_type'] = packsettings["type"]
	params['avg_r'] = packsettings["avg. radius"]
	params['gamma_shape'] = packsettings["gamma shape"]

	# Number of axons, locations (centers) and radii
	params['numberofaxons'] = apsettings["number"]
	params['locations'] = apsettings["locations"]
	params['radii'] = apsettings["radii"]
	params['models'] = apsettings["models"]

	# Minimum and maximum radii
	params['rmin'] = packsettings["min radius"]
	params['rmax'] = packsettings["max radius"]
	# Tolerance for trying to pack circles
	params['circp_tol'] = packsettings["max iterations"]
	# Maximum number of axons per fascicle
	params['max_axs_pf'] = packsettings["max axons per fascicle"]


	# Build tessellated nerve
	# ws.nvt = NerveTess()
	print('generation of ws.nvt')
	ws.nvt.build(params)

	# Axons Start Point Along the z-axis
	# 'start' point. This is used to randomly locate the 
	# axon over the z-axis so that it doesn't necessarily start
	# from a node of Ranvier
	ws.nvt.start_positions = np.zeros_like(ws.nvt.x)
	node_misalignement = ws.anatomy_settings['axons']['myelinated']['nodes misalignement']
	
	for i in range(len(ws.nvt.x)):
		if ws.nvt.cables[i] == "Axon":
			ws.nvt.start_positions[i] = np.random.uniform(0, node_misalignement)