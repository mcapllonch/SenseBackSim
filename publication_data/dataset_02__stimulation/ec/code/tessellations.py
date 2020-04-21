"""
Module to create tessellations

This algorithm successfully crops the diagram and is capable of dealing with
vertical lines (8 June 2018)

"""
import warnings
import numpy as np
import scipy.spatial as sp
import triangle as tgl
from collections import OrderedDict
import matplotlib.patches as mp
from matplotlib.collections import PatchCollection
from matplotlib.ticker import LogFormatter
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt

import algebra as alg
import geometry as geo
import tools
import workspace as ws

# Ignore warnings
warnings.filterwarnings('ignore')

class PowerDiagram():
	"""
	Class for a Power diagram, or a Laguerre-Voronoi diagram
	"""
	def __init__(self, x, y, r, contour=None):
		self.x = x
		self.y = y
		self.r = r
		self.circ_areas = np.pi * r ** 2
		self.nc = x.size
		self.contour = contour
		self.any_errors = False
		# self.build()
		# self.get_polygons()
		
	def build(self):
		"""
		Build it
		"""
		
		ws.log("Building the power diagram")

		x, y, r = self.x, self.y, self.r
		# Generating circles
		gc = GCircles(x, y, r)
		gc.triangulate()

		# First iteration
		pairs, trios, segments, pairs_trios, pairs_oprs, pairs_opts, \
			pairs_allpts, pvpts, lines = build_laguerre(gc, gc.trios)
		# First correction
		new_trios, new_pairs, rmv_trios, rmv_pairs, finished = \
			one_correction(pairs_oprs, pairs_opts, pairs_trios, segments)
		trios = update_trios(new_trios, rmv_trios, trios)
		# Iterate until there's no more bad crossings
		# Iteration counter
		corrcount = 0
		while not finished:
			pairs, trios, segments, pairs_trios, pairs_oprs, pairs_opts, \
				pairs_allpts, pvpts, lines = build_laguerre(gc, trios)
			new_trios, new_pairs, rmv_trios, rmv_pairs, finished = \
				one_correction(pairs_oprs, pairs_opts, pairs_trios, segments)
			trios = update_trios(new_trios, rmv_trios, trios)
			# Go fix the only-two-connections issue
			poss_lost = []
			poss_lost, rmv_trios, new_trios = fix_biconnections(poss_lost, 
				pairs, gc, pairs_trios)
			if len(poss_lost) > 0:
				trios = update_trios(new_trios, rmv_trios, trios)
				pairs, trios, segments, pairs_trios, pairs_oprs, pairs_opts, \
					pairs_allpts, pvpts, lines = build_laguerre(gc, trios)
				# 'This should be fixed now'
				whoslost = poss_lost[0]
			corrcount += 1
			enough = corrcount == gc.nc / 4
			
			# Check if this is finished
			finished = finished or enough
			ws.log('Number of corrections: %i'%corrcount)
			if enough:
				error_msg = 'TESSELLATION ERROR: Something went ' \
				 + 'wrong with this set of circles\n'\
				 + 'This is probably due to a bug in the algorithm\n' \
				 + 'Please try a different set of circles\n'
				ws.log(error_msg)
				self.any_errors = True
				print('TESSELLATION ERROR')
				return
		ws.log('Correction iterations: %i'%corrcount)

		# Check if two connections cross
		for ip, pair1 in enumerate(pairs):
			a = segments[pair1]
			for pair2 in pairs[ip+1:]:
				b = segments[pair2]
				if isinstance(a, geo.Segment) and isinstance(b, geo.Segment):
					crossing, _ = geo.crossing_segments(a, b)
					overlapping = geo.overlapping_segments(a, b)
					if crossing or overlapping:
						ws.log('TESSELLATION ERROR: An error was found for pairs: %s, %s'%(str(pair1), str(pair2)))
						ws.log('This is due to a bug in the algorithm')
						ws.log('Please try a different set of circles')
						self.any_errors = True
						print('TESSELLATION ERROR')
						return
				else:
					pass
		
		# If finished, crop the diagram using the contour
		if self.contour is not None:

			# print('')
			# print('segments:')
			# for k, v in segments.items():
			# 	print(k, v)

			# Crop it
			segments = crop_diagram(gc, pairs, segments, pvpts, lines, self.contour)

		# Finally, remove the segments which are not geo.Segment 
		# instances
		for k, seg in segments.items():
			if not isinstance(seg, geo.Segment):
				if len(seg) != 2:
					segments[k] = None
				if len(seg) == 2:
					segments[k] = geo.Segment(seg)

		self.pairs = pairs
		self.segments = segments
		self.trios = trios

		self.get_polygons()

	def build_preexisting(self, pairs, segments):
		""" Set up a diagram with preset properties 
		Created on 4 October 2018 """
		pairs = list(pairs)
		self.pairs = pairs
		self.segments = segments

		trios = []
		for i, p1 in enumerate(pairs):
			for p2 in pairs[i+1:]:
				s1 = set(p1)
				s2 = set(p2)
				# If there is only three points involved between the two 
				# pairs, this is clearly a potential trio
				union = s1.union(s2)
				if len(union) == 3:
					# But not all trios are real in our diagram, so we need to 
					# check that this trio exists really
					# For this, check that the two non-common points, e.g. 
					# (2, 3) if the pairs were (1, 2) and (1, 3) are a listed 
					# pair
					if tuple(s1.symmetric_difference(s2)) in pairs:
						trios.append(tuple(union))
		self.trios = trios
		self.get_polygons()

	def get_polygons(self):
		"""
		Create the polygons for each cell
		"""
		
		segments = self.segments
		x, y = self.x, self.y
		nc = self.nc
		cellverts = OrderedDict()
		polygons = OrderedDict()
		polg_areas = np.zeros(nc)
		# Add vertices to the cells
		for (i, j), seg in segments.items():
			try:
				a, b = seg.a, seg.b
			except AttributeError:
				# This means that some segment was [None]
				print("None segment for (%i, %i)"%(i, j))
				pass
			else:
				cellverts = tools.append_items(cellverts, i, [a, b])
				cellverts = tools.append_items(cellverts, j, [a, b])
		
		# Set vertices
		for i in range(nc):
			points = np.array(cellverts[i])
			# Reshape
			if points.shape[-1] == 1:
				points = points.reshape(points.shape[:-1])
			# Set them
			points = np.unique(points, axis=0)
			# Sort points by angle
			cellpoint = (x[i], y[i])

			# Sort them by angle. Use some averaged centroid rather than
			# the cellpoint itself. Because if cellpoint is included, it
			# may be confusing and the sorting can go wrong
			# This is OK, since the polygon is always convex
			points = geo.sort_by_angle(points)
			cellverts[i] = points

			# Create polygon
			pol = geo.Polygon(points)

			# Check that the polygon contains the point or this is on 
			# the contour, at least. If not, add it to the polygon
			polccell = pol.contains_point(cellpoint)

			# if False:
			if not polccell:
				points = np.append(points.tolist(), [list(cellpoint)], axis=0)
				points = geo.sort_by_angle(points)
				cellverts[i] = points
				pol = geo.Polygon(points)

			try:
				polg_areas[i] = pol.area
			except AttributeError:
				# This is not a polygon
				pass
			else:
				# If it is, save it
				polygons[i] = pol

		# Store in attributes
		self.cellverts = cellverts
		self.polygons = polygons
		self.polg_areas = polg_areas
		self.free_areas = polg_areas - self.circ_areas

	def draw(
			self, 
			ax, 
			draw_polygons=True, 
			colour='k', 
			line_alpha=1., 
			fill_alpha=1., 
			linestyle='-', 
			linewidth=1.,
			markers=None, 
			values=None, 
			precompv=None, 
			vmin=None, 
			vmax=None, 
			title=None, 
			cmap=None, 
			cblabel=None, 
			logscale=False, 
			allowed_cbticklabels=(None,), 
			extend='neither', 
			zorder=0
		):
		"""
		Draw the diagram
		"""

		drawvalues = values is not None

		if drawvalues:
			patches = []
			colours = np.zeros(self.nc)
	
		for i, pol in self.polygons.items():

			# Colours and filled polygons
			if drawvalues:
				filled_pol = mp.Polygon(self.cellverts[i])
				patches.append(filled_pol)
				if values == 'areas':
					colours[i] = pol.area
				if values == 'inv_areas':
					colours[i] = 1. / self.free_areas[i]
				if values == 'free_areas':
					colours[i] = self.free_areas[i]
				if values == 'precomputed':
					if logscale:
						colours[i] = np.abs(precompv[i])
					else:
						colours[i] = precompv[i]

				# Colour for the polygon contours
				if colour == 'same_as_facecolors':
					pc_ = colours[i]
			else:
				pc_ = colour

			# Draw the polygon
			if draw_polygons:
				pol.draw(ax=ax, colour=pc_, alpha=line_alpha,
					linestyle=linestyle, linewidth=linewidth, markers=markers, 
					zorder=zorder)

		# Line colours
		if drawvalues:
			if vmin is None:
				vmin = colours.min()
			if vmax is None:
				vmax = colours.max()
			if logscale:
				vmin = np.abs(vmin)
				vmax = np.abs(vmax)
				norm = LogNorm(vmin=vmin, vmax=vmax)
			else:
				norm = Normalize(vmin=vmin, vmax=vmax)
			m = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
			line_colours = m.to_rgba(colours)
			print(line_colours)

		if drawvalues:

			collection = PatchCollection(patches, cmap=cmap, norm=norm)
			ax.add_collection(collection)
			collection.set_alpha(fill_alpha)
			collection.set_array(colours)
			if logscale:
				collection.set_clim(vmin=vmin, vmax=vmax)
			else:
				collection.set_clim(vmin=0., vmax=vmax)
			if colour == 'same_as_facecolors':
				collection.set_edgecolors(line_colours)
			else:
				collection.set_edgecolor(colour)
			collection.set_linewidth(linewidth)
	
		if values == 'precomputed':

			collection.set_clim(vmin=vmin, vmax=vmax)

			# Color bar
			if logscale:
				l_f = LogFormatter(10, labelOnlyBase=False)
				cb = plt.colorbar(collection, ax=ax, orientation='vertical', format=l_f, extend=extend)
				# Set minor ticks
				# We need to nomalize the tick locations so that they're in the range from 0-1...
				# minorticks = p.norm(np.arange(1, 10, 2))
				# Minimum and maximum exponents
				min_exp = int(np.log10(min(vmin, vmax)))
				max_exp = int(np.log10(max(vmin, vmax))) + 1
				# Ticks for the colorbar in case it's logscale
				minorticks = 10. ** np.arange(min_exp, max_exp + 1, 1)
				minorticks_ = minorticks.tolist()
				for i in range(2, 10, 1):
					minorticks_ = minorticks_ + (float(i) * minorticks).tolist()
				minorticks_.sort()
				minorticks = np.array(minorticks_)
				minorticks = minorticks[np.where(minorticks >= vmin)]
				minorticks = minorticks[np.where(minorticks <= vmax)]
				print(minorticks)
				cb.set_ticks(minorticks)
				cb.set_ticklabels([str(int(x)) if x in allowed_cbticklabels else '' for x in minorticks])
			else:
				cb = plt.colorbar(collection, ax=ax, orientation='vertical', extend=extend)
			cb.set_label(cblabel)
		if title is not None:
			ax.set_title(title)

		ax.set_xlabel(r'$\rm x$ ($\rm \mu m$)')
		ax.set_ylabel(r'$\rm y$ ($\rm \mu m$)')

class PointsEnsemble():
	""" Cloud of points """
	def __init__(self, x, y):
		self.x = x
		self.y = y
		self.points = np.array([self.x, self.y]).T
		self.nc = len(self.x)
		# Matrices
		self.dx = alg.compl_mat(x)
		self.dy = alg.compl_mat(y)

	def triangulate(self):
		""" Get the Delaunay triangulation of the points """
		# Triangulation
		if self.nc == 2:
			self.pairs = np.array([np.arange(2)])
			self.trios = np.array([np.arange(2)])
			self.triangles = {}
		else:
			tri = sp.Delaunay(self.points)
			trios = tri.simplices.copy()
			self.triangulation = tri
			self.update_triangulation(trios)

	def update_triangulation(self, trios):
		""" Update the triangulation for a new arrangement of trios """

		pairs = []
		# Triangles
		triangles = {}
		for trio in trios:
			trio = tuple(sorted(trio))
			triangle = geo.Triangle(self.points, indices=trio, fas=None)
			triangles[trio] = triangle
			# Compute the resistors for each side of the triangle
			for i, pair in enumerate(triangle.pairs):
				# We do not want to repeat pairs, right?
				pair = tuple(sorted(pair))
				if pair not in pairs:
					pairs.append(pair)
		pairs = np.array(pairs)
		self.pairs = pairs
		self.triangles = triangles
		trios = np.sort(trios)
		self.trios = trios
		
class GCircles(PointsEnsemble):
	""" Ensemble of generating circles in the system """
	def __init__(self, x, y, r):
		PointsEnsemble.__init__(self, x, y)
		self.r = r
		self.nc = self.x.size
		# Vectors
		# Sigmas
		self.sigma = r * r - (x * x + y * y)
		# Matrices
		self.ds = alg.compl_mat(self.sigma)

class LinesEnsemble():
	""" Ensemble of lines """
	def __init__(self):
		pass

	def cutoffpairs(self):
		""" Reduce the size of the matrices according to the pairing in the 
		triangulation """
		pi, pj = self.pairs.T
		self.m = self.m[pi, pj]
		self.n = self.n[pi, pj]
		self.pi, self.pj = pi, pj

	def vertical_to_mbform(self):
		""" Get all the lines in point-slope form, so I can deal with
		vertical lines in other calculations """

		# Provisional b-points (for the point-slope form)
		self.bx = np.zeros_like(self.n)
		self.by = self.n
		# Do the necessary with vertical lines
		pi, pj = self.pi, self.pj
		dxvec = self.pe.dx[pi, pj]
		dyvec = self.pe.dy[pi, pj]
		dsvec = self.pe.ds[pi, pj]
		ivlines = np.where(dyvec == 0.)
		# Find the point bx
		bx = -0.5 * dsvec[ivlines] / dxvec[ivlines]
		# by
		by = self.pe.y[self.pairs[ivlines, 0]]
		self.bx[ivlines] = bx
		self.by[ivlines] = by
				
	def dict(self):
		""" Create line instances in the geometry module so that they can
		be drawn on a figure """
		bx, by = self.bx, self.by
		pairs, n, m = self.pairs, self.n, self.m
		lines = {}
		for i, pair in enumerate(pairs):
			bxi, byi = bx[i], by[i]
			lines[tuple(pair)] = geo.StraightLine((bxi, byi), m[i])
		# Store them as attributes
		self.lines = lines
		self.bxvec = bx
		self.byvec = by
		self.bpoints = np.array([bx, by]).T
	
class LinesConnectingCenters(LinesEnsemble):
	""" Ensemble of lines connecting the centers of an ensemble of 
	circles """
	def __init__(self, gc):
		self.pe = gc
		self.nc = self.pe.nc
		self.pairs = self.pe.pairs
		self.nl = len(self.pairs)
		LinesEnsemble.__init__(self)

	def build(self):
		""" Build the lines """

		x, y = self.pe.x, self.pe.y
		dx, dy = self.pe.dx, self.pe.dy
		m = dy / dx
		n = alg.vectomat(y).T - m * alg.vectomat(x).T

		# Store everything as attributes
		self.m = m
		self.n = n

		self.cutoffpairs()

		# Deal with vertical lines
		self.vertical_to_mbform()

		# Points on the lines
		# Indices of (the first of) the points each line belongs to
		which = self.pairs[:, 0]
		self.bx = x[which]
		self.by = y[which]

		# Create its instances in the geometry module
		self.dict()
				
class PowerLines(LinesEnsemble):
	""" Ensemble of power lines between circles """
	def __init__(self, gc, lcc):
		self.pe = gc
		self.lcc = lcc
		self.nc = gc.nc
		self.pairs = self.pe.pairs
		self.nl = len(self.pairs)
		LinesEnsemble.__init__(self)

	def build(self):
		""" Build the lines """
		dx, dy, ds = self.pe.dx, self.pe.dy, self.pe.ds
		# Slopes
		m = - dx / dy
		# Intercept (y value)
		n = -0.5 * ds / dy
		# Store everything as attributes
		self.m = m
		self.n = n

		self.cutoffpairs()

		# Deal with vertical lines
		self.vertical_to_mbform()

		# Create its instances in the geometry module
		self.dict()

		# Find points in the lines for geo.StraightLine
		# Find the power-line mid-points between pairs
		mps = geo.intersections_many_lines([self, self.lcc])
		xm, ym = mps
		# Only interested in the points that involve both ensembles
		xm = xm[:self.nl, self.nl:]
		ym = ym[:self.nl, self.nl:]
		xm = np.diag(xm)
		ym = np.diag(ym)
		self.bx = xm
		self.by = ym

		# Create its instances in the geometry module
		self.dict()

def remove_repeated_points(x, y):
	""" Make a set of the points so there are no repeated points """

	# Use my geometry module, because it may very well be that points
	# which are the same in theory, are not the same in practice

	# Arrange points
	points = np.array([x, y]).T

	# Distances among all points
	# dists = np.hypot(alg.compl_mat(x), alg.compl_mat(y))
	# pairs, dists = alg.one_side(dists)
	pairs, dists = geo.dist_cloud(x, y)
	# Where are the points nearly identical
	where = np.where(dists <= 1.e-9)
	pairs_matched = pairs[where]

	# Remove the second one from each pair of nearly identical points
	allpts = set(pairs.flatten().tolist())
	allmatched = set(pairs_matched.flatten().tolist())
	removed = set(pairs_matched[:, 1].tolist())
	survivors = allpts - removed
	# survivors = tuple(survivors)
	survivors = np.array(list(survivors))
	points = points[survivors]

	x, y = points.T
	# If there is only one point, make x and y be arrays
	if points.ndim == 1:
		x = np.array([x])
		y = np.array([y])

	return x, y

def build_laguerre(gc, trios):
	""" Build the Laguerre diagram from a cloud of points and its given 
	triangulation """

	# Update the triangulation
	gc.update_triangulation(trios)
	trios = gc.trios

	# Lines connecting centers
	lcc = LinesConnectingCenters(gc)
	lcc.build()
	# Power lines
	pl = PowerLines(gc, lcc)
	pl.build()

	# Possible Vertices (PV)
	# Intersections among power lines
	xi, yi = geo.intersections_many_lines([pl])
	lpairs, xi = alg.one_side(xi)
	lpairs, yi = alg.one_side(yi)
	# Determine the group of circles for each possible vertex
	pvnc = len(yi)
	gcpvc = pl.pairs[lpairs]
	gcpvc = gcpvc.reshape((pvnc, 4))
	gcpvc = np.sort(gcpvc)
	lgp = np.array(list(map(lambda l: len(set(l)), gcpvc)))
	keep_these = np.where(lgp == 3)
	xi = xi[keep_these]
	yi = yi[keep_these]
	lpairs = lpairs[keep_these]
	gcpvc = gcpvc[keep_these]
	pvnc = len(xi)
	gcpvc = np.array(list(map(lambda l: sorted(list(set(l))), gcpvc)))
	# Now I have to iterate over trios
	# Coordinates of the possible vertex for each triangle
	xpv = {}
	ypv = {}
	# Segment for each pair
	segments = {}
	# Number of PV points for each pair
	pvpts = {}
	# One (if the pair is on the boundary) or two trios per pair
	pairs_trios = {}
	# Other pairs of the quad whose diagonal is a certain pair
	pairs_oprs = {}
	# Other points of the quad whose diagonal is a certain pair
	pairs_opts = {}
	# All points of the quad whose diagonal is a certain pair
	pairs_allpts = {}
	# Iterate
	for trio in trios:
		trio_extended = np.tile(trio, pvnc).reshape(gcpvc.shape)
		matches = np.equal(gcpvc, trio_extended)
		matches = np.all(matches, axis=1)
		matches = np.where(matches == True)[0]
		matches = np.unique(matches)
		xim = xi[matches]
		yim = yi[matches]

		# Now reduce the duplicated points for this triangle
		if len(xim) > 1:
			xim, yim = remove_repeated_points(xim, yim)


		# Store them
		if len(xim) > 0:
			tt = tuple(trio) 

			xim, yim = xim.mean(), yim.mean()

			pairs_sorted = [tuple(sorted(item)) for item in gc.triangles[tt].pairs]

			# Pairs
			for pair in pairs_sorted:
				other_pairs = [item for item in pairs_sorted if item != pair]
				try:
					segments[pair].append((xim, yim))
					pvpts[pair] += 1
				except KeyError:
					segments[pair] = [(xim, yim)]
					pvpts[pair] = 1
					# For this pair, store the trio and the other pairs
					pairs_trios[pair] = [tt]
					pairs_oprs[pair] = other_pairs
				else:
					pairs_trios[pair].append(tt)
					pairs_oprs[pair] = pairs_oprs[pair] + other_pairs
				pa = set(np.array(pairs_oprs[pair]).flatten().tolist() \
											+ list(pair))
				pairs_allpts[pair] = list(pa)
				pairs_opts[pair] = sorted(list(pa - set(pair)))

	# List pairs
	pairs = list(pairs_oprs.keys())

	# for k, v in pvpts.items():
	# 	print(k, v)

	# Get far points
	segments = get_far_points(gc, segments, pvpts, pl.lines)

	# Update segments dictionary
	for k, seg in segments.items():
		# if len(seg) != 2:
		if not isinstance(seg, geo.Segment):
			if len(seg) < 2:
				segments[k] = seg
				# pass
			if len(seg) == 2:
				segments[k] = geo.Segment(seg)
	# Return
	return pairs, trios, segments, pairs_trios, pairs_oprs, pairs_opts, pairs_allpts, pvpts, pl.lines

def get_far_points(gc, segments, pvpts, lines):
	"""
	Get the points being far away
	"""


	# Pairs that only have one point (on the boundary)
	iwhobdr = np.where(np.array(list(pvpts.values())) == 1)
	whobdr = np.array(list(pvpts.keys()))[iwhobdr]
	# Give them a second point very far away
	# Calculate limits
	xmin, xmax = gc.x.min(), gc.x.max()
	ymin, ymax = gc.y.min(), gc.y.max()
	dx = xmax - xmin
	dy = ymax - ymin
	xfar_left = xmin - 1e3 * dx
	xfar_rght = xmax + 1e3 * dx
	yfar_bot = ymin - 1e3 * dy
	yfar_top = ymax + 1e3 * dy

	for pair in whobdr:
		# Format the pair to amke it a tuple
		pair = tuple(pair.tolist())
		# Get point and clean it
		try:
			point = tuple([float(x) for x in segments[pair][0]])
		except TypeError:
			point = segments[pair].a
		if point is not None:

			powl = lines[pair]

			# If it's a vertical line, do it differently
			slope = powl.slope
			vline = (slope == np.inf) or (slope == -np.inf) or (np.isnan(slope))

			if vline:
				xfar = point[0]
				ypnt = point[1]
				farpt_up = geo.Segment([(xfar, ypnt), (xfar, yfar_top)])
				farpt_dn = geo.Segment([(xfar, ypnt), (xfar, yfar_bot)])
				possible_segs = [farpt_dn, farpt_up]
			else:

				yfar_left = powl.equation(xfar_left)
				yfar_rght = powl.equation(xfar_rght)
				
				# Only one of the segments does not intersect with any other
				sleft = geo.Segment([point, (xfar_left, yfar_left)])
				srght = geo.Segment([point, (xfar_rght, yfar_rght)])
				possible_segs = [sleft, srght]

			# Try the left-hand segment and see if that's the good one
			valid_segment = np.array([False, False])
			for i, tryseg in enumerate(possible_segs):
				this_one_not_valid = False
				for othpair, othseg in segments.items():
					if (pair != othpair) and isinstance(othseg, geo.Segment):
						crossing, crosspoint = geo.crossing_segments(tryseg, othseg)
						overlapping = geo.overlapping_segments(tryseg, othseg)
						this_one_not_valid = crossing or overlapping
						if this_one_not_valid:
							# Not this one
							break
				if not this_one_not_valid:
					valid_segment[i] = True

			if (valid_segment == True).all():
				# Choose the shortest one
				s0 = possible_segs[0]
				s1 = possible_segs[1]
				if s0.length < s1.length:
					segments[pair] = s0
				else:
					segments[pair] = s1
			elif (valid_segment == True).any():
				segments[pair] = possible_segs[np.where(valid_segment == True)[0][0]]
			else:
				print('No far point could be given to:', pair)

	return segments
				
def crop_diagram(gc, pairs, segments, pvpts, lines, contour):
	"""
	Crop the diagram using a contour
	"""

	segments = get_far_points(gc, segments, pvpts, lines)

	# Calculate limits
	cont_xy = contour['vertices']
	# cont_sg = cont_xy[contour['segments']]
	contour_polygon = geo.Polygon(cont_xy)
	x, y = np.array(cont_xy).T
	xmin, xmax = x.min(), x.max()
	ymin, ymax = y.min(), y.max()
	dx = xmax - xmin
	dy = ymax - ymin

	# Now get the intersection between a pair's segment and the contour
	for pair in pairs:
		
		seg = segments[pair]
		if isinstance(seg, geo.Segment):

			# First off, if the segment is completely outside the contour, 
			# remove it
			a, b = seg.a, seg.b
			discard = (not contour_polygon.contains_point(a)) and (not contour_polygon.contains_point(b))
			if not discard:
				for cs in contour_polygon.segments.values():
					# xyint = geo.intersection_segments(seg, geo.Segment(cs))
					xyint = geo.intersection_segments(seg, cs)
					
					if xyint != (np.nan, np.nan):
						# We found an intersection between the segment and the
						# contour
						# Now the point to be replaced is the one falling outside
						if contour_polygon.contains_point(a):
							change = 1
						elif contour_polygon.contains_point(b):
							change = 0
						# the contour
						seg.change_point(change, xyint)
						# We're done with this loop
						break

				# Sanity Check after cutting it
				# If any segment is even partially outside the contour, remove it
				# I've seen this happen when a segment has one point 
				# outside and one right on the contour, for which no 
				# intersection is found and hence it's left as is

				# Check if any point is detected to be outside the contour (being 
				# on it may yield True, unfortunately), so in that case,
				check = (not contour_polygon.contains_point(a)) or (not contour_polygon.contains_point(b))
				if check:
					# roughly check if it's really outside. In that 
					# case, it may be really really far (thus be a 
					# long segment)
					discard = seg.length > 10. * np.hypot(dx, dy)
					if discard:
						ws.log('Discarding %s'%str(pair))
						ws.log(str(contour_polygon.contains_point(seg.a)))
						ws.log(str(seg.points.flatten()))
						ws.log(str(seg))
						ws.log('length = %s, dx = %s, dy = %s'%(str(seg.length), str(dx), str(dy)))
						ws.log("")
						segments[pair] = [None]
					else:
						# Add the segment
						segments[pair] = seg
				else:
					# Add the segment
					segments[pair] = seg
			else:
				# This is not good doing anymore
				# segments[pair] = [None]
				pass

	# List pairs again
	pairs = list(segments.keys())

	return segments

def one_correction(pairs_oprs, pairs_opts, pairs_trios, segments):
	""" Perform one single trio change """

	# New trios
	new_trios = []
	# New pairs
	new_pairs = []

	# Removed trios
	rmv_trios = []
	# Removed pairs
	rmv_pairs = []

	# Iterate over pairs to check crossing segments
	pairs = list(pairs_trios.keys())[:]
	for pair in pairs:
		other_pairs = pairs_oprs[pair]
		for i, op1 in enumerate(other_pairs):
			a = segments[op1]
			for op2 in other_pairs[i+1:]:
				b = segments[op2]
				# if (len(a) == 2) & (len(b) == 2):
				# if True:
				if isinstance(a, geo.Segment) and isinstance(b, geo.Segment):
					crossing, _ = geo.crossing_segments(a, b)
					if crossing:
						# Change the connection: use the other diagonal of the
						# quad
						# Remove bad pair
						rmv_pairs.append(pair)
						# Update with new pair
						new_pair = tuple(pairs_opts[pair])
						new_pairs.append(new_pair)
						
						# Update its segment
						xyint = geo.intersection_segments(a, b)
						try:
							isgeoSeg = isinstance(segments[new_pair], geo.Segment)
						except KeyError:
							segments = tools.append_items(segments, new_pair, [xyint])
						else:
							if not isgeoSeg:
								segments = tools.append_items(segments, new_pair, [xyint])
						# if not isinstance(segments[new_pair], geo.Segment):
						# try:
						#	segments[new_pair].append(xyint)
						# except KeyError:
						#	segments[new_pair] = [xyint]

						# Update the trios
						# Remove bad trios
						for trio in pairs_trios[pair]:
							rmv_trios.append(trio)
						# Update with new trios
						for p in pair:
							new_trio = tuple(sorted(list(new_pair) + [p]))
							new_trios.append(new_trio)

						# Update segments dictionary
						for k, seg in segments.items():
							if not isinstance(seg, geo.Segment):
								# if len(seg) != 2:
								if len(seg) < 2:
									segments[k] = seg
								if len(seg) == 2:
									segments[k] = geo.Segment(seg)

						# The important thing here
						return new_trios, new_pairs, rmv_trios, rmv_pairs, False

	return new_trios, new_pairs, rmv_trios, rmv_pairs, True

def update_trios(new_trios, rmv_trios, trios):
	""" Update the information according to the correction """
	# Set new_trios
	new_trios = list(set(new_trios))
	# Update the lists of pairs and trios
	trios = [item for item in trios if tuple(item) not in rmv_trios]
	trios = trios + new_trios
	# To tuple
	trios = tuple(map(tuple, trios))
	# Re-set them again
	trios = list(set(trios))
	trios = np.array(trios)
	return trios

def fix_biconnections(poss_lost, pairs, gc, pairs_trios):
	"""
	Tell me which circle has only 2 connections and what to do about it
	"""
	x, y = gc.x, gc.y
	rmv_trios = []
	new_trios = []
	# Dictionary of neighbours
	circs_neighbours = OrderedDict()
	for i in range(gc.nc):
		circs_neighbours[i] = []
	# Find connections
	for (i, j) in pairs:
		circs_neighbours[i].append(j)
		circs_neighbours[j].append(i)
	# Set them
	for i in range(gc.nc):
		ngh = sorted(circs_neighbours[i])
		circs_neighbours[i] = list(set(ngh))
		if len(ngh) < 3:
			# This circle has less than 3 connections
			# Its connections are:
			tngh = tuple(ngh)
			# Their triangles are:
			alltriosngh = pairs_trios[tngh]
			# The triangles that dont include it are:
			notin = [item for item in alltriosngh if i not in item]
			# Find the triangle where it is contained in
			for trio in notin:
				inthisone = gc.triangles[trio].plpol.contains_point((x[i], y[i]))
				if inthisone:
					# If it is, this means that this circle is not on the boundary
					# of the hull, so the issue must be fixed up
					poss_lost.append(i)
					# Find its new connection
					newcon = list(set(trio) - set(tngh))[0]
					# The two new triangles are:
					new_trios = [tuple(sorted([i, newcon, item])) for item in tngh]
					# And this triangle should be removed
					rmv_trios = [trio]
	return poss_lost, rmv_trios, new_trios