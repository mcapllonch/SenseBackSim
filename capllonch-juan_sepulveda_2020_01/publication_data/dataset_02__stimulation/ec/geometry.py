"""
Classes and functions that facilitate handling geometry
"""

import numpy as np
import planar as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from matplotlib.collections import PatchCollection
from collections import OrderedDict

import workspace as ws
import algebra as alg
import tools as tls



def nearly_equal(a, b, tolerance):
	return np.abs(a - b) < tolerance

def dist(p1, p2):
	"""
	Distance between two points
	If one of the points is given as an array of points with the 
	correct shape, this function actually returns the result for the 
	array
	This shape is (2, npoints)
	"""
	return np.hypot(p2[1] - p1[1], p2[0] - p1[0])

def dist_cloud(x, y):
	""" 
	Mutual distances between a cloud of points 
	"""
	dists = np.hypot(alg.compl_mat(x), alg.compl_mat(y))
	pairs, dists = alg.one_side(dists)
	return pairs, dists

def pts_equal(p1, p2, tolerance):
	"""
	Check if p1 and p2 are the same point
	"""
	return nearly_equal(dist(p1, p2), 0, tolerance)

def get_midpoint(a, b):
	"""
	Get the point lying in the middle of points a and b
	"""
	mpx = 0.5 * (a[0] + b[0])
	mpy = 0.5 * (a[1] + b[1])
	return (mpx, mpy)

def points_in_circle(points, circle):
	"""
	Check which of a array of points is inside a circle
	The array of point needs to have the shape (2, npoints)
	IMPORTANT NOTE: Add a 1.e-10 % of the radius to make up for numerical 
	errors
	"""
	return dist(points, circle.c) <= circle.r * (1. + 1.e-10)

def get_slope(p1, p2):
	"""
	Get the slope of a line passing between p1 and p2
	"""
	try:
		slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
	except ZeroDivisionError:
		slope = np.inf
	return slope

def get_orthoslope(slope):
	"""
	Get the slope of a line orthogonal to the line with a given slope 'slope'
	"""
	try:
		oslope = -1. / slope
	except ZeroDivisionError:
		oslope = np.inf
	return oslope

def angnegtopos(a):
	""" Given a negative angle 'a', return 2pi + a """
	twopi = 2. * np.pi
	n = abs(int(a / twopi)) + 1
	if a < 0:
		return n * twopi + a
	return a

def angle_segment(ref, p):
	"""
	Get angle that p has with respect to a reference point ref over the x-axis
	"""
	dx = p[0] - ref[0]
	dy = p[1] - ref[1]
	angle = np.angle(complex(dx, dy))
	# Turn it positive if it is negative
	angle = angnegtopos(angle)
	return angle

def angle_threepts(a, b, c):
	""" Get the angle at point b made by the segments ab and bc.
	Note: the order of a, b and c determines the value of the angle """
	# angle = angpostoneg(angle_segment(b, c)) - angpostoneg(angle_segment(b, a))
	angle = angle_segment(b, c) - angle_segment(b, a)
	# Turn it positive if it is negative
	return angnegtopos(angle)

def pts_angle(ref, p):
	"""
	Get angle that p has with respect to a reference point ref over the x-axis
	"""
	dx = p[0] - ref[0]
	dy = p[1] - ref[1]
	angle = np.angle(complex(dx, dy))
	if angle < 0:
		angle += 2. * np.pi
	return angle

def sort_by_angle(pts, ref=None):
	"""
	Sort a set of points pts by the angle they have with respect to ref
	"""
	# Reference point. If it is None, calculate it as a centroid of the points
	if ref is None:
		ref = pts.mean(axis=0).tolist()

	angles = []
	for p in pts:
		angles.append(pts_angle(ref, p))
	sorder = np.argsort(angles)

	return pts[sorder]

def intersection_two_lines(sl1, sl2):
	"""
	Calculate the intersection point between two straight lines in the 
	point-slope form
	"""
	ps1, ps2 = sl1.ps, sl2.ps
	p1, s1 = ps1
	p2, s2 = ps2
	r1 = sl1.equation
	r2 = sl2.equation
	# In case they are parallel
	if s1 == s2:
		return np.nan, np.nan
	# In case one of them has an infinite slope
	if np.abs(s1) == np.inf:
		x_int = p1[0]
		y_int = r2(x_int)
	elif np.abs(s2) == np.inf:
		x_int = p2[0]
		y_int = r1(x_int)
	else:
		x_int = (p2[1] - p1[1] + s1 * p1[0] - s2 * p2[0]) / (s1 - s2)
		y_int = r1(x_int)
	return x_int, y_int

def intersections_many_lines(lines_):
	""" Calculate all the intersections among one or various ensembles 
	of lines """

	# Get the vectors of number of lines, slopes and yintercpts from the ensembles
	nl = sum([lines.nl for lines in lines_])
	mvec = np.array([lines.m for lines in lines_]).flatten()
	nvec = np.array([lines.n for lines in lines_]).flatten()
	bxvec = np.array([lines.bxvec for lines in lines_]).flatten()
	byvec = np.array([lines.byvec for lines in lines_]).flatten()

	# Create the matrices to solve the problem
	m = alg.diagv(mvec)
	x = - alg.compl_mat(nvec) / alg.compl_mat(mvec)

	ivlines = np.where(np.abs(mvec) == np.inf)[0]
	nvlines = ivlines.size

	ihlines = np.where(np.abs(mvec) == 0.)[0]
	nhlines = ihlines.size

	# Modify values of x
	# Rows
	bx_rep_rows = np.repeat(bxvec[ivlines], nl).reshape((nvlines, nl))
	x[ivlines] = bx_rep_rows[:]

	# Columns
	bx_rep_cols = np.tile(bxvec[ivlines], nl).reshape((nl, nvlines))
	x[:, ivlines] = bx_rep_cols[:]

	# Infinities to nans
	x[np.where(np.abs(x) == np.inf)] = np.nan
	x = alg.make_sym_rem_nans(x)

	# Wherever x = nan, place any dummy number (it won't affect the 
	# theoretical result). This will avoid y being all nan
	xx = x.copy()
	xx[np.where(np.isnan(x))] = 0.

	# Compute y
	y = np.dot(m, xx) + alg.vectomat(nvec)

	# Modify values of y
	by_rep_rows = np.repeat(byvec[ihlines], nl).reshape((nhlines, nl))
	y[ihlines] = by_rep_rows[:]

	# Columns
	by_rep_cols = np.tile(byvec[ihlines], nl).reshape((nl, nhlines))
	y[:, ihlines] = by_rep_cols[:]

	# Infinities to nans
	y[np.where(np.abs(y) == np.inf)] = np.nan
	y = alg.make_sym_rem_nans(y)

	return x, y

def intersection_segments(a, b):
	""" 
	Intersection point between two segments a and b 
	"""

	x_int, y_int = intersection_two_lines(a.line, b.line)

	# Now check if this point belongs to the segments

	if a.contains_point((x_int, y_int), include_ends=False) \
		& b.contains_point((x_int, y_int), include_ends=False):
		return x_int, y_int
	else:
		return np.nan, np.nan

def crossing_segments(a, b):
	""" 
	Determine if two segments a and b are crossing (their intersection 
	point is not on any of their limits
	"""
	x_int, y_int = intersection_segments(a, b)
	if np.isnan(np.array([x_int, y_int])).any():
		return False, (x_int, y_int)

	a0x, a0y = a.a
	a1x, a1y = a.b
	b0x, b0y = b.a
	b1x, b1y = b.b
	x = np.array([a0x, a1x, b0x, b1x])
	y = np.array([a0y, a1y, b0y, b1y])
	_, dists = dist_cloud(x - x_int, y - y_int)
	if (dists <= 1.e-9).any():
		return False, (x_int, y_int)
	return True, (x_int, y_int)

def overlapping_segments(a, b):
	"""
	Determine if two segments overlap. For this: 
	1) check that their lines are the same
	2) check that a point of one of them (any) falls onto the other 
	segment (eg. a point in a contained by b or viceversa)
	"""

	# Get their value of equation(x=0.)
	y0a = a.line.equation(0)
	y0b = b.line.equation(0)

	tolerance = 1.e-9
	same_y0 = nearly_equal(y0a, y0b, tolerance=tolerance)
	same_slope = nearly_equal(a.slope, b.slope, tolerance=tolerance)
	same_lines = same_y0 and same_slope

	# If the lines are not the same, bye bye
	if not same_lines:
		return False
	
	# Check if any point from a falls on b or viceversa
	for p in a.points:
		if b.contains_point(p, include_ends=False):
			return True
	for p in b.points:
		if a.contains_point(p, include_ends=False):
			return True

	# If nothing happened so far, the segments don't overlap
	return False

def intersection_polygons(a, b):
	""" Return the polygon, if it exists, resulting from the 
	intersection between two polygons a and b """
	
	# List the points that are inside one another
	inter_points = []
	for pa in a.verts:
		if b.contains_point(pa):
			inter_points.append(pa)
	for pb in b.verts:
		if a.contains_point(pb):
			inter_points.append(pb)

	# Find the intersections between the segments

	for sa in a.segments.values():
		for sb in b.segments.values():
			xi, yi = intersection_segments(sa, sb)
			if (not np.isnan(xi)) and (not np.isnan(yi)):
				inter_points.append(np.array([xi, yi]))

	# Remove repeated points
	avoid_these_points = []
	for i, p1 in enumerate(inter_points):
		for j, p2 in enumerate(inter_points[i+1:]):
			j += i + 1
			# print(i, j)
			if pts_equal(p1, p2, 1.e-9):
				avoid_these_points.append(j)

	# Remove redundant points
	inter_points = np.array([p for (i, p) in enumerate(inter_points) if not i in avoid_these_points])

	# Finally, create the polygon/segment/point and return it if it has any points
	if len(inter_points) > 1:
		return Polygon(sort_by_angle(inter_points))
	else:
		return None

###############################################################################


class Point():
	"""
	Point.
	NOTE: The majority of the code here is written to treat points as just 
	pairs of values, not to use this class. But this class is defined here 
	to be used in some circumstances and also to use it in a future version 
	of the code.
	"""

	def __init__(self, x, y):
		self.x = x
		self.y = y
		self.array = np.array([x, y])

	def draw(self, ax, colour='k', marker='o', markersize=1, alpha=1.):
		ax.plot(
				self.x, 
				self.y, 
				color=colour, 
				marker=marker, 
				markersize=markersize, 
				alpha=alpha, 
			)

class StraightLine():
	"""
	Straight Line
	"""

	def __init__(self, point, slope):
		self.point = point
		self.slope = slope
		self.ps = (point, slope)
		self.set_equation()
		self.vertical = np.abs(self.slope) == np.inf
		self.horizontal = nearly_equal(np.abs(self.slope), 0., tolerance=1.e-9)

	def set_equation(self):
		"""
		Equation for the straight line
		"""
		p, s = self.point, self.slope
		self.y_ord = p[1] - s * p[0]
		y = self.y_ord
		self.equation = lambda x: y + s * x

	def draw(self, ax, colour='k', linewidth=1, alpha=1.):
		"""
		Draw the line on ax
		"""
		xleft, xright = ax.get_xlim()
		yleft, yright = self.equation(xleft), self.equation(xright)
		if (self.slope == np.inf) or (self.slope == -np.inf) or (np.isnan(self.slope)):
			xleft = xright = self.ps[0][0]
			yleft, yright = ax.get_ylim()
		ax.plot(
				[xleft, xright], 
				[yleft, yright], 
				color=colour, 
				lw=linewidth, 
				alpha=alpha, 
			)

	def contains_point(self, p):
		"""
		Check if a point p is on this line
		"""
		xp, yp = [tls.list_to_number(x) for x in p]
		if self.vertical:
			return nearly_equal(xp, self.point[0], tolerance=1.e-9)
		return nearly_equal(yp, self.equation(xp), tolerance=1.e-9)


class Segment():
	"""
	Segment
	"""

	def __init__(self, points):

		# Validation
		if len(points) != 2:
			msg = 'WARNING: Defined a segment with %i points:\n%s\nSegment not created. Returned None'%(len(points), str(points))
			try:
				ws.log(msg)
			except AttributeError:
				print(msg)
			return None

		# Clean the points. Give them a nice and clean format
		points = (
			tuple([tls.list_to_number(x) for x in points[0]]), 
			tuple([tls.list_to_number(x) for x in points[1]])
		)
		self.points_list = list(points)
		# Add attributes
		self.update(points)

	def __str__(self):
		""" Print the segment with its values """
		return 'Segment. a: (%f, %f); b: (%f, %f)'%tuple([x for x in self.points.flatten()])

	def update(self, points):
		""" Update the values of the segment """
		self.points = np.array(list(points))
		self.a = points[0]
		self.b = points[1]
		# Get the straight line onto which the segment is
		self.slope = get_slope(*self.points.tolist())
		self.line = StraightLine(self.a, self.slope)
		self.vertical = self.line.vertical
		self.horizontal = self.line.horizontal
		# Length
		self.length = dist(self.a, self.b)

	def change_point(self, index, point):
		""" Change one point of the segment; hence change it """
		# Clean the point
		point = tuple([tls.list_to_number(x) for x in point])
		# Add it to the list
		self.points_list[index] = point
		# Update the segment
		self.update(self.points_list)

	def contains_point(self, p, include_ends=True):
		""" Check if a point is on the segment """

		# Format point correctly
		p = tuple([tls.list_to_number(x) for x in p])

		# Now check if the line contains the point
		# If not, return False
		if not self.line.contains_point(p):
			return False

		# Now check that the point falls in the rectangle where seg is the 
		# diagonal
		xp, yp = p
		a, b = self.a, self.b

		# Note: if the segment is vertical, include_ends = True
		include_ends = include_ends or self.vertical or self.horizontal

		if include_ends:
			inx = (a[0] <= xp <= b[0]) | (a[0] >= xp >= b[0])
			iny = (a[1] <= yp <= b[1]) | (a[1] >= yp >= b[1])
		else:
			inx = (a[0] < xp < b[0]) | (a[0] > xp > b[0])
			iny = (a[1] < yp < b[1]) | (a[1] > yp > b[1])

		if inx and iny:
			return True
		else:
			return False

	def draw(self, ax, colour='k', linewidth=1, alpha=1.):
		""" Draw the segment on ax """
		ax.plot(
				(self.a[0], self.b[0]), 
				(self.a[1], self.b[1]), 
				c=colour, 
				lw=linewidth, 
				alpha=alpha, 
			)


class Circle():
	"""
	Circle
	"""

	def __init__(self, c, r):
		self.c = c
		self.r = r

	def draw(self, ax, colour='k', linewidth=1, alpha=1.):
		""" 
		Draw the circle 
		"""
		plt_circle = plt.Circle(self.c, self.r, color=colour, 
			linewidth=linewidth, alpha=alpha, fill=False)
		ax.add_artist(plt_circle)

		# If ax limits are smaller than the circle, expand them
		self.xmin = self.c[0] - self.r
		self.xmax = self.c[0] + self.r
		self.ymin = self.c[1] - self.r
		self.ymax = self.c[1] + self.r
		ax_xmin, ax_xmax = ax.get_xlim()
		ax_ymin, ax_ymax = ax.get_ylim()
		margin = 0.05 * max([self.xmax - self.xmin, self.ymax - self.ymin])
		# x-axis
		if self.xmin < ax_xmin:
			ax.set_xlim(left=self.xmin - margin)
		if self.xmax > ax_xmax:
			ax.set_xlim(right=self.xmax + margin)
		# y-axis
		if self.ymin < ax_ymin:
			ax.set_ylim(bottom=self.ymin - margin)
		if self.ymax > ax_ymax:
			ax.set_ylim(top=self.ymax + margin)


class Polygon(object):
	"""
	Class for any polygon that is defined with a list of verts
	*args and **kwargs account for any arguments and keyword arguments that may enter in instances of childs of this class
	"""

	def __new__(cls, verts, *args, **kwargs):

		lv = len(verts)

		# Validation
		if lv == 0:
			msg = 'ERROR: Point list or array is empty.'
			try:
				ws.log(msg)
			except AttributeError:
				print(msg)
			else:
				ws.terminate()

		elif lv == 1:
			msg = 'WARNING: Defined polygon with only 1 point:\n%s. Returing point'%str(verts)
			try:
				ws.log(msg)
			except AttributeError:
				print(msg)
			return Point(*verts[0])

		elif lv == 2:
			msg = 'WARNING: Defined polygon with only 2 points:\n%s. Returing Segment'%str(verts)
			try:
				ws.log(msg)
			except AttributeError:
				print(msg)
			return Segment(verts)

		# elif lv == 3:
			# This is a triangle
			# return super(Triangle, cls).__new__(cls)
			# return super(Triangle, cls).__new__(cls)
			# I didn't get this to work...
			pass

		elif lv >= 3:
			# Now, this is a proper polygon. So, proceed.
			return super(Polygon, cls).__new__(cls)
		
	def __init__(self, verts, *args):
		# print('2: Creating polygon')

		# Attributes
		self.verts = verts
		# If it's not an array, make it be
		try:
			verts.shape
		except AttributeError:
			self.verts = np.array(verts)
		self.nverts = len(verts)
		self.x, self.y = self.verts.T
		# Compute attributes
		self.set_segments()
		self.set_angles()
		self.get_area()
		# self.get_centroid()
		self.get_circumcircle()
		# planar.Polygon. I need it for planar.Polygon.contains_point
		self.plpol = pl.Polygon(self.verts.tolist())

	def set_segments(self):
		"""
		Set the segments of this polygon. Also their lengths and slopes
		"""

		segments = OrderedDict()
		sides = []
		slopes = []

		# Iterate over self.nverts - 2 points
		for i in range(self.nverts - 1):
			seg = Segment(self.verts[i:i+2])
			sides.append(seg.length)
			slopes.append(seg.slope)
			segments[(i, i+1)] = seg

		# Last segment
		seg = Segment((self.verts[i+1], self.verts[0]))
		sides.append(seg.length)
		slopes.append(seg.slope)
		segments[(i+1, 0)] = seg

		# Add as attribute
		self.segments = segments
		self.sides = np.array(sides)
		self.slopes = np.array(slopes)

	def set_angles(self):
		"""
		Get the angles for each one of the vertices
		Do it COUNTER-CLOCKWISE, which is how the points are oriented
		"""
		self.angles = np.zeros(self.nverts)
		for i in np.arange(0, self.nverts - 1, 1):
			c, b, a = self.verts[i - 1], self.verts[i], self.verts[i + 1]
			self.angles[i] = angle_threepts(a, b, c)
		# Last vertex
		c, b, a = self.verts[-2], self.verts[-1], self.verts[0]
		self.angles[-1] = angle_threepts(a, b, c)

		# Extra:
		# Tell if the polygon is convex
		self.is_convex = (self.angles <= np.pi).all()

		# Angles to degrees
		self.angles_deg = (180. / np.pi) * self.angles

	def get_area(self):
		"""
		Calculate the area of this polygon
		"""
		x, y = self.x, self.y
		self.signed_area = 0.5 * (np.dot(x, np.roll(y, 1)) \
			- np.dot(y, np.roll(x, 1)))
		# NOTE: This only works for polygons with no crossing segments, right?
		# WARNING: True, this is correct only for polygons with non-intersecting segments
		self.area = np.abs(self.signed_area)

	def get_centroid(self):
		"""
		Get the centroid of the polygon
		ERROR: SOMETHING IS WRONG WITH THIS
		"""
		oosixa = 1. / (6. * self.signed_area)
		x, y = self.x, self.y
		a = x[:-1] * y[1:] - x[1:] * y[:-1]
		cx = oosixa * ((x[:-1] + x[1:]) * a).sum()
		cy = oosixa * ((y[:-1] + y[1:]) * a).sum()
		self.centroid = (cx, cy)

	def get_circumcircle(self):
		""" 
		Create the circumcircle with its properties 
		"""

		# First and foremost, we have to try with the two most distant 
		# points. They will give us THE SMALLEST POSSIBLE CIRCLE
		# If this circle does not contain all the points (it does not
		# meet the condition), then we use the algorithm of the 
		# triangles (which always works)

		# Find the two most distant points
		xv, yv = self.verts.T
		pairs, distances = dist_cloud(xv, yv)
		self.most_distant_points = self.verts[pairs[np.where(distances == distances.max())][0]]
		# Get the centre
		a, b = self.most_distant_points
		circumcentre = get_midpoint(a, b)
		circumradius = dist(circumcentre, a)
		circumcircle = Circle(circumcentre, circumradius)

		# Provisionally save the previous circumcircle for comparison
		self.circumcircle_prov = circumcircle

		# Check if it meets the condition
			# self.circumcircle = circumcircle
		if not (points_in_circle(self.verts.T, circumcircle)).all():
		# else:
			# It didn't meet the condition, so go with the triangles

			# First off, get all the possible triangles
			# And get their circumcircles
			self.triangles = {}
			possible_circumcircles = []
			vertices = range(self.nverts)
			# Loop over all possible triangles
			for i in vertices:
				for j in vertices[i+1:]:
					for k in vertices[j+1:]:
						tverts = (
								self.verts[i], 
								self.verts[j], 
								self.verts[k]
							)

						# Create triangle
						# print('tverts:', tverts)
						triangle = Triangle(tverts)

						# Check if the circumcircle contains all points
						if (points_in_circle(
								self.verts.T, triangle.circumcircle)
							).all():
							pc = triangle.circumcircle

							if len(possible_circumcircles) == 0:
								# If this circumcircle is the smallest one so far, keep it
								circumcircle = pc
								circumcentre = pc.c
								circumradius = pc.r
							elif pc.r < circumcircle.r:
								# If we have more possibilities, choose it if it's the smallest one
								circumcircle = pc
								circumcentre = pc.c
								circumradius = pc.r

							# Even if this was positive, add it to the pc's 
							# that meet the condition
							possible_circumcircles.append(pc)
						# Save triangle
						self.triangles[(i, j, k)] = triangle
			# Save the possible circumcircles
			self.possible_circumcircles = possible_circumcircles

		# Save attributes
		self.circumcircle = circumcircle
		self.circumcentre = circumcentre
		self.circumradius = circumradius
		self.circumdiameter = 2. * circumradius

		# print('Polygon.')
		# print('Centroid:', self.centroid)
		# print('Circumcentre:', self.circumcentre)
		
	def contains_point(self, p, include_contour=True):
		"""
		Check if the polygon contains a point p
		"""

		# First, check if its planar.Polygon has it
		if self.plpol.contains_point(p):
			return True

		# If it doesn't, check if it's on the contour
		if include_contour:
			for seg in self.segments.values():
				if seg.contains_point(p, include_ends=True):
					return True

		# If we're here, the point is not in or on this polygon
		return False

	def draw(self, ax, colour='k', alpha=1., linestyle='-', linewidth=1., 
		markers=None, facecolor='none', dashes=(1, 1), zorder=0):
		"""
		Draw the polygon on ax
		"""
		# Add the last point
		ptstoplot = self.verts.tolist()
		ptstoplot.append(ptstoplot[0])
		ptstoplot = np.array(ptstoplot)
		if facecolor == 'none':
			try:
				if linestyle == '--':
					ax.plot(ptstoplot[:, 0], ptstoplot[:, 1], c=colour, alpha=alpha,
						ls=linestyle, dashes=dashes, lw=linewidth, marker=markers, zorder=zorder)
				else:
					ax.plot(ptstoplot[:, 0], ptstoplot[:, 1], c=colour, alpha=alpha,
						ls=linestyle, lw=linewidth, marker=markers, zorder=zorder)
			except IndexError:
				pass
		else:
			patches = []
			filled_pol = mp.Polygon(self.verts)
			patches.append(filled_pol)
			collection = PatchCollection(patches, facecolors=facecolor)
			ax.add_collection(collection)
			collection.set_alpha(alpha)
			collection.set_edgecolor(colour)
			collection.set_linewidth(linewidth)
			collection.set_linestyle(linestyle)


class Triangle(Polygon):
	"""
	Class for a Triangle that contains all the important parts
	I need
	"""
	# def __init__(self, verts, *args, indices=(0, 1, 2), fas=None):
	def __init__(self, verts, indices=(0, 1, 2), fas=None):
		# print('1: Creating triangle')
		self.indices = indices
		self.verts = np.array([verts[i] for i in indices])
		Polygon.__init__(self, self.verts)
		# Get relevant attributes
		# print('Here we have a triangle with %i vertices'%self.nverts)
		self.stablish_pairs()
		self.stablish_lines()
		# self.get_sides()
		self.get_incentre()
		self.get_medians()
		self.get_altitudes()
		self.get_circumcircle()
		self.get_orthodistances()
		# Check if it is inside a fascicle
		# Depending on that, the pseudo_incentre will be calculated or not
		if fas is not None:
			self.fas = fas
			self.fibers = fas.fibers
			self.get_pseudo_incentre()

	def stablish_pairs(self):
		self.pairs = []
		self.pairs.append((self.indices[1], self.indices[2]))
		self.pairs.append((self.indices[2], self.indices[0]))
		self.pairs.append((self.indices[0], self.indices[1]))

	def stablish_lines(self):
		# NOTE: Lines are actually 'segments', to be more precise
		# Lines
		# The index should be the same as the index of the opposite 
		# vertex
		self.lines = []
		self.lines.append([self.verts[1], self.verts[2]])
		self.lines.append([self.verts[2], self.verts[0]])
		self.lines.append([self.verts[0], self.verts[1]])

	# def get_sides(self):
	# 	# Line lengths or sides and slopes
	# 	sides = []
	# 	slopes = []
	# 	for line in self.lines:
	# 		sides.append(dist(line[0], line[1]))
	# 		# Slopes for the lines
	# 		slopes.append(get_slope(line[0], line[1]))
	# 	self.sides = np.array(sides)
	# 	self.slopes = np.array(slopes)

	def get_incentre(self):
		# Position of the incentre
		self.incentre = np.array([self.sides[i] \
		 * np.array(self.verts[i]) for i, p in \
		  enumerate(self.verts)]).sum(axis=0) / self.sides.sum()		
		# print(self.sides)
		a, b, c = self.sides
		s = 0.5 * self.sides.sum()
		self.inradius = np.sqrt((s - a) * (s - b) * (s - c) / s)		

	def get_pseudo_incentre(self):
		"""
		If the incentre is inside a cell, move it outwards
		"""
		# First, the pseudo_incentre is a copy of incentre
		self.pseudo_incentre = list(self.incentre)
		# Check if incentre is inside any of the triangle's cells
		for i in self.indices:
			cell = self.fibers[i]
			circ = cell.circumph
			c, r = circ
			# circ = sg.Circle(tup2sympoint(c), r)
			if inside_circle(self.incentre, circ):
				# Move pseudo_incentre
				# Place it on the cell's boundary, on the line 
				# perpendicular to the triangle's wall opposite to its
				# center. So first, find the equation of that line

				# Get the indices of the cells pairing on the side 
				# I want
				othercells = [j for j in self.indices if j != i]
				slope = get_slope(self.fibers[othercells[0]].center, 
					self.fibers[othercells[1]].center)
				orthoslope = get_orthoslope(slope)
				ortholine = StraightLine(c, orthoslope)
				# ortholine = sg.Line(cell.center, slope=orthoslope)
				inters = intersection_sl_circ(ortholine, circ)
				# Of the two points from inters, I need the closest 
				# one to the incentre
				self.pseudo_incentre = find_closest(inters, self.incentre)

	def get_medians(self):
		""" Calculate the medians of the triangle """
		a, b, c = self.sides
		self.medians = 0.5 * np.array([
			np.sqrt(2. * (b*b + c*c) - a*a),
			np.sqrt(2. * (a*a + c*c) - b*b), 
			np.sqrt(2. * (b*b + a*a) - c*c)
			])
		xa, xb, xc = self.verts

	def get_altitudes(self):
		""" Get altitudes from the sides """
		a, b, c = self.sides
		s = 0.5 * self.sides.sum()
		self.altitudes = \
			2. * np.sqrt(s * (s - a) * (s - b) * (s - c)) / self.sides

	def get_circumcircle(self):
		""" Get and save the circumcircle of this triangle """
		midpoints = []
		side_bisectors = []
		for slope, segment in zip(self.slopes, self.segments.values()):
			a, b = segment.a, segment.b
			mp = get_midpoint(a, b)
			midpoints.append(mp)
			side_bisectors.append(StraightLine(mp, get_orthoslope(slope)))
		
		# Find the circumcentre, which is the intersection between two
		# of these lines
		circumcentre = intersection_two_lines(
			side_bisectors[0], 
			side_bisectors[1]
			)
		
		# Circumradius
		circumradius = dist(circumcentre, self.verts[0])

		# Get circumcircle
		self.circumcircle = Circle(circumcentre, circumradius)

		# Save attributes
		self.circumcentre = circumcentre
		self.circumradius = circumradius
		self.midpoints = midpoints
		self.side_bisectors = side_bisectors

	def get_orthodistances(self):
		""" Get distances from the sides to the orthocenter """
		self.orthodistances = np.empty(3)
		for i in range(3):
			# sign = np.sign(np.cos(self.angles[i]))
			sign = -np.sign(np.cos(self.angles[i]))
			# vh = sign * np.sqrt(4 * self.circumradius ** 2 - self.sides[i] ** 2)
			vh = np.sqrt(4 * self.circumradius ** 2 - self.sides[i] ** 2)
			self.orthodistances[i] = sign * (self.altitudes[i] - vh)
