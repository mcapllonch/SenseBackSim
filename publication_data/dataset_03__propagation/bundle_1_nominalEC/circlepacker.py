"""
Circle packing algorithm(s)
"""
import numpy as np
import planar as pl


def fill(contour, rmin=0., rmax=1.e99, min_sep=0., nmaxpf=None, 
	tolerance=1e4, tol_ta=1, distribution="uniform", distr_params=None):

	# In case it's not an array:
	try:
		_ = contour.size
	except:
		contour = np.array(contour)

	# Polygon
	pts = [tuple(item) for item in contour.tolist()]
	polygon = pl.Polygon(pts)
	center = (contour[:-1, 0].mean(), contour[:-1, 1].mean())
	diameter = 2 * np.hypot(contour[:-1, 0] - center[0], 
		contour[:-1, 1] - center[1]).max()
	radius = 0.5 * diameter

	# List of fibers
	
	rlist = []
	positions = []

	tryanother = 0
	nonstop = True
	# Cell counter
	cc = 0

	# minsep = min_sep * 0.5

	while nonstop:
		# Random radius
		if distribution == "uniform":
			r = np.random.uniform(rmin, rmax, 1)
		elif distribution == "gamma":
			valid = False
			while not valid:
				# Work with diameters first, then back to radius
				mean = 2. * distr_params["mean"]
				shape = distr_params["shape"]
				scale = mean / shape
				r = 0.5 * np.random.gamma(shape, scale, 1)
				valid = (r >= rmin) & (r <= rmax)
 
		allowed_r = r + min_sep
		# Go trying different positions
		k = 0
		while (k < tolerance):
			# Random position inside the polygon
			correct = False
			while not correct:
				# Try a random position
				dx = np.random.uniform(-radius, -radius + diameter)
				dy = np.random.uniform(-radius, -radius + diameter)
				x = center[0] + dx
				y = center[1] + dy
				# 1: Check if the center is inside the polygon
				correct = polygon.contains_point((x, y))
				# 2: Check that the whole fiber is inside the polygon
				# For this, I actually check that all the points of the 
				# polygon are farther than the fiber radius
				cx, cy = contour[:, 0], contour[:, 1]
				correct = correct & (np.hypot(cx - x, cy - y) > allowed_r).all()
				k += 1

			# Check for clashes with other existing sub-units
			clash = False
			for i, xy in enumerate(positions):
				xc, yc = xy
				if (np.hypot(x - xc, y - yc) < allowed_r + rlist[i]):
					# Clash. This sub-unit doesn't fit here
					clash = True
					break

			if clash:
				k += 1
				# If I tried too many times with a sub-unit and it 
				# didn't fit, stop
				if k >= tolerance:
					tryanother += 1
					if tryanother >= tol_ta:
						nonstop = False
				# Else, keep trying
			else:
				# No clash at all. Good position
				xy = (x, y)
				positions.append(xy)
				rlist.append(r)
				cc += 1
				k = tolerance
				if cc >= nmaxpf:
					nonstop = False

	positions, rlist = np.array(positions), np.array(rlist)
	r = rlist[..., 0]
	x, y = positions.T
	return x, y, r