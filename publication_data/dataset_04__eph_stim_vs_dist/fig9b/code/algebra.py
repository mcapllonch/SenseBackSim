"""
Algebraic operations
"""

import csv
import numpy as np
from collections import OrderedDict
from scipy.stats import linregress

import algebra as alg
import geometry as geo
import workspace as ws


def parameter_regressions(filename):
	""" Open and read a file with parameters over which to make 
	regressions and obtain new expressions and values """
	# Read data
	with open(filename, 'r') as f:
		frl = list(csv.reader(f))
		keys = frl[0]
		data = np.array(frl[1:], dtype=np.float64)

	n_fields = len(keys)
	n_entries = len(data)
	data_dict = OrderedDict()
	for i in range(n_fields):
		data_dict[keys[i]] = data[:, i]

	# Independent variable
	indepvar = keys[0]

	# Dictionary containing the regressions
	regressions = OrderedDict()
	# Of course, indicate which one is the independent variable
	regressions['independent variable'] = indepvar

	for i, (k, v) in enumerate(list(data_dict.items())[1:]):

		x = data_dict[indepvar]

		# Linear regression
		slope, intercept, r_value, p_value, std_err = linregress(x, v)
		# Check if it's valid
		if intercept < 0:
			# Negative values, unacceptable. Use a polynomial fit
			# Quadratic fit
			parabolic_coeffs = np.polyfit(x, v, 2)
			if parabolic_coeffs[-1] < 0:
				# Not even a quadratic fit worked
				msg = 'ERROR: Could not obtain positive values for ' + k + ' for small ' + indepvar
				ws.log(msg)
				ws.terminate()
			# Give it a x-array that covers the whole domain
			x_ = np.linspace(0, x.max(), 100)
			ypol = polynomial(x_, parabolic_coeffs)
			# return polynomial regression
			a, b, c = parabolic_coeffs
			regressions[k] = lambda x: a * x ** 2 + b * x + c
		else:
			# Plot linear regression
			# Create straight line object and draw it
			sl = geo.StraightLine((0, intercept), slope)
			regressions[k] = sl.equation
	return regressions

def polynomial(x, c):
	""" Get the values of y from a polynomial of degree len(c) - 1 with 
	coefficients c, using abscissa x 
	The coefficients c must start from the higher order, just like in 
	np.polyfit (for consistency) """
	y = np.zeros(len(x), dtype=np.float64)
	for i, cc in enumerate(c[::-1]):
		y += cc * x ** i
	return y

def nan_to_inf(m):
	""" Turn all nan elements of a matrix m to infinity """
	m[np.where(np.isnan(m))] = np.inf
	return m

def diagv(v):
	""" Put the vector v on the diagonal of a square matrix """
	eye = np.eye(v.size)
	r = v * eye
	r[np.where(eye == 0.)] = 0.
	return r
	
def vectomat(v):
	""" Create a matrix from a vector by repeating it """
	n = v.size
	return v.repeat(n).reshape((n, n))

def compl_mat(v):
	""" Obtain a particular matrix from a vector """
	vm = vectomat(v)
	return vm - vm.T

def one_side(m):
	""" Return a vector containing one side of the off-diagonal terms of a 
	matrix and their indices """

	# Create list of pairs
	nc = m.shape[0]
	ids = np.arange(nc)
	idvm = vectomat(ids)
	# pairs = np.array([idvm, idvm.T]).reshape((nc, nc, 2))
	idvm2 = idvm.T

	# Create mask
	x = np.arange(nc)
	y = np.arange(nc)
	x, y = np.meshgrid(x, y)
	mask = x - y

	# Mask the values of m and pairs
	# pairs = np.ma.masked_where(mask <= 0, pairs)
	ids1 = np.ma.masked_where(mask <= 0, idvm)
	ids2 = np.ma.masked_where(mask <= 0, idvm2)
	m = np.ma.masked_where(mask <= 0, m)
	ids1 = ids1.compressed()
	ids2 = ids2.compressed()
	pairs = np.array([ids1, ids2]).T
	m = m.compressed()

	return pairs, m

def make_sym_rem_nans(m):
	""" Given a matrix m which is supposed to be symmetric but it's 
	not because there are nans, remove the non-due nans to make 
	the matrix symmetric """
	mmasked = np.ma.masked_where(np.isnan(m), m)
	mask = mmasked.mask
	nonan = np.where(mask.T == False)
	m[nonan] = mmasked.T[nonan]
	return m