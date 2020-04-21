"""
This is just a module containing handy tools for lists, dictionaries, etc.
"""

import os
import numpy as np

import workspace as ws



def flattenlist(l):
	""" Flatten a list of lists """
	return [item for sublist in l for item in sublist]

def list_to_number(l):
	""" Given a list or a list of lists with ONLY ONE NUMBER INSIDE, 
	retrieve that number.
	If there's more than one number, yield an error and quit """


	# Check if it's a list or a number
	try:
		ll = len(l)
	except TypeError:
		# It's probably a number, so just return it as it is
		return float(l)

	l = np.array(l).flatten()
	if len(l) != 1:
		msg = 'ERROR: tools.list_to_number got a list with more than one number'
		ws.log(msg)
		ws.terminate()

	return float(l[0])

def append_items(d, k, items):
	""" Append a list of items to a list which is a value of a dictionary 
	entry and create the entry if it does not exist yet """
	try:
		d[k].append(items[0])
	except KeyError:
		# Create entry
		d[k] = items[:]
	else:
		d[k] = d[k] + items[1:]
	return d

def arrtonum(x):
	""" Meant only for arrays with one size 1. 
	Return the element as a float """
	try:
		a = x.mean()
	except AttributeError:
		a = x
	return a

def str2bool(s):
	""" Very simple tool to parse the strings 'True' or 'False' into
	their corresponding booleans """
	str2booldict = {
		"True": True, 
		"False": False
	}
	try:
		return str2booldict[s]
	except KeyError:
		return s

def explore_folder(path):
	""" List all the items of a folder and return its sub-folders """
	items = []
	dirs = []
	for item in os.listdir(path):
		items.append(item)
		if os.path.isdir(item):
			dirs.append(item)
	return items, dirs