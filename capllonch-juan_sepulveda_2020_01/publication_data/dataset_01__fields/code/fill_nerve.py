"""
This program fills the nerve with all the components of the resistive 
network (axons, NAELC...), including the tesselation
"""

import os
import csv
import matplotlib.pyplot as plt

import anatomy
import geometry as geo
import analysis as als
import workspace as ws
import simcontrol



###############################################################################
# Data management

# Prepare the folders that the code will use
simcontrol.get_paths()

# Remove results from previous simulation(s)
als.remove_previous_results()

###############################################################################

# Prepare the necessary variables from the configuration files
simcontrol.prepare_workspace()

###############################################################################
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

# Number of axons, locations (centers) and radii
params['numberofaxons'] = apsettings["number"]
params['locations'] = apsettings["locations"]
params['radii'] = apsettings["radii"]

# Minimum and maximum radii
params['rmin'] = packsettings["min radius"]
params['rmax'] = packsettings["max radius"]
# Tolerance for trying to pack circles
params['circp_tol'] = packsettings["max iterations"]
# Maximum number of axons per fascicle
params['max_axs_pf'] = packsettings["max axons per fascicle"]

# Build tessellated nerve
nvt = anatomy.NerveTess()
nvt.build(params)

# Save it
spath = os.path.join(ws.folders["data/saved"], 'created_nerve_internal_topology.csv')
nvt.save_to_file(spath)

# Create figure to visualise the tessellation process
fig, ax = plt.subplots()
ws.ax = ax
nvt.pd.draw(ax)
nvt.draw_axons_and_points(ax)
nvt.draw_contours(ax, c='r')
ax.set_aspect('equal')
plt.show()