"""
Code modified from 13 November 2018. Work in progress
"""

import os
from neuron import h, gui
from collections import OrderedDict

import workspace as ws
import simcontrol
import electrodes
import anatomy
import analysis as als
import biophysics as bio
import visualisation as vis



# Prepare the necessary stuff for the simulation
simcontrol.prepare_simulation()

# Run the simulation
simcontrol.run()

# Postprocessing
als.postprocessing()