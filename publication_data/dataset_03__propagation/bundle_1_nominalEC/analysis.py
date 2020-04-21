"""
Miguel Capllonch Juan
This module performs analysis of the results
"""

import numpy as np
import os
import csv
from neuron import h


import workspace as ws



def remove_previous_results():
	""" Remove all previous results, either text or figures """
	resultsdir = ws.folders["data/results"]
	# Remove all .csv in the results
	for item in os.listdir(resultsdir):
		if item.endswith(".csv"):
			os.remove(os.path.join(resultsdir, item))
	# Remove all images
	imgdir = os.path.join(resultsdir, 'images')
	for item in os.listdir(imgdir):
		os.remove(os.path.join(imgdir, item))

def set_time_recording():
	""" Set the recording for the time vector of the simulation """

	# Record time
	ws.tvec = h.Vector()
	ws.tvec.record(h._ref_t)

def set_recording(seg, var):
	""" Set up a recording vector in hoc for any given segment 
	in the system """
	recording = h.Vector()
	# recording.buffer_size(ws.nt)
	exec('recording.record(seg._ref_%s)'%var)
	ws.recordings.append(recording)
	ws.counters['recordings'] += 1

def record_cables():
	""" Set up the recording vectors for the desired cables and the
	desired variables """
	dr = ws.settings["data records"]
	rc = dr["record cables"]
	ra = dr["record axons"]
	if rc == "all":
		# Record all types of cables
		for cable in ws.cables:
			if cable.cabletype == "NAELC":
				cable.set_recordings()
			elif cable.cabletype == "Axon":
				if ra == "all":
					# Record all axons
					cable.set_recordings()
	elif rc == "NAELC":
		# Record only NAELC
		for cable in ws.cables:
			if cable.cabletype == "NAELC":
				cable.set_recordings()
	elif rc == "axons":
		# Record only axons
		if ra == "all":
			# Record all axons
			for cable in ws.cables:
				if cable.cabletype == "Axon":
					cable.set_recordings()
	elif rc == "None":
		# Nothing to record
		pass

def remove_artefacts(time, recs):
	""" Remove the stimulation artefacts in the electrode recordings """

	# Start the analysis

	# Turn time and recordings into arrays to be able to work with them
	time = np.array(time)
	recs = np.array(recs)
	# Make everything positive
	rpos = np.abs(recs)
	# Average
	ravg = rpos.mean()
	# Identify where rpos is higher than its mean: there should be the artefacts
	rhigh = np.where(rpos > ravg)

	# Identify the separation between pulses
	dt = np.gradient(time).mean()
	# The 1.0000001 factor is just to make sure we get > and not >=
	try:
		cuts = np.where(np.gradient(time[rhigh]) > 1.000001 * dt)[0][0::2]
	except ValueError:
		# This may happen, for instance, if np.gradient can't be 
		# performed. In case of error, just say there's no cuts
		cuts = []
	# Split the pulses
	segments = []
	for cut in cuts:
		ttrr = []
		ttrr.append(time[rhigh][:cut + 1])
		ttrr.append(recs[rhigh][:cut + 1])
		segments.append(ttrr)

	# Last pulse
	ttrr = []
	try:
		ttrr.append(time[rhigh][cut + 1:])
		ttrr.append(recs[rhigh][cut + 1:])
	except UnboundLocalError:
		# There are no cuts
		pass
	else:
		segments.append(ttrr)

	# Identify the jumps at the cuts

	# Choose which data are outisde the pulses. Use rhigh as a mask
	toutsp = np.ma.masked_where(rpos > ravg, time)
	routsp = np.ma.masked_where(rpos > ravg, recs)
	# Get the jumps from that
	jumps = np.where(routsp.mask[1:] != routsp.mask[:-1])[0]
	# Pulse counter
	k = 0
	# (Average) jump heights of each pulse
	pulseheights = []
	# New values for the recordings during the pulses
	rnews = []
	# Iterate over jumps to get the avg. jump height for each pulse
	for i, jump in enumerate(jumps):
		jheight = abs(recs[jump + 1] - recs[jump])
		if i % 2 == 0:
			jumpup = jheight
		elif len(segments) > 0:
			rawrec = segments[k][1]
			avgjh = 0.5 * (jumpup + jheight) * np.sign(rawrec)
			# Subtract this value from the corresponding segment in recs
			rnew = rawrec - avgjh
			rnews.append(rnew)
			# Add one to the pulse counter
			k += 1
		else:
			rnews = recs.copy()

	# Finally, get the entire thing
	finalrecs = recs.copy()
	if len(segments) > 0:
		finalrecs[rhigh] = [item for sublist in rnews for item in sublist]
	return finalrecs

def save_vmem(data):
	""" Save the membrane potentials """
	f_sti = open(os.path.join(ws.folders["data/results"], 'vmemb_stm_site.csv'), 'w')
	f_rec = open(os.path.join(ws.folders["data/results"], 'vmemb_rec_site.csv'), 'w')
	writers = {
		'Sti': csv.writer(f_sti),
		'Rec': csv.writer(f_rec)
	}
	for i, vecs in enumerate(data):
		for k, vv in vecs.items():
			writers[k].writerow([i] + vv.as_numpy().tolist())
	f_sti.close()
	f_rec.close()

def save_electrode_recordings():
	""" Save the recordings from all electrodes into files """

	# Time
	tvarr = ws.tvec.as_numpy()

	# Iterate over electrodes:
	for name, electrode in ws.electrodes.items():
		if electrode.role == 'recording':
			# Iterate over rings on the electrode
			for ring in electrode.rings:
				ring_number = ring.ring_number

				# Iterate over pads on the ring
				for pad, rec in ring.recordings.items():
					# Recording ID
					rec_id = rec['recording id']
					# Fetch the vector with the data from ws
					v = ws.recordings[rec_id]

					# Save with artefacts
					with open(os.path.join(ws.folders["data/results"], 
						'recordings_E_%s_R%iP%i_withartefacts.csv'\
						%(name.replace(' ', '_'), ring_number, pad)), 'w') as f:
						fw = csv.writer(f)
						fw.writerow(['Time (ms)', 'V (mV)'])
						for tt, vv in zip(tvarr, v):
							fw.writerow([tt, vv])

					# Remove the artefacts
					try:
						v = remove_artefacts(tvarr, v.as_numpy())
					except Exception as e:
						message = "ERROR in %s.save_electrode_recordings:\n"%__name__ + type(e).__name__ + ':\n' + str(e)
						ws.log(message)
						# Don't terminate, just go on
						# ws.terminate()

					# Save
					with open(os.path.join(ws.folders["data/results"], 
						'recordings_E_%s_R%iP%i_noartefacts.csv'\
						%(name.replace(' ', '_'), ring_number, pad)), 'w') as f:
						fw = csv.writer(f)
						fw.writerow(['Time (ms)', 'V (mV)'])
						for tt, vv in zip(tvarr, v):
							fw.writerow([tt, vv])

def save_cable_recordings():
	""" Save the recordings from cables into files """

	# Time
	tvarr = ws.tvec.as_numpy()

	# Iterate over cables
	for cable in ws.cables:

		if cable.has_recordings:

			# Open file
			with open(os.path.join(ws.folders["data/results"], '%s%04i.csv'\
				%(cable.cabletype, cable.specific_number)), 'w') as f:
				fw = csv.writer(f)

				# Write the most basic properties of this cable:
				# Position and radius
				fw.writerow([cable.x, cable.y, cable.r])

				# Write the anatomy of the axon along the z-axis
				lp = ["segment positions:"]
				for i, seg in enumerate(ws.segments[cable.cable_number]):
					lp.append(str(seg))
					lp.append(cable.properties['zprofile'][i])
				fw.writerow(lp)

				# Iterate over recorded variables
				for var, recs in cable.recordings.items():

					# Write the name of the variable
					fw.writerow([var])
					
					# Iterate over recording vectors (or segments)
					for rec in recs:

						# Get numerical data
						data = ws.recordings[rec["number"]].as_numpy()

						# Store results
						# fw.writerow([ws.sectypes[cable.cable_number][i]] + [x for x in data])
						fw.writerow([rec["segment"]] + [x for x in data])

def save_memcurrents(i_mem_, i_mem_ly2):
	""" Save the membrane currents into files """
	x, y = ws.nvt.x, ws.nvt.y
	with open(os.path.join(ws.folders["data/results"], 'membranecurrents.csv'), 'w') as f:
		fw = csv.writer(f)
		fw.writerow(['Fiber', 'secname', 'x', 'y', 'z', 'i_membrane time series', 'i_membrane on layer 2'])
		for i, fiber in i_mem_.items():
			for j, segment in enumerate(fiber):
				fw.writerow([i, ws.segments[i][j].sec.name(), x[i], 
					y[i], ws.zprofs[i][j]] + segment.tolist() + \
					i_mem_ly2[i][j].tolist())

def store_i_mem(i_t):
	""" Store the membrane currents of all the cells.
	So far, this only works for axons with 2 extracellular layers """
	for i in np.arange(ws.nNAELC, ws.nc, 1):
		for j, seg in enumerate(ws.segments[i]):
			if 'node' in seg.sec.name():
				ws.i_mem_[i][j, i_t] = seg.i_membrane * seg.area()
				# Outwards is positive
				ws.i_mem_ly2[i][j, i_t] = (seg.vext[0] - seg.vext[1]) * seg.xg[0] * seg.area()
			else:
				ws.i_mem_[i][j, i_t] = 0.
				ws.i_mem_[i][j, i_t] = seg.i_membrane * seg.area()
				# Outwards is positive
				ws.i_mem_ly2[i][j, i_t] = (seg.vext[0] - seg.vext[1]) * seg.xg[0] * seg.area()

def postprocessing():
	""" Do all the post-processing tasks after a simulation """

	# Save the recordings from cables into files
	save_cable_recordings()

	# Save the electrode recordings into files
	# if that's what we want
	if ws.settings["electrode recordings"]:
		# and if the nerve model we're using is compatible with that
		if ws.settings["nerve model"] != "simple":
			save_electrode_recordings()