

## ^^^ Leave two lines blank at top which will be filled by CMake BASIS

##############################################################################
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2016 Imperial College London
# Copyright 2016 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#	 http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############################################################################

"""Command-line tool for reconstruction of neonatal cortex

This command implements the deformable surfaces method for the reconstruction of
the neonatal cortex as detailed in the conference paper submission to ISBI 2017.

See -help output for details.
"""

import os
import re
import sys
import csv
import argparse
import traceback

try:
	from contextlib import ExitStack  # Python 3
except:
	from contextlib2 import ExitStack  # Python 2 backport

try:
	from configparser import SafeConfigParser  # Python 3
except:
	from ConfigParser import SafeConfigParser  # Python 2

import mirtk.deformable.neonatal_cortex as neoctx


# ==============================================================================
# neonatal cortex reconstruction pipeline
# ==============================================================================

# ------------------------------------------------------------------------------
def get_default_config(work_dir='.', section='recon-neonatal-cortex'):
	"""Get default configuration."""
	# directories
	session_dir = os.path.join('%(WorkDir)s', '%(SubjectId)s-%(SessionId)s')
	input_dir   = os.path.join(session_dir, 'input')
	temp_dir	= os.path.join(session_dir, 'temp')
	mesh_dir	= os.path.join(session_dir, 'meshes')
	logs_dir	= os.path.join(session_dir, 'logs')
	# configuration
	config  = SafeConfigParser(defaults={'work_dir': work_dir, 'temp_dir': temp_dir})
	section = args.section
	config.add_section(section)
	config.set(section, 'logs_dir', logs_dir)
	# input file paths
	config.set(section, 'input_t1w_image',	 os.path.join(input_dir, 't1w-image.nii.gz'))
	config.set(section, 'input_t2w_image',	 os.path.join(input_dir, 't2w-image.nii.gz'))
	config.set(section, 'input_brain_mask',	os.path.join(input_dir, 'brain-mask.nii.gz'))
	config.set(section, 'input_labels_image',  os.path.join(input_dir, 'brain-labels.nii.gz'))
	config.set(section, 'input_tissues_image', os.path.join(input_dir, 'tissue-labels.nii.gz'))
	# intermediate file paths
	config.set(section, 't1w_image',			 os.path.join(temp_dir, 't1w-image.nii.gz'))
	config.set(section, 't2w_image',			 os.path.join(temp_dir, 't2w-image.nii.gz'))
	config.set(section, 'brain_mask',			os.path.join(temp_dir, 'brain-mask.nii.gz'))
	config.set(section, 'white_matter_mask',	 os.path.join(temp_dir, 'white-matter-mask.nii.gz'))
	config.set(section, 'gray_matter_mask',	  os.path.join(temp_dir, 'gray-matter-mask.nii.gz'))
	config.set(section, 'deep_gray_matter_mask', os.path.join(temp_dir, 'deep-gray-matter-mask.nii.gz'))
	config.set(section, 'corpus_callosum_mask',  os.path.join(temp_dir, 'corpus-callosum-mask.nii.gz'))
	config.set(section, 'ventricles_mask',	   os.path.join(temp_dir, 'ventricles-mask.nii.gz'))
	config.set(section, 'ventricles_dmap',	   os.path.join(temp_dir, 'ventricles-dmap.nii.gz'))
	config.set(section, 'regions_mask',		  os.path.join(temp_dir, 'regions.nii.gz'))
	config.set(section, 'cortical_hull_dmap',	os.path.join(temp_dir, 'cortical-hull-dmap.nii.gz'))
	# output file paths
	config.set(section, 'brain_mesh',			os.path.join(mesh_dir, 'brain.vtp'))
	config.set(section, 'bs_cb_mesh',			os.path.join(mesh_dir, 'bs+cb.vtp'))
	config.set(section, 'internal_mesh',		 os.path.join(mesh_dir, 'internal.vtp'))
	config.set(section, 'cerebrum_mesh',		 os.path.join(mesh_dir, 'cerebrum.vtp'))
	config.set(section, 'right_cerebrum_mesh',   os.path.join(mesh_dir, 'cerebrum-rh.vtp'))
	config.set(section, 'left_cerebrum_mesh',	os.path.join(mesh_dir, 'cerebrum-lh.vtp'))
	config.set(section, 'white_mesh',			os.path.join(mesh_dir, 'white.vtp'))
	config.set(section, 'right_white_mesh',	  os.path.join(mesh_dir, 'white-rh.vtp'))
	config.set(section, 'left_white_mesh',	   os.path.join(mesh_dir, 'white-lh.vtp'))
	config.set(section, 'pial_mesh',			 os.path.join(mesh_dir, 'pial.vtp'))
	config.set(section, 'right_pial_mesh',	   os.path.join(mesh_dir, 'pial-rh.vtp'))
	config.set(section, 'left_pial_mesh',		os.path.join(mesh_dir, 'pial-lh.vtp'))
	# parameters of subdivide-brain-image step to create regions mask
	config.set(section, 'subcortex_closing', '5')
	config.set(section, 'brainstem_closing', '5')
	config.set(section, 'cerebellum_closing', '5')
	# default labels used when none specified are those of Draw-EM (all_labels)
	config.set(section, 'white_matter_labels', '51..82')
	config.set(section, 'gray_matter_labels', '5..16,20..39')
	config.set(section, 'deep_gray_matter_labels', '1..4,40..47,85..87')
	config.set(section, 'lateral_ventricles_labels', '49,50')
	config.set(section, 'corpus_callosum_labels', '48')
	config.set(section, 'inter_hemisphere_labels', '40..47,85..87')
	config.set(section, 'brainstem_labels', '19')
	config.set(section, 'cerebellum_labels', '17,18')
	rh_labels = []  # right hemisphere structures
	rh_labels.extend(range(2,  48, 2))
	rh_labels.extend(range(50, 63, 2))
	rh_labels.extend(range(63, 82, 2))
	rh_labels.append(86)
	config.set(section, 'right_hemisphere_labels', ','.join([str(x) for x in rh_labels]))
	lh_labels = []  # left hemisphere structures
	lh_labels.extend(range(1,  18, 2))
	lh_labels.extend(range(21, 62, 2))
	lh_labels.extend(range(64, 83, 2))
	lh_labels.append(87)
	config.set(section, 'left_hemisphere_labels', ','.join([str(x) for x in lh_labels]))
	return config


# ------------------------------------------------------------------------------
def get_labels(config, section, option):
	"""Get labels list from configuration."""
	return config.get(section, option).split(',')


# ------------------------------------------------------------------------------
def require_input_brain_mask(config, section, config_vars, stack, verbose=0):
	"""Create brain mask from segmentation if none provided."""
	input_brain_mask = config.get(section, 'input_brain_mask', vars=config_vars)
	if os.path.isfile(input_brain_mask):
		return input_brain_mask
	if verbose > 0:
		print("Creating brain mask from segmentation image")
	segmentation = config.get(section, 'input_labels_image', vars=config_vars)
	if not segmentation or not os.path.isfile(segmentation):
		segmentation = config.get(section, 'input_tissues_image', vars=config_vars)
	if not segmentation or not os.path.isfile(segmentation):
		raise Exception("Input brain segmentation required")
	neoctx.binarize(
		name=input_brain_mask,
		segmentation=segmentation,
		temp=config.get(section, 'temp_dir', vars=config_vars)
	)
	return neoctx.push_output(stack, input_brain_mask)


# ------------------------------------------------------------------------------
def require_regions_mask(config, section, config_vars, stack, verbose=0):
	"""Create regions label image from segmentations if none provided."""
	regions_mask = config.get(section, 'regions_mask', vars=config_vars)
	if os.path.isfile(regions_mask):
		return regions_mask
	require_input_brain_mask(config, section, config_vars, stack)
	if verbose > 0:
		print("Dividing brain volume into disjoint regions")
	segmentation = config.get(section, 'input_labels_image', vars=config_vars)
	if not segmentation or not os.path.isfile(segmentation):
		raise Exception("Input brain segmentation required")
	tissues_image = config.get(section, 'input_tissues_image', vars=config_vars)
	if not os.path.isfile(tissues_image):
		tissues_image = None
	white_labels = get_labels(config, section, 'white_matter_labels')
	white_labels.extend(get_labels(config, section, 'deep_gray_matter_labels'))
	white_labels.extend(get_labels(config, section, 'corpus_callosum_labels'))
	white_labels.extend(get_labels(config, section, 'lateral_ventricles_labels'))
	subcortex_labels = get_labels(config, section, 'inter_hemisphere_labels')
	subcortex_labels.extend(get_labels(config, section, 'corpus_callosum_labels'))
	cortical_hull_dmap = config.get(section, 'cortical_hull_dmap', vars=config_vars)
	if cortical_hull_dmap == '':
		cortical_hull_dmap = None
	neoctx.subdivide_brain(
		name=regions_mask,
		segmentation=segmentation,
		tissues=tissues_image,
		white_labels=white_labels,
		cortex_labels=get_labels(config, section, 'gray_matter_labels'),
		right_labels=get_labels(config, section, 'right_hemisphere_labels'),
		left_labels=get_labels(config, section, 'left_hemisphere_labels'),
		subcortex_labels=subcortex_labels,
		subcortex_closing=config.getint(section, 'subcortex_closing'),
		brainstem_labels=get_labels(config, section, 'brainstem_labels'),
		brainstem_closing=config.getint(section, 'brainstem_closing'),
		cerebellum_labels=get_labels(config, section, 'cerebellum_labels'),
		cerebellum_closing=config.getint(section, 'cerebellum_closing'),
		brain_mask=config.get(section, 'input_brain_mask', vars=config_vars),
		cortical_hull_dmap=cortical_hull_dmap,
		temp=config.get(section, 'temp_dir', vars=config_vars),
		merge_bs_cb=True
	)
	if cortical_hull_dmap:
		neoctx.push_output(stack, cortical_hull_dmap)
	return neoctx.push_output(stack, regions_mask)


# ------------------------------------------------------------------------------
def require_brain_mask(config, section, config_vars, stack, verbose=0, keep_regions_mask=False):
	"""Resample input brain mask to standard RAS space."""
	brain_mask = config.get(section, 'brain_mask', vars=config_vars)
	if os.path.isfile(brain_mask):
		return brain_mask
	require_input_brain_mask(config, section, config_vars, stack, verbose)
	if keep_regions_mask:
		require_regions_mask(config, section, config_vars, None, verbose)
	else:
		require_regions_mask(config, section, config_vars, stack, verbose)
	if verbose > 0:
		print("Resampling brain mask to standard RAS space")
	neoctx.binarize(
		name=brain_mask,
		segmentation=config.get(section, 'input_brain_mask', vars=config_vars),
		image=config.get(section, 'regions_mask', vars=config_vars),
		threshold=0,
		temp=config.get(section, 'temp_dir', vars=config_vars)
	)
	if stack:
		brain_mask = neoctx.push_output(stack, brain_mask)
	neoctx.run('close-image', args=[brain_mask, brain_mask], opts={'iterations': 5})
	return brain_mask


# ------------------------------------------------------------------------------
def optional_corpus_callosum_mask(config, section, config_vars, stack, verbose=0):
	"""Create corpus callosum mask from segmentation if none provided."""
	corpus_callosum_mask = config.get(section, 'corpus_callosum_mask', vars=config_vars)
	if os.path.isfile(corpus_callosum_mask):
		return corpus_callosum_mask
	if verbose > 0:
		print("Creating corpus callosum mask from segmentation image")
	segmentation = config.get(section, 'input_labels_image', vars=config_vars)
	if segmentation and os.path.isfile(segmentation):
		require_brain_mask(config, section, config_vars, stack, verbose)
		neoctx.binarize(
			name=corpus_callosum_mask,
			segmentation=segmentation,
			labels=get_labels(config, section, 'corpus_callosum_labels'),
			image=config.get(section, 'brain_mask', vars=config_vars),
			temp=config.get(section, 'temp_dir', vars=config_vars)
		)
		if stack:
			return neoctx.push_output(stack, corpus_callosum_mask)
		else:
			return corpus_callosum_mask
	else:
		return None


# ------------------------------------------------------------------------------
def require_white_matter_mask(config, section, config_vars, stack, verbose=0):
	"""Create white matter mask from segmentation if none provided."""
	white_matter_mask = config.get(section, 'white_matter_mask', vars=config_vars)
	if os.path.isfile(white_matter_mask):
		return white_matter_mask
	if verbose > 0:
		print("Creating white matter mask from segmentation image")
	segmentation = config.get(section, 'input_tissues_image', vars=config_vars)
	if not os.path.isfile(segmentation):
		segmentation = config.get(section, 'input_labels_image', vars=config_vars)
		if not segmentation or not os.path.isfile(segmentation):
			raise Exception("Input brain segmentation required")
	neoctx.binarize(
		name=white_matter_mask,
		segmentation=segmentation,
		labels=get_labels(config, section, 'white_matter_labels'),
		image=config.get(section, 'brain_mask', vars=config_vars),
		temp=config.get(section, 'temp_dir', vars=config_vars)
	)
	if stack:
		return neoctx.push_output(stack, white_matter_mask)
	else:
		return white_matter_mask


# ------------------------------------------------------------------------------
def require_gray_matter_mask(config, section, config_vars, stack, verbose=0):
	"""Create gray matter mask from segmentation if none provided."""
	gray_matter_mask = config.get(section, 'gray_matter_mask', vars=config_vars)
	if os.path.isfile(gray_matter_mask):
		return gray_matter_mask
	if verbose > 0:
		print("Creating gray matter mask from segmentation image")
	segmentation = config.get(section, 'input_tissues_image', vars=config_vars)
	if not os.path.isfile(segmentation):
		segmentation = config.get(section, 'input_labels_image', vars=config_vars)
		if not segmentation or not os.path.isfile(segmentation):
			raise Exception("Input brain segmentation required")
	neoctx.binarize(
		name=gray_matter_mask,
		segmentation=segmentation,
		labels=get_labels(config, section, 'gray_matter_labels'),
		image=config.get(section, 'brain_mask', vars=config_vars),
		temp=config.get(section, 'temp_dir', vars=config_vars)
	)
	if stack:
		return neoctx.push_output(stack, gray_matter_mask)
	else:
		return gray_matter_mask


# ------------------------------------------------------------------------------
def require_deep_gray_matter_mask(config, section, config_vars, stack, verbose=0):
	"""Create deep gray matter mask from segmentation if none provided."""
	deep_gray_matter_mask = config.get(section, 'deep_gray_matter_mask', vars=config_vars)
	if os.path.isfile(deep_gray_matter_mask):
		return deep_gray_matter_mask
	if verbose > 0:
		print("Creating deep gray matter mask from segmentation image")
	segmentation = config.get(section, 'input_tissues_image', vars=config_vars)
	if not os.path.isfile(segmentation):
		segmentation = config.get(section, 'input_labels_image', vars=config_vars)
		if not segmentation or not os.path.isfile(segmentation):
			raise Exception("Input brain segmentation required")
	neoctx.binarize(
		name=deep_gray_matter_mask,
		segmentation=segmentation,
		labels=get_labels(config, section, 'deep_gray_matter_labels'),
		image=config.get(section, 'brain_mask', vars=config_vars),
		temp=config.get(section, 'temp_dir', vars=config_vars)
	)
	if stack:
		return neoctx.push_output(stack, deep_gray_matter_mask)
	else:
		return deep_gray_matter_mask


# ------------------------------------------------------------------------------
def require_ventricles_mask(config, section, config_vars, stack, verbose=0):
	"""Create lateral ventricles mask from segmentation if none provided."""
	ventricles_mask = config.get(section, 'ventricles_mask', vars=config_vars)
	if os.path.isfile(ventricles_mask):
		return ventricles_mask
	if verbose > 0:
		print("Creating lateral ventricles mask from segmentation image")
	segmentation = config.get(section, 'input_tissues_image', vars=config_vars)
	if not os.path.isfile(segmentation):
		segmentation = config.get(section, 'input_labels_image', vars=config_vars)
		if not segmentation or not os.path.isfile(segmentation):
			raise Exception("Input brain segmentation required")
	neoctx.binarize(
		name=ventricles_mask,
		segmentation=segmentation,
		labels=get_labels(config, section, 'lateral_ventricles_labels'),
		image=config.get(section, 'brain_mask', vars=config_vars),
		temp=config.get(section, 'temp_dir', vars=config_vars)
	)
	if stack:
		return neoctx.push_output(stack, ventricles_mask)
	else:
		return ventricles_mask


# ------------------------------------------------------------------------------
def require_ventricles_dmap(config, section, config_vars, stack, verbose=0):
	"""Compute lateral ventricles distance image if none provided."""
	ventricles_dmap = config.get(section, 'ventricles_dmap', vars=config_vars)
	if os.path.isfile(ventricles_dmap):
		return ventricles_dmap
	require_ventricles_mask(config, section, config_vars, stack)
	if verbose > 0:
		print("Computing lateral ventricles distance map")
	ventricles_mask = config.get(section, 'ventricles_mask', vars=config_vars)
	neoctx.calculate_distance_map(iname=ventricles_mask, oname=ventricles_dmap)
	if stack:
		return neoctx.push_output(stack, ventricles_dmap)
	else:
		return ventricles_dmap


# ------------------------------------------------------------------------------
def require_cortical_hull_dmap(config, section, config_vars, stack, verbose=0):
	"""Compute distance map of cortical hull."""
	require_regions_mask(config, section, config_vars, stack, verbose)
	return config.get(section, 'cortical_hull_dmap', vars=config_vars)


# ------------------------------------------------------------------------------
def recon_neonatal_cortex(config, section, config_vars,
						  with_brain_mesh=False,
						  with_bs_cb_mesh=True,
						  with_cerebrum_mesh=False,
						  with_white_mesh=True,
						  with_pial_mesh=True,
						  keep_t1w_image=False,
						  keep_t2w_image=False,
						  keep_regions_mask=False,
						  pial_outside_white_surface=False,
						  join_internal_mesh=False,
						  join_bs_cb_mesh=False,
						  cut=True,
						  force=False,
						  check=True,
						  verbose=0,
						  join_tol = 1,
						  use_fast_collision=False):
	"""Reconstruct surfaces of neonatal cortex."""

	# working directory
	temp_dir = config.get(section, 'temp_dir', vars=config_vars)

	# intermediate image files in standard RAS space
	t1w_image			 = config.get(section, 't1w_image',			 vars=config_vars)
	t2w_image			 = config.get(section, 't2w_image',			 vars=config_vars)
	brain_mask			= config.get(section, 'brain_mask',			vars=config_vars)
	wm_mask			   = config.get(section, 'white_matter_mask',	 vars=config_vars)
	gm_mask			   = config.get(section, 'gray_matter_mask',	  vars=config_vars)
	deep_gray_matter_mask = config.get(section, 'deep_gray_matter_mask', vars=config_vars)
	regions_mask		  = config.get(section, 'regions_mask',		  vars=config_vars)
	cortical_hull_dmap	= config.get(section, 'cortical_hull_dmap',	vars=config_vars)
	ventricles_dmap	   = config.get(section, 'ventricles_dmap',	   vars=config_vars)

	# output mesh files
	brain_mesh		  = config.get(section, 'brain_mesh',		  vars=config_vars)
	bs_cb_mesh		  = config.get(section, 'bs_cb_mesh',		  vars=config_vars)
	internal_mesh	   = config.get(section, 'internal_mesh',	   vars=config_vars)
	cerebrum_mesh	   = config.get(section, 'cerebrum_mesh',	   vars=config_vars)
	right_cerebrum_mesh = config.get(section, 'right_cerebrum_mesh', vars=config_vars)
	left_cerebrum_mesh  = config.get(section, 'left_cerebrum_mesh',  vars=config_vars)
	white_mesh		  = config.get(section, 'white_mesh',		  vars=config_vars)
	right_white_mesh	= config.get(section, 'right_white_mesh',	vars=config_vars)
	left_white_mesh	 = config.get(section, 'left_white_mesh',	 vars=config_vars)
	pial_mesh		   = config.get(section, 'pial_mesh',		   vars=config_vars)
	right_pial_mesh	 = config.get(section, 'right_pial_mesh',	 vars=config_vars)
	left_pial_mesh	  = config.get(section, 'left_pial_mesh',	  vars=config_vars)

	internal_base = os.path.splitext(os.path.basename(internal_mesh))[0]

	if not with_brain_mesh:
		brain_mesh = None
	if not with_bs_cb_mesh:
		bs_cb_mesh = None
	if not with_pial_mesh:
		pial_mesh = None

	if bs_cb_mesh and join_bs_cb_mesh:
		bs_cb_mesh_1 = bs_cb_mesh
		bs_cb_mesh_2 = None
	else:
		bs_cb_mesh_1 = None
		bs_cb_mesh_2 = bs_cb_mesh

	# remove intermediate files that did not exist before upon exit
	with ExitStack() as stack:

		recon_pial = with_pial_mesh and (force or not os.path.isfile(pial_mesh))
		recon_white = (with_white_mesh or recon_pial) and (force or not os.path.isfile(white_mesh))
		recon_cerebrum = (with_cerebrum_mesh or recon_white) and (force or not os.path.isfile(cerebrum_mesh))
		recon_brain = with_brain_mesh and (force or not os.path.isfile(brain_mesh))
		recon_bs_cb_mesh = ((with_bs_cb_mesh or (recon_cerebrum and bs_cb_mesh_1) or (recon_white and bs_cb_mesh_2)) and
							(force or not os.path.isfile(bs_cb_mesh)))

		# the surface reconstruction relies on a resampling of the intensity
		# images to the standard RAS space defined by the regions_mask / brain_mask
		if recon_brain or recon_bs_cb_mesh or recon_cerebrum or recon_white or recon_pial:
			require_brain_mask(config, section, config_vars, stack, verbose,
							   keep_regions_mask=keep_regions_mask)
		elif keep_regions_mask:
			# create regions mask only, skip surface reconstruction
			require_regions_mask(config, section, config_vars, None, verbose)

		if recon_white or recon_pial:
			if not os.path.isfile(t2w_image):
				input_t2w_image = config.get(section, 'input_t2w_image', vars=config_vars)
				if os.path.isfile(input_t2w_image):
					if verbose > 0:
						print("Resampling T2-weighted image to standard RAS space")
					neoctx.makedirs(t2w_image)
					neoctx.run(
						'transform-image',
						args=[
							input_t2w_image,
							t2w_image
						],
						opts={
							'interp': 'fast cubic bspline with padding',
							'Sp': 0,
							'dofin': 'Id',
							'target': brain_mask,
							'type': 'float'
						}
					)
					if keep_t2w_image:
						neoctx.push_output(stack, t2w_image)
				else:
					raise Exception("Input T2-weighted image required")

			if not os.path.isfile(t1w_image):
				input_t1w_image = config.get(section, 'input_t1w_image', vars=config_vars)
				if os.path.isfile(input_t1w_image):
					if verbose > 0:
						print("Resampling T1-weighted image to standard RAS space")
					neoctx.makedirs(t1w_image)
					neoctx.run(
						'transform-image',
						args=[
							input_t1w_image,
							t1w_image
						],
						opts={
							'interp': 'fast cubic bspline with padding',
							'Sp': 0,
							'dofin': 'Id',
							'target': brain_mask,
							'type': 'float'
						}
					)
					if not keep_t1w_image:
						neoctx.push_output(stack, t1w_image)
				else:
					if verbose > 0:
						print("No input T1-weighted image found, using only T2-weighted image")
					t1w_image = None

		# reconstruct boundary of brain mask
		if recon_brain:
			if verbose > 0:
				print("Reconstructing boundary of brain mask")
			neoctx.recon_brain_surface(name=brain_mesh, mask=brain_mask, temp=temp_dir)

		# reconstruct brainstem plus cerebellum surface
		if recon_bs_cb_mesh:
			if verbose > 0:
				print("Reconstructing brainstem plus cerebellum surface")
			neoctx.recon_brainstem_plus_cerebellum_surface(name=bs_cb_mesh, regions=regions_mask, temp=temp_dir)

		# reconstruct inner-cortical surface from segmentation
		if recon_cerebrum:

			# at the moment already ensured by require_brain_mask above...
			if keep_regions_mask:
				require_regions_mask(config, section, config_vars, None, verbose)
			else:
				require_regions_mask(config, section, config_vars, stack, verbose)

			# reconstruct inner-cortical surfaces of right and left hemispheres
			if force or not os.path.isfile(right_cerebrum_mesh):
				corpus_callosum_mask = optional_corpus_callosum_mask(config, section, config_vars, stack, verbose)
				if verbose > 0:
					print("Reconstructing boundary of right cerebral hemisphere segmentation")
				neoctx.recon_cortical_surface(name=right_cerebrum_mesh,
											  regions=regions_mask, hemisphere=neoctx.Hemisphere.Right,
											  corpus_callosum_mask=corpus_callosum_mask, temp=temp_dir, use_fast_collision=use_fast_collision)
			if force or not os.path.isfile(left_cerebrum_mesh):
				corpus_callosum_mask = optional_corpus_callosum_mask(config, section, config_vars, stack, verbose)
				if verbose > 0:
					print("Reconstructing boundary of left cerebral hemisphere segmentation")
				neoctx.recon_cortical_surface(name=left_cerebrum_mesh,
											  regions=regions_mask, hemisphere=neoctx.Hemisphere.Left,
											  corpus_callosum_mask=corpus_callosum_mask, temp=temp_dir, use_fast_collision=use_fast_collision)

			# join cortical surfaces of right and left hemispheres
			if verbose > 0:
				print("Joining surfaces of right and left cerebral hemispheres")
			neoctx.join_cortical_surfaces(name=cerebrum_mesh, regions=regions_mask,
										  right_mesh=right_cerebrum_mesh,
										  left_mesh=left_cerebrum_mesh,
										  bs_cb_mesh=bs_cb_mesh_1,
										  internal_mesh=internal_mesh,
										  temp=temp_dir, check=check, join_tol = join_tol)

			# remove cortical surfaces of right and left hemispheres
			if not with_cerebrum_mesh:
				os.remove(right_cerebrum_mesh)
				os.remove(left_cerebrum_mesh)

		# insert internal mesh into into initial inner-cortical surface
		if with_cerebrum_mesh and join_internal_mesh:
			cerebrum_prefix, cerebrum_ext = os.path.splitext(cerebrum_mesh)
			cerebrum_plus_internal_mesh = cerebrum_prefix + '+' + internal_base + cerebrum_ext
			if force or not os.path.isfile(cerebrum_plus_internal_mesh):
				if verbose > 0:
					print("Merging initial surface with internal mesh")
				neoctx.append_surfaces(cerebrum_plus_internal_mesh, surfaces=[cerebrum_mesh, internal_mesh], merge=True, tol=0)

		# reconstruct inner-cortical surface
		if recon_white:

			require_white_matter_mask(config, section, config_vars, stack, verbose)
			require_gray_matter_mask(config, section, config_vars, stack, verbose)
			require_deep_gray_matter_mask(config, section, config_vars, stack, verbose)
			require_ventricles_dmap(config, section, config_vars, stack, verbose)
			require_cortical_hull_dmap(config, section, config_vars, stack, verbose)

			if verbose > 0:
				print("Reconstructing inner-cortical surface")
			neoctx.recon_white_surface(name=white_mesh,
									   t1w_image=t1w_image, t2w_image=t2w_image,
									   wm_mask=wm_mask, gm_mask=gm_mask,
									   cortex_mesh=cerebrum_mesh, bs_cb_mesh=bs_cb_mesh_2,
									   subcortex_mask=deep_gray_matter_mask,
									   cortical_hull_dmap=cortical_hull_dmap,
									   ventricles_dmap=ventricles_dmap,
									   temp=temp_dir, check=check, use_fast_collision=use_fast_collision)

			# remove initial surface mesh
			if not with_cerebrum_mesh:
				os.remove(cerebrum_mesh)

		# insert internal mesh and cut surface at medial plane
		split_white = (cut and (force or not os.path.isfile(right_white_mesh) or not os.path.isfile(left_white_mesh)))
		white_prefix, white_ext = os.path.splitext(white_mesh)
		white_plus_internal_mesh = white_prefix + '+' + internal_base + white_ext
		if (split_white or join_internal_mesh) and (force or not os.path.isfile(white_plus_internal_mesh)):
			if verbose > 0:
				print("Merging inner-cortical surface with internal mesh")
			neoctx.append_surfaces(white_plus_internal_mesh, surfaces=[white_mesh, internal_mesh], merge=True, tol=0)
			if not join_internal_mesh:
				neoctx.push_output(stack, white_plus_internal_mesh)
		if split_white:
			if verbose > 0:
				print("Cutting inner-cortical surface at medial cutting plane")
			neoctx.split_cortical_surfaces(joined_mesh=white_plus_internal_mesh,
										   right_name=right_white_mesh,
										   left_name=left_white_mesh,
										   temp=temp_dir)

		# reconstruct outer-cortical surface
		if recon_pial:

			require_white_matter_mask(config, section, config_vars, stack, verbose)
			require_gray_matter_mask(config, section, config_vars, stack, verbose)

			if verbose > 0:
				print("Reconstructing outer-cortical surface")
			neoctx.recon_pial_surface(name=pial_mesh, t2w_image=t2w_image,
									  wm_mask=wm_mask, gm_mask=gm_mask, brain_mask=brain_mask,
									  white_mesh=white_mesh, bs_cb_mesh=bs_cb_mesh_2,
									  outside_white_mesh=pial_outside_white_surface,
									  temp=temp_dir, check=check, use_fast_collision=use_fast_collision)

			# remove inner-cortical surface
			if not with_white_mesh:
				os.remove(white_mesh)

		# insert internal mesh and cut surface at medial plane
		if with_pial_mesh:
			split_pial = (cut and (force or not os.path.isfile(right_pial_mesh) or not os.path.isfile(left_pial_mesh)))
			pial_prefix, pial_ext = os.path.splitext(pial_mesh)
			pial_plus_internal_mesh = pial_prefix + '+' + internal_base + pial_ext
			if (split_pial or join_internal_mesh) and (force or not os.path.isfile(pial_plus_internal_mesh)):
				if verbose > 0:
					print("Merging pial surface with internal mesh")
				neoctx.append_surfaces(pial_plus_internal_mesh, surfaces=[pial_mesh, internal_mesh], merge=True, tol=0)
				if not join_internal_mesh:
					neoctx.push_output(stack, pial_plus_internal_mesh)
			if split_pial:
				if verbose > 0:
					print("Cutting outer-cortical surface at medial cutting plane")
				neoctx.split_cortical_surfaces(joined_mesh=pial_plus_internal_mesh,
											   right_name=right_pial_mesh,
											   left_name=left_pial_mesh,
											   temp=temp_dir)


# ==============================================================================
# SLURM
# ==============================================================================


# ------------------------------------------------------------------------------
def sbatch(job_name, log_dir, session, args, config_vars):
	"""Submits SLURM jobs to run this script, one for each subject."""
	from subprocess import Popen, PIPE
	if not os.path.isdir(log_dir):
		os.makedirs(log_dir)
	outlog = os.path.join(log_dir, job_name + '-%j.out')
	errlog = os.path.join(log_dir, job_name + '-%j.err')
	p = Popen(
		[
			'sbatch', '--mem=4G', '-n', '1', '-c', str(args.threads),
			'-p', args.queue, '-o', outlog, '-e', errlog, '-J', job_name
		],
		stdout=PIPE, stderr=PIPE, stdin=PIPE
	)
	args_map = {
		'interpreter': sys.executable,
		'script': __file__,
		'config': args.config,
		'section': args.section,
		'work_dir': args.work_dir,
		'session': session,
		'threads': args.threads,
		'verbose': ' '.join(['-v'] * args.verbose),
		'debug': ' '.join(['-d'] * args.debug)
	}
	script = "#!/bin/sh\nexec {interpreter} {script} --threads={threads} {verbose} {debug}"
	script += " --work-dir='{work_dir}' --config='{config}' --section='{section}' --session='{session}'"
	if args.brain:
		script += ' --brain'
	if args.cerebrum:
		script += ' --cerebrum'
	if args.white:
		script += ' --white'
	if args.pial:
		script += ' --pial'
	if args.force:
		script += ' --force'
	if not args.cut:
		script += ' --nocut'
	if not args.check:
		script += ' --nocheck'
	if args.join_internal_mesh:
		script += ' --join-with-internal-mesh'
	if args.join_bs_cb_mesh:
		script += ' --join-with-brainstem-and-cerebellum'
	if args.pial_outside_white:
		script += ' --ensure-pial-is-outside-white-surface'
	if args.keep_t1w_image:
		script += ' --keep-t1w-image'
	if args.keep_t2w_image:
		script += ' --keep-t2w-image'
	if args.keep_regions_mask:
		script += ' --keep-regions-mask'
	for name, value in config_vars.items():
		if "'" in value:
			value = value.replace("'", "\\'")
		name = name.replace('-', '_')
		script += ' --' + name + "='" + value + "'"
	(out, err) = p.communicate(input=script.format(**args_map).encode('utf-8'))
	if p.returncode != 0:
		raise Exception(err)
	m = re.match('Submitted batch job ([0-9]+)', out)
	if m:
		return int(m.group(1))
	return out


# ==============================================================================
# main
# ==============================================================================


def split_config_args(args):
	"""Split -[-]name=value configuration command-line arguments."""
	res = []
	for arg in args:
		if arg.startswith('-'):
			for part in arg.split('=', 1):
				res.append(part)
		else:
			res.append(arg)
	return res


# parse arguments
parser = argparse.ArgumentParser(description='Reconstruct neonatal cortex from MR brain scan and Draw-EM segmentation')
parser.add_argument('-r', '-root', '--root', '-work-dir', '--work-dir', dest='work_dir', default=os.getcwd(),
					help='Root working directory')
parser.add_argument('-c', '-config', '--config', default='',
					help='Optional custom configuration file')
parser.add_argument('-section', '--section', default='recon-neonatal-cortex',
					help='Configuration section name')
parser.add_argument('-s', '-sessions', '--sessions', default=[], nargs='+',
					help="Either list of '{SubjectId}[-{SessionId}]' strings or path of CSV file", required=True)
parser.add_argument('-b', '-brain', '--brain', action='store_true',
					help='Reconstruct surface of brain mask')
parser.add_argument('-w', '-white', '--white', action='store_true',
					help='Reconstruct white surface')
parser.add_argument('-cerebrum', '--cerebrum', action='store_true',
					help='Reconstruct/keep initial white surface')
parser.add_argument('-p', '-pial', '--pial', action='store_true',
					help='Reconstruct pial surface')
parser.add_argument('-ensure-pial-is-outside-white-surface', '--ensure-pial-is-outside-white-surface',
					dest='pial_outside_white', action='store_true',
					help='Ensure that pial surface is strictly outside the white surface')
parser.add_argument('-join-with-internal-mesh', '--join-with-internal-mesh',
					dest='join_internal_mesh', action='store_true',
					help='Join final mesh with internal (hemispheres) dividier mesh')
parser.add_argument('-join-with-brainstem-and-cerebellum', '--join-with-brainstem-and-cerebellum',
					dest='join_bs_cb_mesh', action='store_true',
					help="Merge cerebrum surface mesh with brainstem and cerebellum surface mesh")
parser.add_argument('-nocut', '-nosplit', '--nocut', '--nosplit', dest='cut', action='store_false',
					help='Save individual (closed) genus-0 surfaces for each hemisphere')
parser.add_argument('-use-fast-collision', dest='fastcollision', action='store_true',
					help='Use the fast collision test')
parser.add_argument('-nocheck', '--nocheck', action='store_false', dest='check',
					help='Disable consistency and self-intersection checks of (intermediate) surface meshes')
parser.add_argument('-keep-t1w-image', '--keep-t1w-image', action='store_true',
					help="Keep resampled T1-weighted image even when no -debug option given")
parser.add_argument('-keep-t2w-image', '--keep-t2w-image', action='store_true',
					help="Keep resampled T2-weighted image even when no -debug option given")
parser.add_argument('-keep-regions-mask', '--keep-regions-mask', action='store_true',
					help="Keep regions label image even when no -debug option given")
parser.add_argument('-f', '-force', '--force', action='store_true',
					help='Overwrite existing output files')
parser.add_argument('-v', '-verbose', '--verbose', action='count', default=0,
					help='Increase verbosity of output messages')
parser.add_argument('-d', '-debug', '--debug', action='count', default=0,
					help='Keep/write debug output in temp_dir')
parser.add_argument('-t', '-threads', '--threads', default=0,
					help='No. of cores to use for multi-threading')
parser.add_argument('-j', '-jointol', '--jointol', dest='join_tol', default=1, type = float,
					help='Join tolerance')
parser.add_argument('-q', '-queue', '--queue', default='',
					help='SLURM partition/queue')

[args, config_args] = parser.parse_known_args()
#print config_args
#print args
args.work_dir = os.path.abspath(args.work_dir)
if not args.cerebrum and not args.white and not args.pial and not args.keep_regions_mask:
	args.white = True
	args.pial = True
elif args.pial:
	args.white = True

config_vars = {}
config_args = split_config_args(config_args)
if len(config_args) % 2 != 0:
	raise Exception("Custom configuration options must come in pairs of -[-]<name> <value>:\n{}".format(config_args))
for i in range(0, len(config_args), 2):
	name = config_args[i]
	if name.startswith('--'):
		name = name[2:]
	elif name.startswith('-'):
		name = name[1:]
	else:
		raise Exception("Custom configuration options must start with either one or two dashes")
	name = name.replace('-', '_')
	config_vars[name] = config_args[i + 1]

# read configuration
config = get_default_config(work_dir=args.work_dir, section=args.section)
config.read(os.path.join(args.work_dir, 'recon-neonatal-cortex.cfg'))
if args.config:
	with open(args.config, 'r') as config_file:
		config.readfp(config_file)

# set global flags
neoctx.verbose = max(0, args.verbose - 1)
neoctx.showcmd = max(0, args.verbose - 1)
neoctx.debug = max(0, args.debug)
neoctx.force = args.force

# read subject and session IDs from CSV file
if len(args.sessions) == 1 and os.path.isfile(args.sessions[0]):
	csv_name = args.sessions[0]
	sessions = []
	with open(csv_name) as f:
		reader = csv.DictReader(f)
		for row in reader:
			if 'SubjectID' in row:
				session = row['SubjectID']
			elif 'SubjectId' in row:
				session = row['SubjectId']
			else:
				raise Exception("Missing 'SubjectID' or 'SubjectId' column in CSV file")
			if 'SessionID' in row:
				session += '-' + row['SessionID']
			elif 'SessionId' in row:
				session += '-' + row['SessionId']
			sessions.append(session)
else:
	sessions = args.sessions

# for each session...
failed = 0
for session in sessions:
	match = re.match('^(.*)-([^-]+)$', session)
	if match:
		subject_id = str(match.group(1))
		session_id = match.group(2)
	else:
		subject_id = str(session)
		session_id = '0'
		session = session + '-0'
	info = {
		'subid': subject_id,
		'subject_id': subject_id,
		'subjectid': subject_id,
		'SubjectID': subject_id,
		'SubjectId': subject_id,
		'sesid': session_id,
		'session_id': session_id,
		'sessionid': session_id,
		'SessionID': session_id,
		'SessionId': session_id,
		'WorkDir': args.work_dir,
		'work_dir': args.work_dir
	}
	try:
		if args.queue:
			sys.stdout.write("Submitting SLURM job for {SubjectId} session {SessionId}: ".format(**info))
			job_name = 'rec-{SubjectId}-{SessionId}'.format(**info)
			log_dir  = config.get(args.section, 'logs_dir', vars=info)
			job_id   = sbatch(job_name, log_dir, session, args, config_vars)
			sys.stdout.write('Job ID = {}\n'.format(job_id))
		else:
			sys.stdout.write("\nReconstructing cortical surfaces of {SubjectId} session {SessionId}\n".format(**info))
			config_vars.update(info)
			#print list(config.options('recon-neonatal-cortex'))
			#print config_vars

			recon_neonatal_cortex(config=config, section=args.section, config_vars=config_vars,
								  with_brain_mesh=args.brain,
								  with_cerebrum_mesh=args.cerebrum,
								  with_white_mesh=args.white,
								  with_pial_mesh=args.pial,
								  keep_t1w_image=args.keep_t1w_image,
								  keep_t2w_image=args.keep_t2w_image,
								  keep_regions_mask=args.keep_regions_mask,
								  pial_outside_white_surface=args.pial_outside_white,
								  join_internal_mesh=args.join_internal_mesh,
								  join_bs_cb_mesh=args.join_bs_cb_mesh,
								  verbose=args.verbose,
								  check=args.check,
								  join_tol = args.join_tol,
								  use_fast_collision = args.fastcollision)
	except Exception as e:
		failed += 1
		if args.queue:
			sys.stdout.write("failed\n")
		sys.stdout.write("\n")
		if args.verbose > 0 or args.debug > 0:
			exc_type, exc_value, exc_traceback = sys.exc_info()
			traceback.print_exception(exc_type, exc_value, exc_traceback)
		else:
			sys.stderr.write('Exception: {}\n'.format(str(e)))
if failed > 0:
	sys.exit(1)
