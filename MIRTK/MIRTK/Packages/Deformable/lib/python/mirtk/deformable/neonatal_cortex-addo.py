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

"""Python module for reconstruction of neonatal cortex

This module implements the deformable surfaces method for the reconstruction of
the neonatal cortex as detailed in the conference paper submission to ISBI 2017.
The recon-neonatal-cortex script provides a command-line tool for executing
the functions of this module to obtain explicit surface representations
of the inner and outer cortical surfaces given a Draw-EM segmentation
and the bias corrected T2-weighted (and T1-weighted if available) image(s).

Note: The merge-surfaces tool failed when MIRTK was built with VTK 6.2.
	  Use a more recent version of VTK if you encounter such issue.

"""

# Tested with MIRTK master branch commit f2f9bbb (Nov 14, 2016)
#
# Previous MIRTK versions (may not work with current script):
# - schuhschuh/MIRTK: develop branch commit 91e77be (Oct 10, 2016)
# - schuhschuh/MIRTK: develop branch commit 517e1e9 (Oct 19, 2016)
# - schuhschuh/MIRTK: develop branch commit e85a3eb (Oct 28, 2016)
# - schuhschuh/MIRTK: develop branch commit 115b40e (Nov  4, 2016)

import os
import re
import sys

from contextlib import contextmanager

try:
	from contextlib import ExitStack  # Python 3
except:
	from contextlib2 import ExitStack  # Python 2 backport

from mirtk.subprocess import flatten, check_call, check_output
from mirtk.subprocess import run as _run

# ==============================================================================
# global settings
# ==============================================================================

verbose = 0	# verbosity of output messages
showcmd = 0	# whether to print binary path and arguments of subprocesses
threads = 0	# maximum number of allowed threads for subprocess execution
			   # 0: use all available cores
debug   = 0	# debug level, keep intermediate files when >0
force   = True # whether to overwrite existing output files

_cortex_mask_array	= 'CortexMask'
_region_id_array	  = 'RegionId'
_collision_mask_array = 'CollisionMask'
_collision_type_array = 'CollisionType'

# ==============================================================================
# enumerations
# ==============================================================================

# ------------------------------------------------------------------------------
class Hemisphere(object):
	"""Enumeration of brain hemispheres"""
	Unspecified = 0
	Right	   = 1
	Left		= 2
	Both		= 3

# ------------------------------------------------------------------------------
def hemi2str(hemi):
	"""Convert Hemisphere enumeration value to string"""
	if   hemi == Hemisphere.Right: return 'rh'
	elif hemi == Hemisphere.Left:  return 'lh'
	else:						  return ''

# ==============================================================================
# path utilities
# ==============================================================================

# ------------------------------------------------------------------------------
def splitext(name):
	"""Split file name into name and extension."""
	(base, ext) = os.path.splitext(name)
	if ext.lower() == '.gz':
		(base, fext) = os.path.splitext(base)
		ext = fext + ext
	return (base, ext)

# ------------------------------------------------------------------------------
def splitname(name):
	"""Split temporary output file name into base name, incremental output ID, and extension."""
	(base, ext) = splitext(name)
	m = re.match('^(.*)-([0-9]+)$', base)
	if m:
		base = m.group(1)
		num  = int(m.group(2))
	else:
		num = 0
	return (base, num, ext)

# ------------------------------------------------------------------------------
def joinname(base, num, ext):
	"""Join temporary file name parts."""
	return '{}-{}{}'.format(base, num, ext)

# ------------------------------------------------------------------------------
def nextname(name, temp=None):
	"""Get next incremental output file path."""
	(base, num, ext) = splitname(name)
	if temp:
		base = os.path.join(temp, os.path.basename(base))
	return joinname(base, num+1, ext)

# ------------------------------------------------------------------------------
def makedirs(name):
	"""Make directories for output file if not existent."""
	path = os.path.dirname(name)
	if not os.path.isdir(path):
		os.makedirs(path)

# ------------------------------------------------------------------------------
def rename(src, dst):
	"""Rename file."""
	makedirs(dst)
	os.rename(src, dst)

# ------------------------------------------------------------------------------
def try_remove(name):
	"""Try to remove file, do not throw exception if it fails."""
	if name:
		try:
			os.remove(name)
		except:
			return False
	return True

# ==============================================================================
# pipeline utilities
# ==============================================================================

# ------------------------------------------------------------------------------
def iterable_or_int_to_list(arg):
	"""Given either an int or an iterable, returns list of items."""
	if isinstance(arg, int):
		return [arg]
	else:
		return list(arg)

# ------------------------------------------------------------------------------
def run(tool, args=[], opts={}):
	"""Run MIRTK command with global `showcmd` flag and maximum allowed number of `threads`."""
	_run(tool, args=args, opts=opts, showcmd=showcmd, threads=threads)

# ------------------------------------------------------------------------------
@contextmanager
def output(name_or_func, delete=False):
	"""Context with (temporary) output file path.

	This context is used to ensure that a temporary file is deleted on
	exit or when an exception occurs. Moreover, it is used to guarantee that an
	output file which is the result of a number of intermediate processes which
	write to the same file corresponds to the final output. When an error occurs
	before the final output was written to the specified file, the output file
	is removed when this context is exited to ensure that every output file has
	a known final state. Alternatively, use different unique output file name
	for each intermediate file and remove these intermediate files on exit.

	Parameters
	----------
	name_or_func : str, def
		File path or function which creates the output file and returns its path.
	delete : bool
		Whether the output file is temporary and should be deleted when done.

	Yields
	------
	aname : str
		Absolute path of output file.

	"""
	if isinstance(name_or_func, str): path = name_or_func
	elif isinstance(name_or_func, unicode): path = str(name_or_func) 
	else:							 path = name_or_func()
	if path:
		try:
			yield os.path.abspath(path)
		except BaseException as e:
			if debug < 1:
				try_remove(path)
			raise e
		if debug == 0 and path and delete:
			try_remove(path)
	else:
		raise Exception("Invalid file path")

# ------------------------------------------------------------------------------
def push_output(stack, name_or_func, delete=True):
	"""Perform interim processing step with auto-file removal upon exit."""
	return stack.enter_context(output(name_or_func, delete=delete))

# ------------------------------------------------------------------------------
def get_voxel_size(image):
	"""Get voxel size of image file."""
	info  = check_output(['info', image])
	match = re.search('Spacing:\s+([0-9][0-9.]*)\s*x\s*([0-9][0-9.]*)\s*x\s*([0-9][0-9.]*)', info)
	try:
		dx = float(match.group(1))
		dy = float(match.group(2))
		dz = float(match.group(3))
	except:
		raise Exception("Failed to determine voxel size of image: {}".format(image))
	return (dx, dy, dz)

# ------------------------------------------------------------------------------
def get_surface_property(name_or_info, props, mesh=False, topology=False, opts={}):
	"""Get surface property values"""
	values = []
	if os.path.splitext(name_or_info)[1] in ['.vtp', '.vtk', '.stl', '.obj']:
		info = evaluate_surface(name_or_info, mesh=mesh, topology=topology, opts=opts)
	else:
		info = name_or_info
	if isinstance(props, str):
		props = [props]
	for prop in props:
		regex = re.compile('^\s*' + re.escape(prop) + r'\s*=\s*(.+)\s*$', re.MULTILINE)
		match = regex.search(info)
		if not match:
			if debug > 0:
				sys.stderr.write(info)
			raise Exception("Missing surface property: " + prop)
		values.append(match.group(1))
	if len(values) == 1:
		return values[0]
	return values

# ------------------------------------------------------------------------------
def get_num_components(name_or_info):
	"""Get number of connected surface components"""
	return int(get_surface_property(name_or_info, 'No. of components', mesh=True))

# ------------------------------------------------------------------------------
def get_num_boundaries(name_or_info):
	"""Get number of surface boundaries"""
	return int(get_surface_property(name_or_info, 'No. of boundary segments', mesh=True))

# ------------------------------------------------------------------------------
def get_euler_characteristic(name_or_info):
	"""Get Euler characteristic of surface mesh"""
	return int(get_surface_property(name_or_info, 'Euler characteristic', topology=True))

# ------------------------------------------------------------------------------
def get_genus(name_or_info):
	"""Get Genus of surface mesh"""
	return float(get_surface_property(name_or_info, 'Genus', topology=True))

# ==============================================================================
# subroutines
# ==============================================================================

# ------------------------------------------------------------------------------
def invert_mask(iname, oname=None):
	"""Invert a binary mask."""
	if oname:
		makedirs(oname)
	else:
		oname = nextname(iname)
	run('calculate-element-wise', args=[iname], opts=[('binarize', [0, 0]), ('out', oname, 'binary')])
	return oname

# ------------------------------------------------------------------------------
def close_image(iname, oname=None, iterations=1, connectivity=18):
	"""Close (binary) image using morphological dilation followed by erosion."""
	if oname:
		makedirs(oname)
	else:
		oname = nextname(iname)
	run('close-image', args=[iname, oname], opts={'iterations': iterations, 'connectivity': connectivity})
	return oname

# ------------------------------------------------------------------------------
def del_mesh_attr(iname, oname=None, **opts):
	"""Delete surface mesh attributes."""
	if not oname:
		oname = iname
	run('delete-pointset-attributes', args=[iname, oname], opts=opts)

# ------------------------------------------------------------------------------
def evaluate_surface(name, oname=None, mesh=False, topology=False, intersections=False, collisions=0, opts={}):
	"""Evaluate properties of surface mesh"""
	argv = ['evaluate-surface-mesh', name]
	if oname:
		argv.extend([oname, '-v'])
	argv.extend(['-threads', str(threads)])
	if mesh:	 argv.append('-attr')
	if topology: argv.append('-topology')
	if collisions > 0 or intersections:
		argv.extend(['-collisions', collisions])
	if len(opts) > 0:
		if isinstance(opts, list):
			for item in opts:
				if isinstance(item, (tuple, list)):
					opt = item[0]
					arg = flatten(item[1:])
				else:
					opt = item
					arg = None
				if not opt.startswith('-'):
					opt = '-' + opt
				argv.append(opt)
				if not arg is None:
					argv.extend(flatten(arg))
		else:
			for opt, arg in opts.items():
				if not opt.startswith('-'):
					opt = '-' + opt
				argv.append(opt)
				if not arg is None:
					argv.extend(flatten(arg))
	return check_output(argv, verbose=showcmd)

# ------------------------------------------------------------------------------
def smooth_surface(iname, oname=None, iterations=1, lambda_value=1, mu=None, mask=None, weighting='combinatorial', excl_node=False):
	if not oname: oname = nextname(iname)
	argv = ['smooth-surface', iname, oname, '-threads', str(threads), '-iterations', iterations, '-' + weighting]
	if mask: argv.extend(['-mask', mask])
	if excl_node: argv.append('-exclnode')
	else:		 argv.append('-inclnode')
	argv.extend(['-lambda', lambda_value])
	if mu: argv.extend(['-mu', mu])
	check_call(argv, verbose=showcmd)
	return oname

# ------------------------------------------------------------------------------
def check_intersections(iname, oname=None):
	"""Check surface mesh for triangle/triangle intersections."""
	argv = ['evaluate-surface-mesh', iname]
	if oname:
		argv.extend([oname, '-v'])
	argv.extend(['-threads', threads, '-collisions', 0])
	info  = check_output(argv, verbose=showcmd)
	match = re.search('No\. of self-intersections\s*=\s*(\d+)', info)
	return int(match.group(1))

# ------------------------------------------------------------------------------
def remove_intersections(iname, oname=None, mask=None, max_attempt=10, smooth_iter=5, smooth_lambda=1):
	"""Remove intersections of surface mesh triangles.

	.. seealso:: MRISremoveIntersections function of FreeSurfer (dev/utils/mrisurf.c)
	"""
	if not oname:
		oname = nextname(iname)
	with output(oname):
		itr = nbr = 1
		cur = check_intersections(iname, oname)
		while cur > 0:
			itr += 1
			if itr == 1 and verbose > 0:
				print(("Trying to resolve {} remaining self-intersections of {}".format(cur, iname)))
			if itr > max_attempt:
				raise Exception("Failed to resolve self-intersections of {}".format(iname))
			run('dilate-scalars', args=[oname, oname], opts={'array': _collision_mask_array, 'iterations': nbr})
			if mask:
				run('calculate-element-wise', args=[oname], opts=[('scalars', _collision_mask_array), ('mul', mask), ('out', oname)])
			smooth_surface(oname, oname, iterations=smooth_iter, lambda_value=smooth_lambda,
						   weighting='combinatorial', mask=_collision_mask_array, excl_node=False)
			pre = cur
			cur = check_intersections(oname, oname)
			if cur >= pre: nbr += 1
		del_mesh_attr(oname, oname, pointdata=_collision_mask_array, celldata=_collision_type_array)
	return oname

# ------------------------------------------------------------------------------
def remove_white_pial_intersections(iname, oname=None, mask=None, max_attempt=10, smooth_iter=1, smooth_lambda=1, region_id_array=_region_id_array):
	"""Smooth white and pial surfaces in an attempt to remove intersections between triangles of white and pial surface."""
	iname = os.path.abspath(iname)
	if oname:
		oname = os.path.abspath(oname)
	else:
		oname = nextname(iname)
	with output(oname):
		# constants
		weighting		 = 'combinatorial'
		smooth_mask_array = 'SmoothMask'
		curvature_array   = 'Curvature'
		# check if triangles intersect
		cur = check_intersections(iname, oname)
		if cur > 0:
			# compute surface curvature
			run('calculate-surface-attributes', args=[oname, oname], opts={'H': curvature_array, 'vtk-curvatures': None})
			run('copy-pointset-attributes', args=[oname, oname], opts={'celldata-as-pointdata': region_id_array, 'unanimous': None})
		itr = nbr = 1
		while cur > 0:
			itr += 1
			if itr == 1 and verbose > 0:
				print(("Trying to resolve {} remaining self-intersections of {}".format(cur, iname)))
			if itr > max_attempt:
				raise Exception("Failed to resolve self-intersections of {}".format(iname))
			# smooth intersected convex regions (sulci) of pial surface
			run('calculate-element-wise', args=[oname], opts=[('point-data', region_id_array), ('label', 3, 4),
															  ('mask', _collision_mask_array), ('set', 1), ('pad', 0),
															  ('mul', curvature_array), ('threshold-ge', 0), ('pad', 0),
															  ('out', oname, 'binary', smooth_mask_array)])
			run('dilate-scalars', args=[oname, oname], opts={'array': smooth_mask_array, 'iterations': nbr})
			if mask:
				run('calculate-element-wise', args=[oname], opts=[('scalars', smooth_mask_array), ('mul', mask), ('out', oname)])
			smooth_surface(oname, oname, iterations=smooth_iter, lambda_value=smooth_lambda,
						   weighting=weighting, mask=smooth_mask_array, excl_node=False)
			if check_intersections(oname, oname) > 0:
				# smooth intersected concave regions (gyri) of white surface
				run('calculate-element-wise', args=[oname], opts=[('point-data', region_id_array), ('label', 1, 2),
																  ('mask', _collision_mask_array), ('set', 1), ('pad', 0),
																  ('mul', curvature_array), ('threshold-le', 0), ('pad', 0),
																  ('out', oname, 'binary', smooth_mask_array)])
				run('dilate-scalars', args=[oname, oname], opts={'array': smooth_mask_array, 'iterations': nbr})
				if mask:
					run('calculate-element-wise', args=[oname], opts=[('scalars', smooth_mask_array), ('mul', mask), ('out', oname)])
				smooth_surface(oname, oname, iterations=smooth_iter, lambda_value=smooth_lambda,
							   weighting=weighting, mask=smooth_mask_array, excl_node=False)
			# re-evaluate intersections
			pre = cur
			cur = check_intersections(oname, oname)
			if cur >= pre: nbr += 1
		# remove auxiliary attributes
		del_mesh_attr(oname, pointdata=[region_id_array, curvature_array, smooth_mask_array, _collision_mask_array], celldata=_collision_type_array)
	return oname

# ------------------------------------------------------------------------------
def project_mask(iname, oname, mask, name, dilation=0, invert=False, fill=False):
	"""Project binary mask onto surface."""
	makedirs(oname)
	run('project-onto-surface', args=[iname, oname],
		opts={'labels': mask, 'fill': fill, 'dilation-radius': dilation, 'name': name})
	if invert:
		run('calculate-element-wise', args=[oname],
			opts=[('scalars', name), ('threshold-inside', [0, 0]), ('set', 0), ('pad', 1), ('out', oname)])

# ------------------------------------------------------------------------------
def calculate_distance_map(iname, oname=None, temp=None, distance='euclidean', isotropic=1):
	"""Calculate signed distance map from binary object mask."""
	if not oname:
		if not temp:
			temp = os.path.dirname(iname)
		oname = splitname(os.path.basename(iname))[0]
		if oname.endswith('-mask'):
			oname = oname[0:-5]
		oname = os.path.join(temp, '{}-dmap.nii.gz'.format(oname))
	makedirs(oname)
	run('calculate-distance-map', args=[iname, oname], opts={'distance': distance, 'isotropic': isotropic})
	return oname

# ------------------------------------------------------------------------------
def extract_isosurface(iname, oname=None, temp=None, isoval=0, blur=0, close=True):
	"""Extract iso-surface from image."""
	if not oname:
		if not temp:
			temp = os.path.dirname(iname)
		oname = splitname(os.path.basename(iname))[0]
		if oname.endswith('-dmap') or oname.endswith('-mask'):
			oname = oname[0:-5]
		oname = '{}-iso.vtp'.format(oname)
		oname = os.path.join(temp, oname)
	makedirs(oname)
	run('extract-surface', args=[iname, oname], opts={'isovalue': isoval, 'blur': blur, 'close': close})
	return oname

# ------------------------------------------------------------------------------
def remesh_surface(iname, oname=None, temp=None, edge_length=0, min_edge_length=0, max_edge_length=10):
	"""Remesh surface."""
	if not oname:
		if not temp:
			temp = os.path.dirname(iname)
		oname = os.path.join(temp, nextname(os.path.basename(iname)))
	makedirs(oname)
	if edge_length > 0:
		min_edge_length = edge_length
		max_edge_length = edge_length
	run('remesh-surface', args=[iname, oname], opts={'min-edge-length': min_edge_length, 'max-edge-length': max_edge_length})
	return oname

# ------------------------------------------------------------------------------
def get_convex_hull(iname, oname=None, temp=None, remeshing=10, edge_length=1, smoothing=100):
	"""Get convex hull of surface mesh."""
	if not oname:
		oname = splitname(iname)[0]
		if oname.endswith('-iso'):
			oname = oname[0:-4]
		oname = '{}-hull.vtp'.format(oname)
	if force:
		try_remove(oname)
	if not os.path.isfile(oname):
		with ExitStack() as stack:
			# compute convex hull of isosurface
			hull = push_output(stack, nextname(oname, temp=temp))
			run('extract-pointset-surface', opts={'input': iname, 'output': hull, 'hull': None})
			# iteratively remesh and smooth surface to remove small self-intersections
			# which may have been introduced by the hull tesselation / local remeshing operations
			for i in range(0, remeshing):
				name = push_output(stack, nextname(hull))
				run('remesh-surface', args=[hull, name], opts={'edge-length': edge_length})
				run('smooth-surface', args=[name, name],
					opts={'combinatorial': None, 'exclnode': None,
						  'iterations': smoothing, 'lambda': .33, 'mu': -.34})
				if debug < 2:
					try_remove(hull)
				hull = name
			rename(hull, oname)
	return oname

# ------------------------------------------------------------------------------
def extract_convex_hull(iname, oname=None, temp=None, isoval=0, blur=0):
	"""Calculate convex hull of implicit surface."""
	if not temp:
		temp = os.path.dirname(iname)
	base = splitname(os.path.basename(iname))[0]
	if base.endswith('-dmap') or base.endswith('-mask'):
		base = base[0:-5]
	if not oname:
		oname = os.path.join(temp, '{}-hull.vtp'.format(base))
	if force:
		try_remove(oname)
	if not os.path.isfile(oname):
		makedirs(oname)
		with output(extract_isosurface(iname, temp=temp, isoval=isoval, blur=blur), delete=True) as iso:
			get_convex_hull(iso, oname, temp=temp)
	return oname

# ------------------------------------------------------------------------------
def add_corpus_callosum_mask(iname, mask, oname=None):
	"""Add ImplicitSurfaceFillMask point data array to surface file.
	
	This mask is considered by the ImplicitSurfaceDistance force when deforming
	the convex hull towards the WM segmentation boundary. Distance values of points
	with a zero mask value are excluded from the clustering based hole filling.
	This is to allow the surface to deform into the sulcus next to corpus callosum.
	"""
	if not oname:
		oname = nextname(iname)
	if force:
		try_remove(oname)
	if not os.path.isfile(oname):
		makedirs(oname)
		with output(oname):
			project_mask(iname, oname, mask, name='ImplicitSurfaceFillMask', dilation=20, invert=True)
	return oname

# ------------------------------------------------------------------------------
def del_corpus_callosum_mask(iname, oname=None):
	"""Remove ImplicitSurfaceFillMask again."""
	if not oname:
		oname = nextname(iname)
	if force:
		try_remove(oname)
	if not os.path.isfile(oname):
		makedirs(oname)
		del_mesh_attr(iname, oname, pointdata='ImplicitSurfaceFillMask')
	return oname

# ------------------------------------------------------------------------------
# TODO: Option to require minimum distance from medial cutting plane.
def add_cortex_mask(iname, mask, name=_cortex_mask_array, region_id_array=_region_id_array, oname=None):
	"""Add a CortexMask cell data array to the surface file."""
	if not oname:
		oname = nextname(iname)
	if debug > 0:
		assert os.path.realpath(iname) != os.path.realpath(oname), "iname != oname"
	if force or not os.path.isfile(oname):
		makedirs(oname)
		with output(oname):
			run('project-onto-surface', args=[iname, oname],
				opts={'labels': mask, 'name': name, 'dilation-radius': .5,
					  'point-data': False, 'cell-data': True, 'fill': True,
					  'max-hole-size': 100, 'min-size': 1000, 'smooth': 2})
			run('erode-scalars',  args=[oname, oname], opts={'cell-data': name, 'iterations': 2})
			run('calculate-element-wise', args=[oname],
				opts=[('cell-data', region_id_array), ('label', 7), ('set', 0), ('pad', 1),
					  'reset-mask', ('mul', name), ('out', oname, 'binary', name)])
			run('calculate-surface-attributes', args=[oname, oname],
				opts={'cell-labels': region_id_array, 'border-mask': 'RegionBorder'})
			run('calculate-element-wise', args=[oname],
				opts=[('cell-data', 'RegionBorder'),
					  ('threshold', 1), ('set', 0), ('pad', 1), 'reset-mask',
					  ('mul', region_id_array), ('mul', name), ('binarize', 1),
					  ('out', oname, 'char', name)])
			run('evaluate-surface-mesh', args=[oname, oname], opts={'where': name, 'gt': 0})
			run('calculate-element-wise', args=[oname],
				opts=[('cell-data', 'ComponentId'), ('binarize', 1, 2),
					  ('mul', name), ('out', oname, 'binary', name)])
			run('calculate-element-wise', args=[oname],
				opts=[('cell-data', region_id_array),
						  ('label', 1), ('mul', name), ('threshold-gt', 0), ('add', 5),
					  'reset-mask',
						  ('label', 2), ('mul', name), ('threshold-gt', 0), ('add', 6),
					  ('out', oname)])
			del_mesh_attr(oname, oname, name=['BoundaryMask', 'ComponentId', 'RegionBorder', 'DuplicatedMask'])
	return oname

# ------------------------------------------------------------------------------
def append_surfaces(name, surfaces, merge=True, tol=1e-6):
	"""Append surface meshes into single mesh file.

	Parameters
	----------
	name : str
		File path of output mesh.
	surfaces : list
		List of file paths of surface meshes to combine.
	merge : bool
		Whether to merge points of surface meshes.
	tol : float
		Maximum distance between points to be merged.

	Returns
	-------
	aname : str
		Absolute file path of output mesh.

	"""
	name = os.path.abspath(name)
	makedirs(name)
	args = surfaces
	args.append(name)
	opts = {}
	if merge: opts['merge'] = tol
	run('convert-pointset', args=args, opts=opts)
	return name

# ------------------------------------------------------------------------------
def white_refinement_mask(name, subcortex_mask):
	"""Create image foreground mask used for reconstruction of WM/cGM surface."""
	with output(invert_mask(subcortex_mask, nextname(name)), delete=True) as mask:
		close_image(mask, name, iterations=2, connectivity=18)
	return name

# ------------------------------------------------------------------------------
def deform_mesh(iname, oname=None, temp=None, opts={}):
	"""Deform mesh with the specified options."""
	if not temp:
		temp = os.path.dirname(iname)
	if not oname:
		oname = os.path.join(temp, os.path.basename(nextname(iname)))
	if force:
		try_remove(oname)
	if not os.path.isfile(oname):
		base = splitext(os.path.basename(oname))[0]
		debug_prefix = os.path.join(temp, base + '-')
		if debug > 1:
			opts['debug'] = (debug-2 if debug > 3 else 1)
			opts['debug-interval'] = (1 if debug > 2 else 10)
			opts['debug-prefix'] = debug_prefix
			opts['level-prefix'] = False
		fname_prefix = debug_prefix + 'output_'
		for fname in os.listdir(temp):
			if fname.startswith(fname_prefix):
				try_remove(os.path.join(temp, fname))
		makedirs(oname)
		run('deform-mesh', args=[iname, oname], opts=opts)
	return oname

# ------------------------------------------------------------------------------
def extract_surface(iname, oname, labels, array=_region_id_array):
	"""Extract sub-surface from combined surface mesh."""
	makedirs(oname)
	opts = [('normals', True), ('where', array), 'or']
	opts.extend([('eq', label) for label in labels])
	run('extract-pointset-cells', args=[iname, oname], opts=opts)

# ==============================================================================
# image segmentation
# ==============================================================================

# ------------------------------------------------------------------------------
def binarize(name, segmentation, labels=[], image=None, interp='linear', threshold=.5, temp=None):
	"""Make binary mask from label image.

	Parameters
	----------
	name : str
		Path of output image file.
	segmentation : str
		Path of segmentation label image file.
	labels : list
		List of segmentation labels.
	image : str, optional
		Path of reference image file. When specified, the binary image is
		resampled on this reference image grid using the specified
		interpolation followed by thresholding using the given `threshold`
		value unless `interp='nn'`.
	threshold : float
		Threshold used to binarize a non-nearest neighbor interpolated
		resampled binary image when reference `image` is given.

	Returns
	-------
	mask : str,
		Absolute path of output image file.

	"""
	if debug > 0:
		assert name, "Invalid 'name' argument"
		assert segmentation, "Invalid 'segmentation' argument"
	name = os.path.abspath(name)
	if temp:
		temp = os.path.abspath(temp)
	else:
		temp = os.path.dirname(name)
	with ExitStack() as stack:
		if image:
			mask = push_output(stack, nextname(name))
		else:
			mask = name
		makedirs(mask)
		# write binary image
		if not isinstance(labels, int) and len(labels) == 0:
			opts=[('mask', 0)]
		else:
			opts=[('label', labels)]
		opts.extend([('set', 1), ('pad', 0), ('out', mask, 'binary')])
		run('calculate-element-wise', args=[segmentation], opts=opts)
		# resample mask on target image grid
		if image:
			if interp == 'nn':
				run('transform-image', args=[mask, name], opts={'target': image, 'dofin': 'Id', 'interp': interp, 'type': 'binary'})
			else:
				resampled_mask = push_output(stack, nextname(mask))
				run('transform-image', args=[mask, resampled_mask], opts={'target': image, 'dofin': 'Id', 'interp': interp, 'type': 'float'})
				run('calculate-element-wise', args=[resampled_mask], opts=[('threshold-le', threshold), ('set', 1), ('out', [name, 'binary'])])
	return name

# ------------------------------------------------------------------------------
def binarize_cortex(regions, name=None, temp=None):
	"""Make binary cortex mask from regions label image."""
	if debug > 0:
		assert regions, "Invalid 'regions' argument"
	if not name:
		name = 'cortex-mask.nii.gz'
	if not temp:
		temp = os.path.dirname(regions)
	name = os.path.join(temp, name)
	return binarize(name=name, segmentation=regions, labels=1)

# ------------------------------------------------------------------------------
def binarize_white_matter(regions, hemisphere=Hemisphere.Unspecified, name=None, temp=None):
	"""Make binary white matter mask from regions label image.

	Parameters
	----------
	name : str
		Relative or absolute path of output mask.
	regions : str
		File path of regions image with right and left cerebrum labelled.
	hemisphere : Hemisphere
		Which hemisphere of the cerebral cortex to binarize.

	Returns
	-------
	aname : str
		Absolute path of binary white matter mask.

	"""
	if debug > 0:
		assert regions, "Invalid 'regions' argument"
	if not name:
		suffix = hemi2str(hemisphere)
		if suffix != '':
			suffix = '-' + suffix
		name = 'white{}-mask.nii.gz'.format(suffix)
	if not temp:
		temp = os.path.dirname(regions)
	name = os.path.join(temp, name)
	if   hemisphere == Hemisphere.Right: labels = 2
	elif hemisphere == Hemisphere.Left:  labels = 3
	else:								labels = [2, 3]
	return binarize(name=name, segmentation=regions, labels=labels)

# ------------------------------------------------------------------------------
def binarize_brainstem_plus_cerebellum(regions, name=None, temp=None):
	"""Make binary brainstem plus cerebellum mask from regions label image."""
	if debug > 0:
		assert regions, "Invalid 'regions' argument"
	if not name:
		name = 'brainstem+cerebellum-mask.nii.gz'
	if not temp:
		temp = os.path.dirname(regions)
	name = os.path.join(temp, name)
	makedirs(name)
	return binarize(name=name, segmentation=regions, labels=[4, 5])

# ------------------------------------------------------------------------------
def subdivide_brain(name, segmentation, white_labels, cortex_labels, right_labels, left_labels,
					subcortex_labels=[], subcortex_closing=5,
					brainstem_labels=[], brainstem_closing=5,
					cerebellum_labels=[], cerebellum_closing=5,
					tissues=None, brain_mask=None, merge_bs_cb=True,
					cortical_hull_dmap=None, temp=None):
	"""Subdivide brain into major building blocks, nicely cut disjoint regions.
	
	Parameters
	----------
	name : str
		Path of output brain regions image file.
	segmentation : str
		Path of brain segmentation label image file. May be a structural or tissue
		segmentation. If a separate tissue segmentation label image should be used
		for delineating the white matter and gray matter, provide also a `tissues` image.
	white_labels : int, str, list
		Label(s) of `segmentation` corresponding to interior of white matter surface.
		When a `tissues` image is specified, these labels are considered to
		be corresponding to the white matter tissue label instead.
	cortex_labels : int, str, list
		Label(s) of `segmentation` corresponding to cortical gray matter.
		Can contain strings of label ranges such as '1..4' instead of int's only.
		When a `tissues` image is specified, these labels are considered to
		be corresponding to the gray matter tissue label instead.
	right_labels : int, str, list
		Label(s) of `segmentation` corresponding to structures of right hemisphere.
	left_labels : int, str, list
		Label(s) of `segmentation` corresponding to structures of left hemisphere.
	subcortex_labels : int, str, list, optional
		Label(s) of `segmentation` corresponding to subcortical structures.
	subcortex_closing : int
		Number of iterations of morphological closing applied to subcortical structures mask.
	brainstem_labels : int, str, list, optional
		Label(s) of `segmentation` corresponding to brainstem.
	brainstem_closing : int
		Number of iterations of morphological closing applied to brainstem mask.
	cerebellum_labels : int, str, list, optional
		Label(s) of `segmentation` corresponding to cerebellum.
	cerebellum_closing : int
		Number of iterations of morphological closing applied to cerebellum mask.
	tissues : str
		Path of tissue segmentation used for white and gray matter
		instead of the (structural) brain `segmentation`.
	brain_mask : str
		Path of brain extraction mask used to ensure that the resampled
		image grid contains the whole of the brain.
	merge_bs_cb : bool
		Whether to merge brainstem and cerebellum labels into a single output region.
	cortical_hull_dmap : str, optional
		Path of output cortical hull distance map image file.
	temp : str
		Path of temporary directory for segmentation masks of white or gray
		matter tissue, respectively, when `tissues` segmentation provided.

	"""
	name = os.path.abspath(name)
	if temp:
		temp = os.path.abspath(temp)
	else:
		temp = os.path.dirname(name)
	white_labels	  = iterable_or_int_to_list(white_labels)
	cortex_labels	 = iterable_or_int_to_list(cortex_labels)
	right_labels	  = iterable_or_int_to_list(right_labels)
	left_labels	   = iterable_or_int_to_list(left_labels)
	subcortex_labels  = iterable_or_int_to_list(subcortex_labels)
	brainstem_labels  = iterable_or_int_to_list(brainstem_labels)
	cerebellum_labels = iterable_or_int_to_list(cerebellum_labels)
	with ExitStack() as stack:
		if tissues:
			white_mask  = os.path.join(temp, 'inner-surface-mask.nii.gz')
			cortex_mask = os.path.join(temp, 'outer-surface-mask.nii.gz')
			white_labels  = push_output(stack, binarize(name=white_mask,  segmentation=tissues, labels=white_labels))
			cortex_labels = push_output(stack, binarize(name=cortex_mask, segmentation=tissues, labels=cortex_labels))
		opts = {
			'rh': right_labels, 'lh': left_labels,
			'wm': white_labels, 'gm': cortex_labels,
			'sb': subcortex_labels,
			'bs': brainstem_labels,
			'cb': cerebellum_labels,
			'bs+cb': merge_bs_cb,
			'subcortical-closing': subcortex_closing,
			'brainstem-closing': brainstem_closing,
			'cerebellum-closing': cerebellum_closing
		}
		if brain_mask:
			opts['fg'] = os.path.abspath(brain_mask)
		if len(subcortex_labels)  > 0:
			opts['sb'] = subcortex_labels
		if len(brainstem_labels)  > 0:
			opts['bs'] = brainstem_labels
		if len(cerebellum_labels) > 0:
			opts['cb'] = cerebellum_labels
		if cortical_hull_dmap:
			cortical_hull_dmap = os.path.abspath(cortical_hull_dmap)
			opts['output-inner-cortical-distance'] = cortical_hull_dmap
			makedirs(cortical_hull_dmap)
		makedirs(name)
		run('subdivide-brain-image', args=[segmentation, name], opts=opts)
	return name

# ==============================================================================
# surface reconstruction
# ==============================================================================

# ------------------------------------------------------------------------------
def recon_boundary(name, mask, blur=1, edge_length=1, temp=None):
	"""Reconstruct surface of mask."""
	if debug > 0:
		assert name, "Invalid 'name' argument"
		assert mask, "Invalid 'mask' argument"
	makedirs(name)
	with ExitStack() as stack:
		dmap = push_output(stack, calculate_distance_map(mask, temp=temp))
		surf = push_output(stack, extract_isosurface(dmap, blur=blur))
		return remesh_surface(surf, oname=name, edge_length=edge_length)

# ------------------------------------------------------------------------------
def recon_brain_surface(name, mask, temp=None):
	"""Reconstruct surface of brain mask."""
	return recon_boundary(name, mask, blur=2, temp=temp)

# ------------------------------------------------------------------------------
def recon_brainstem_plus_cerebellum_surface(name, mask=None, regions=None,
											region_id_array=_region_id_array,
											cortex_mask_array=_cortex_mask_array,
											temp=None):
	"""Reconstruct surface of merged brainstem plus cerebellum region."""
	if debug > 0:
		assert mask or regions, "Either 'regions' or 'mask' argument required"
	name = os.path.abspath(name)
	if temp:
		temp = os.path.abspath(temp)
	else:
		temp = os.path.dirname(name)
	with ExitStack() as stack:
		if not mask:
			mask = push_output(stack, binarize_brainstem_plus_cerebellum(regions, temp=temp))
		mesh = push_output(stack, nextname(name, temp=temp))
		recon_boundary(mesh, mask, blur=1, temp=temp)
		if region_id_array:
			run('project-onto-surface', args=[mesh, mesh],
				opts={'constant': 7, 'name': region_id_array, 'type': 'short',
					  'point-data': False, 'cell-data': True})
		if cortex_mask_array:
			run('project-onto-surface', args=[mesh, mesh],
				opts={'constant': 0, 'name': cortex_mask_array, 'type': 'uchar',
					  'point-data': False, 'cell-data': True})
		rename(mesh, name)
	return name

# ------------------------------------------------------------------------------
def recon_cortical_surface(name, mask=None, regions=None,
						   hemisphere=Hemisphere.Unspecified,
						   corpus_callosum_mask=None, temp=None, use_fast_collision=False):
	"""Reconstruct initial surface of right and/or left cerebrum.
 
	The input is either a regions label image generated by `subdivide-brain-image`
	or a custom binary mask created from a given brain segmentation. The boundary of
	the mask must be sufficiently close to the WM/cGM interface.

	Attention: Order of arguments may differ from the order of parameter help below!
			   Pass parameter values as keyword argument, e.g., regions='regions.nii.gz'.

	Parameters
	----------
	name : str
		File path of output surface file.
	mask : str, optional
		File path of binary mask to use instead of `regions` segments.
		When this argument is given, both `regions` and `hemisphere` are ignored.
	regions : str
		File path of regions image with right and left cerebrum labelled.
		Ignored when a custom `mask` image is given instead.
	hemisphere : Hemisphere
		Which hemisphere of the cerebral cortex to reconstruct. If Hemisphere.Both, 
		the surfaces of both hemispheres are reconstructed and the two file paths
		returned as a 2-tuple. When Hemisphere.Unspecified, the hemisphere is
		determined from the output file `name` if possible.
	corpus_callosum_mask : str, optional
		File path of binary mask of segmented corpus callosum.
		Used to disable clustering based hole filling of implicit
		surface distance for points nearby the corpus callosum.
	temp : str
		Path of temporary working directory. Intermediate files are written to
		this directory and deleted on exit unless the global `debug` flag is set.
		When not specified, intermediate files are written to the same directory
		as the output mesh file, i.e., the directory of `name`.

	Returns
	-------
	path : str, tuple
		Absolute path of output surface mesh file.
		A 2-tuple of (right, left) surface mesh file paths is returned when
		no `mask` but a `regions` image was given with `hemisphere=Hemisphere.Both`.

	"""
	if debug > 0:
		assert name, "Invalid 'name' argument"
		assert mask or regions, "Either 'regions' or 'mask' argument required"
	name = os.path.abspath(name)
	(base, ext) = splitext(name)
	if not ext:
		ext  = '.vtp'
		name = base + ext
	if temp:
		temp = os.path.abspath(temp)
	else:
		temp = os.path.dirname(name)
	if base.endswith('-lh') or base.endswith('.lh') or base.startswith('lh.'):
		if base.startswith('lh.'):
			base = base[3:]
		else:
			base = base[0:-3]
		if hemisphere == Hemisphere.Right:
			sys.stderr.write("Warning: Cortical surface file name suggests left hemisphere, but right hemisphere is being reconstructed!\n")
		elif hemisphere == Hemisphere.Unspecified:
			hemisphere = Hemisphere.Left
	elif base.endswith('-rh') or base.endswith('.rh') or base.startswith('rh.'):
		if base.startswith('rh.'):
			base = base[3:]
		else:
			base = base[0:-3]
		if hemisphere == Hemisphere.Left:
			sys.stderr.write("Warning: Cortical surface file name suggests right hemisphere, but left hemisphere is being reconstructed!\n")
		elif hemisphere == Hemisphere.Unspecified:
			hemisphere = Hemisphere.Right
	if not mask and hemisphere == Hemisphere.Both:
		rname = recon_cortical_surface('{}-{}{}'.format(base, hemi2str(Hemisphere.Right), ext),
									   regions=regions, hemisphere=Hemisphere.Right,
									   corpus_callosum_mask=corpus_callosum_mask, temp=temp)
		lname = recon_cortical_surface('{}-{}{}'.format(base, hemi2str(Hemisphere.Left), ext),
									   regions=regions, hemisphere=Hemisphere.Left,
									   corpus_callosum_mask=corpus_callosum_mask, temp=temp)
		return (rname, lname)
	base = os.path.basename(base)
	hemi = hemi2str(hemisphere)
	with ExitStack() as stack:
		if not mask:
			mask = '{}-{}-mask.nii.gz'.format(base, hemi)
			mask = push_output(stack, binarize_white_matter(regions, name=mask, hemisphere=hemisphere, temp=temp))
		dmap = push_output(stack, calculate_distance_map(mask, temp=temp))
		hull = push_output(stack, extract_convex_hull(dmap, temp=temp))

		opts = {
			'implicit-surface': dmap,
			'distance': 1,
				'distance-measure': 'normal',
				'distance-threshold': 2,
				'distance-max-depth': 5,
				'distance-hole-filling': True,
			'curvature': .8,
			'gauss-curvature': .4,
				'gauss-curvature-outside': 1,
				'gauss-curvature-minimum': .1,
				'gauss-curvature-maximum': .5,
			'repulsion': 4,
				'repulsion-distance': .5,
				'repulsion-width': 1,
			'optimizer': 'EulerMethod',
				'step': [.5, .1],
				'steps': [300, 200],
				'epsilon': 1e-8,
				'extrinsic-energy': True,
				'delta': 1e-2,
				'min-active': '1%',
				'min-width': .2,
				'min-distance': .5,
				'fast-collision-test': use_fast_collision,
				'non-self-intersection': True,
			'remesh': 1,
				'min-edge-length': .5,
				'max-edge-length': 1,
				'triangle-inversion': True,
				'reset-status': True
		}

		if corpus_callosum_mask:
			out1 = nextname(name, temp=temp)
			out1 = push_output(stack, add_corpus_callosum_mask(hull, mask=corpus_callosum_mask, oname=out1))
			out2 = nextname(out1)
		else:
			out1 = hull
			out2 = nextname(name, temp=temp)
		out2 = push_output(stack, deform_mesh(out1, out2, opts=opts))
		if corpus_callosum_mask:
			out3 = push_output(stack, del_corpus_callosum_mask(out2))
		else:
			out3 = out2
		remove_intersections(out3, oname=name, max_attempt=10)
	return name

# ------------------------------------------------------------------------------
def join_cortical_surfaces(name, regions, right_mesh, left_mesh, bs_cb_mesh=None,
						   region_id_array=_region_id_array, cortex_mask_array=_cortex_mask_array,
						   internal_mesh=None, temp=None, check=True, join_tol = 1):
	"""Join cortical surfaces of right and left hemisphere at medial cut.

	Optionally, the brainstem plus cerebellum surface can be joined with the
	resulting surface mesh as well at the brainstem cut. This, however, currently
	fails for a significant number of cases to find a cutting plane which intersects
	the joined surface mesh with one closed curve near the superior brainstem cut
	and is therefore not recommended.

	Attention: Order of arguments may differ from the order of parameter help below!
			   Pass parameter values as keyword argument, e.g., regions='regions.nii.gz'.

	Parameters
	----------
	name : str
		Path of output surface mesh file.
	regions : str
		Path of regions image file with right and left cerebrum labelled.
		The boundary between right and left cerebrum segments defines the cut
		where the two surfaces are merged. Similarly for the brainstem cut
		when `bs_cb_mesh` surface given.
	right_mesh : str
		Path of right surface mesh file.
	left_mesh : str
		Path of left surface mesh file.
	bs_cb_mesh : str, optional
		Path of brainstem plus cerebellum surface mesh file.
	region_id_array : str, optional
		Name of cell data array with labels corresponding to the volumetric `regions`
		label of which the surface is the boundary of. Required for `cortex_mask_array`.
		Note: Required by subsequent cortical surface reconstruction steps.
	cortex_mask_array : str, optional
		Name of cortex mask cell data array. If None or an empty string,
		no cortex mask data array is added to the output surface file.
		Note: Required by subsequent cortical surface reconstruction steps.
	internal_mesh : str, optional
		Path of output surface mesh file with cutting plane cross sections
		that divide the interior of the merged surfaces into right, left,
		and brainstem plus cerebellum (when `bs_cb_mesh` surface given).
	temp : str
		Path of temporary working directory. Intermediate files are written to
		this directory and deleted on exit unless the global `debug` flag is set.
		When not specified, intermediate files are written to the same directory
		as the output mesh file, i.e., the directory of `name`.
	check : bool
		Check topology and consistency of output surface mesh.

	Returns
	-------
	path : str
		Absolute path of output surface mesh file.

	"""

	if debug > 0:
		assert name, "Output file 'name' required"
		assert regions, "Input 'regions' label image required"
		assert right_mesh, "Input 'right_mesh' required"
		assert left_mesh, "Input 'left_mesh' required"

	join_bs_cb = True # deprecated option, always True when bs_cb_mesh not None
 
	name = os.path.abspath(name)
	if not temp:
		temp = os.path.dirname(name)
	else:
		temp = os.path.abspath(temp)

	if region_id_array:
		region_id_array_name = region_id_array
	else:
		region_id_array_name = _region_id_array

	with ExitStack() as stack:
		# merge surface meshes
		joined = push_output(stack, nextname(name, temp=temp))
		if force or not os.path.isfile(joined):
			surfaces = [right_mesh, left_mesh]
			if bs_cb_mesh and join_bs_cb:
				surfaces.append(bs_cb_mesh)
			run('merge-surfaces',
				opts={'input': surfaces, 'output': joined, 'labels': regions, 'source-array': region_id_array_name,
					  'tolerance': join_tol, 'largest': True, 'dividers': (internal_mesh != None), 'snap-tolerance': .1,
					  'smoothing-iterations': 100, 'smoothing-lambda': 1})
			if bs_cb_mesh and join_bs_cb:
				run('calculate-element-wise', args=[joined], opts=[('cell-data', region_id_array_name),
																   ('map', (-1, -3), (-2, -1), (-3, -2), (3, 7)),
																   ('out', joined)])
			del_mesh_attr(joined, pointdata=region_id_array_name)
		# check topology of joined surface mesh
		if check:
			info = evaluate_surface(joined, mesh=True, topology=True)
			num_boundaries = get_num_boundaries(info)
			if num_boundaries != 0:
				raise Exception('Merged surface is non-closed, no. of boundaries: {}'.format(num_boundaries))
			euler = get_euler_characteristic(info)
			if internal_mesh:
				if bs_cb_mesh: expect_euler = 4
				else:		  expect_euler = 3
			else:			  expect_euler = 2
			if euler != expect_euler:
				raise Exception('Merged surface with dividers has unexpected Euler characteristic: {} (expected {})'.format(euler, expect_euler))
		# ensure there are no self-intersections of the joined surface mesh
		checked = push_output(stack, nextname(joined))
		if force or not os.path.isfile(checked):
			remove_intersections(joined, checked, max_attempt=10, smooth_iter=10)
		joined = checked
		# when brainstem+cerebellum surface given, but it should not be joined with the
		# cerebrum surface, remove triangles near brainstem cut and any triangles of the
		# brainstem+cerebellum surface that intersect with the cerebrum surface
		if bs_cb_mesh and not join_bs_cb:
			joined_with_bscb = push_output(stack, nextname(joined))
			run('merge-surfaces',
				opts={'input': [joined, bs_cb_mesh], 'output': joined_with_bscb, 'labels': regions,
					  'source-array': region_id_array_name, 'tolerance': .5, 'join': False})
			run('calculate-element-wise', args=[joined_with_bscb], opts=[('cell-data', region_id_array_name), ('map', (3, 7)), ('out', joined_with_bscb)])
			modified_bscb = push_output(stack, nextname(joined_with_bscb))
			check_intersections(joined_with_bscb, modified_bscb)
			run('calculate-element-wise', args=[modified_bscb],
				opts=[('cell-data', region_id_array_name), ('label', 7), ('set', 1), ('pad', 0),
					  ('mul', 'CollisionType'), ('binarize', 0, 0), ('out', modified_bscb, 'binary', 'SelectionMask')])
			run('extract-pointset-cells', args=[modified_bscb, modified_bscb], opts=[('where', 'SelectionMask'), ('eq', 0)])
			del_mesh_attr(modified_bscb, pointdata='CollisionMask', celldata=['SelectionMask', 'CollisionType'])
			joined_with_modified_bscb = push_output(stack, nextname(modified_bscb))
			append_surfaces(joined_with_modified_bscb, surfaces=[joined, modified_bscb], merge=True, tol=0)
			joined = joined_with_modified_bscb
		# optionally, add cortex mask highlighting cells which are nearby cGM
		# this mask contains exactly two components, a right and a left cortex
		if cortex_mask_array:
			with output(binarize_cortex(regions, temp=temp), delete=True) as mask:
				joined = push_output(stack, add_cortex_mask(joined, mask, name=cortex_mask_array, region_id_array=region_id_array_name))
			if check:
				info = evaluate_surface(joined, mesh=True, opts=[('where', cortex_mask_array), ('gt', 0)])
				num_components = get_num_components(info)
				if num_components != 2:
					raise Exception("No. of cortex mask components: {} (expected 2)".format(num_components))
		# save divider(s) as separate surface mesh such that the inside of the
		# white surface is clearly defined when converting it to a binary mask
		# during the image edge-based refinement step below
		if internal_mesh:
			joined_without_dividers = push_output(stack, nextname(joined))
			internal_mesh = os.path.abspath(internal_mesh)
			run('extract-pointset-cells', args=[joined, internal_mesh],		   opts=[('where', region_id_array_name), ('lt', 0)])
			run('extract-pointset-cells', args=[joined, joined_without_dividers], opts=[('where', region_id_array_name), ('gt', 0)])
			joined = joined_without_dividers
		# remove RegionId array if not desired by caller
		if not region_id_array:
			if internal_mesh:
				del_mesh_attr(internal_mesh, celldata=region_id_array_name)
			del_mesh_attr(joined, celldata=region_id_array_name)
		rename(joined, name)
	if internal_mesh:
		return (name, internal_mesh)
	return name

# ------------------------------------------------------------------------------
def recon_white_surface(name, t2w_image, wm_mask, gm_mask, cortex_mesh, bs_cb_mesh=None,
						cortex_mask_array=_cortex_mask_array, region_id_array=_region_id_array,
						t1w_image=None, cortical_hull_dmap=None, segmentation=None,
						subcortex_labels =[], subcortex_mask =None,
						ventricles_labels=[], ventricles_mask=None, ventricles_dmap=None,
						cerebellum_labels=[], cerebellum_mask=None, cerebellum_dmap=None,
						temp=None, check=True, use_fast_collision=False):
	"""Reconstruct white surface based on WM/cGM image edge distance forces.

	This step refines the initial surface mesh. In areas where the image
	edges are weak, no change of position should be enforced. Otherwise, move
	the surface points a limited distance from their initial position to a
	location along the normal direction where the image gradient is stronger.
	The resulting surface should delinate the WM/cGM interface at least as
	good as the initial surface mesh, but likely better in many places.

	Attention: Order of arguments may differ from the order of parameter help below!
			   Pass parameter values as keyword argument, e.g., wm_mask='wm.nii.gz'.

	Parameters
	----------
	name : str
		Path of output surface mesh file.
	cortex_mesh : str
		Path of initial cortical surface mesh file. This surface is reconstructed
		from a given brain segmentation with possible errors such as in particular
		CSF mislabelled as either WM or GM.
	bs_cb_mesh : str, optional
		Path of brainstem plus cerebellum mesh file.
	cortex_mask_array : str
		Name of cortex mask cell data array. Only nodes adjacent to only cortical
		faces are deformed. Non-cortical surface nodes remain unchanged.

	t1w_image : str, optional
		Path of T1-weighted intensity image file.
	t2w_image : str
		Path of T2-weighted intensity image file.
	wm_mask : str
		Path of binary WM segmentation image file.
	gm_mask : str
		Path of binary cGM segmentation image file.

	segmentation : str, optional
		Structural brain segmentation label image used to extract subcortical
		structures, ventricles, and cerebellum when structure labels instead
		of a pre-computed binary mask or structure distance image, respectively,
		is given. Otherwise, this image path is unused.
	subcortex_labels : int, list, optional
		Labels of `segmentation` labels corresponding to subcortical structures.
		Used to create a binary mask image when no `subcortex_mask` is given.
	subcortex_mask : str, optional
		Path of subcortical structures mask file. These structures are excluded
		from the image foreground such that the image-based edge distance
		force is not mislead by a WM/dGM edge of a subcortical structure.
	ventricles_labels : int, list, optional
		Labels of `segmentation` labels corresponding to lateral ventricles.
		Used to create a binary mask image when no `ventricles_mask` is given.
	ventricles_mask : str, optional
		Path of binary mask image file of lateral ventricles.
	ventricles_dmap : str, optional
		Path of distance map computed from lateral ventricles segment.
	cerebellum_labels : int, list, optional
		Labels of `segmentation` labels corresponding to cerebellum.
		Used to create a binary mask image when no `cerebellum_mask` is given.
	cerebellum_mask : str, optional
		Path of binary mask image file of cerebellum.
	cerebellum_dmap : str, optional
		Path of distance map computed from cerebellum segment.
	cortical_hull_dmap : str, optional
		Path of distance image from each voxel to the cortical hull with
		positive values inside the cortical hull. This image is an optional
		output of the `subdivide-brain-image` tool.

	temp : str
		Path of temporary working directory. Intermediate files are written to
		this directory and deleted on exit unless the global `debug` flag is set.
		When not specified, intermediate files are written to the same directory
		as the output mesh file, i.e., the directory of `name`.
	check : bool
		Check final surface mesh for self-intersections and try to remove these.

	Returns
	-------
	path : str
		Absolute path of output surface mesh file.

	"""

	if debug > 0:
		assert name, "Output file 'name' required"
		assert t2w_image, "T2-weighted intensity image required"
		assert wm_mask, "White matter segmentation mask required"
		assert gm_mask, "Gray matter segmentation mask required"
		assert cortex_mesh, "Initial cortical surface mesh required"

	name = os.path.abspath(name)
	(base, ext) = splitext(name)
	if not ext:
		ext  = '.vtp'
		name = base + ext
	if not temp:
		temp = os.path.dirname(name)
	else:
		temp = os.path.abspath(temp)

	with ExitStack() as stack:

		init_mesh = push_output(stack, nextname(name, temp=temp))

		# optionally, append brainstem plus cerebellum
		if bs_cb_mesh:
			append_surfaces(init_mesh, surfaces=[cortex_mesh, bs_cb_mesh], merge=False)
			cortex_mesh = init_mesh

		# initialize node status
		status_array = 'InitialStatus'
		run('copy-pointset-attributes', args=[cortex_mesh, cortex_mesh, init_mesh],
			opts={'celldata-as-pointdata': [cortex_mask_array, status_array, 'other', 'binary'], 'unanimous': None})
		run('erode-scalars', args=[init_mesh, init_mesh], opts={'array': status_array, 'iterations': 8})

		# deform surface towards WM/cGM image edges
		opts={'image': t2w_image,
				  'wm-mask': wm_mask,
				  'gm-mask': gm_mask,
			  'edge-distance': 1,
				  'edge-distance-type': 'Neonatal T2-w WM/cGM',
				  'edge-distance-max-depth': 5,
				  'edge-distance-median': 1,
				  'edge-distance-smoothing': 1,
				  'edge-distance-averaging': [4, 2, 1],
			  'curvature': 2,
			  'gauss-curvature': .5,
				  'gauss-curvature-outside': .5,
				  'gauss-curvature-minimum': .1,
				  'gauss-curvature-maximum': .2,
			  'optimizer': 'EulerMethod',
				  'step': .2,
				  'steps': [50, 100],
				  'extrinsic-energy': True,
				  'epsilon': 1e-6,
				  'delta': .01,
				  'min-active': '1%',
				  'reset-status': True,
				  'non-self-intersection': True,
				  'fast-collision-test': use_fast_collision,
				  'min-width': .1,
			  'repulsion': 4,
				  'repulsion-distance': .5,
				  'repulsion-width': 1,
			  'remesh': 1,
				  'min-edge-length': .5,
				  'max-edge-length': 1}
		if t1w_image:
			opts['t1w-image'] = t1w_image
		if cortical_hull_dmap:
			opts['inner-cortical-distance-image'] = cortical_hull_dmap
		if (segmentation and ventricles_labels) or ventricles_mask or ventricles_dmap:
			if not ventricles_dmap:
				if not ventricles_mask:
					ventricles_mask = os.path.join(temp, 'ventricles-mask.nii.gz')
					ventricles_mask = push_output(stack, binarize(name=ventricles_mask, segmentation=segmentation, labels=ventricles_labels))
				ventricles_dmap = push_output(stack, calculate_distance_map(ventricles_mask, temp=temp))
			opts['ventricles-distance-image'] = ventricles_dmap
		if (segmentation and cerebellum_labels) or cerebellum_mask or cerebellum_dmap:
			if not cerebellum_dmap:
				if not cerebellum_mask:
					cerebellum_mask = os.path.join(temp, 'cerebellum-mask.nii.gz')
					cerebellum_mask = push_output(stack, binarize(name=cerebellum_mask, segmentation=segmentation, labels=cerebellum_labels))
				cerebellum_dmap = push_output(stack, calculate_distance_map(cerebellum_mask, temp=temp))
			opts['cerebellum-distance-image'] = cerebellum_dmap
		if (segmentation and subcortex_labels) or subcortex_mask:
			if not subcortex_mask:
				subcortex_mask = os.path.join(temp, 'subcortex-mask.nii.gz')
				subcortex_mask = push_output(stack, binarize(name=subcortex_mask, segmentation=segmentation, labels=subcortex_labels))
			mask = os.path.join(temp, os.path.basename(base) + '-foreground.nii.gz')
			opts['mask'] = push_output(stack, white_refinement_mask(mask, subcortex_mask))
		mesh = push_output(stack, deform_mesh(init_mesh, opts=opts))
		if bs_cb_mesh:
			run('extract-pointset-cells', args=[mesh, mesh], opts=[('where', region_id_array), ('ne', 7)])

		# smooth white surface mesh
		smooth = push_output(stack, smooth_surface(mesh, iterations=100, lambda_value=.33, mu=-.34, weighting='combinatorial', excl_node=True, mask=status_array))

		# remove intersections if any
		if check:
			remove_intersections(smooth, oname=name, mask=status_array)

		# write output mesh
		del_mesh_attr(smooth, pointdata=status_array)
		rename(smooth, name)

	return name

# ------------------------------------------------------------------------------
def recon_pial_surface(name, t2w_image, wm_mask, gm_mask, white_mesh,
					   bs_cb_mesh=None, brain_mask=None, remesh=0, outside_white_mesh=True,
					   region_id_array=_region_id_array, cortex_mask_array=_cortex_mask_array,
					   temp=None, check=True, use_fast_collision = False):
	"""Reconstruct pial surface based on cGM/CSF image edge distance forces.

	When the pial surface is not allowed to intersect the white surface mesh
	(i.e., `outside_white_surface=True`), the following steps are performed:
	1. Deform cortical nodes outside up to a maximum distance from the white surface
	   - Hard non-self-intersection constraint disabled
	2. Blend between initial and outside white surface node position
	3. Merge white and pial surfaces into single connected surface mesh
	   - If white and pial surfaces intersect afterwards, try to
		 remove those by smoothing the surface near the intersections
	4. Deform pial surface nodes towards cGM/CSF image edges
	   - Hard non-self-intersection constraint enabled
	5. Remove white surface mesh from output mesh file

	Otherwise, only step 4 is performed to deform the input white surface mesh
	towards the outer cortical surface based on image edge forces. In this case,
	the resulting pial surface mesh may potentially intersect the `white_mesh`.

	Attention: Order of arguments may differ from the order of parameter help below!
			   Pass parameter values as keyword argument, e.g., wm_mask='wm.nii.gz'.

	Parameters
	----------
	name : str
		Path of output surface mesh file.
	t2w_image : str
		Path of T2-weighted intensity image file.
	wm_mask : str
		Path of binary WM segmentation image file.
	gm_mask : str
		Path of binary cGM segmentation image file.
	white_mesh : str
		Path of white surface mesh file.
	bs_cb_mesh : str, optional
		Path of brainstem plus cerebellum surface mesh file. If specified,
		this surface mesh is added to the initial pial surface mesh for step 4
		of the pial surface reconstruction such that the pial surface does
		not deform into the brainstem plus cerebellum region.
	remesh : int
		Number of Euler steps between local remeshing is performed.
		Generally, it is better to not remesh (`remesh=0`) such that ID based
		one-to-one vertex correspondences between white and pial surface meshes
		are preserved. These vertex correspondences simplify further analysis.
	outside_white_mesh : bool
		Whether pial surface mesh is required to be strictly outside the white
		surface mesh, i.e., the pial surface may not intersect the `white_mesh`.
		This is recommended, but may currently fail for some cases when the
		blended offset surface intersects the white surface mesh after step 3.
	brain_mask : str, optional
		Path of brain extraction mask file.
	temp : str, optional
		Path of temporary working directory. Intermediate files are written to
		this directory and deleted on exit unless the global `debug` flag is set.
		When not specified, intermediate files are written to the same directory
		as the output mesh file, i.e., the directory of `name`.
	check : bool
		Check topology and consistency of intermediate surface meshes.

	Returns
	-------
	path : str
		Absolute path of output surface mesh file.

	"""

	if debug > 0:
		assert name, "Output file 'name' required"
		assert t2w_image, "T2-weighted intensity image required"
		assert wm_mask, "White matter segmentation mask required"
		assert gm_mask, "Gray matter segmentation mask required"
		assert white_mesh, "White surface mesh required"

	name = os.path.abspath(name)
	(base, ext) = splitext(name)
	if not ext:
		ext  = '.vtp'
		name = base + ext
	if not temp:
		temp = os.path.dirname(name)
	else:
		temp = os.path.abspath(temp)

	with ExitStack() as stack:

		# create foreground mask
		mask = push_output(stack, os.path.join(temp, os.path.basename(base) + '-foreground.nii.gz'))
		run('extract-pointset-surface', args=[], opts={'input': white_mesh, 'mask': mask, 'reference': t2w_image, 'outside': True})
		if brain_mask:
			run('calculate-element-wise', args=[mask], opts=[('mul', brain_mask), ('out', mask)])

		# initialize deformable surface mesh
		init_mesh = push_output(stack, nextname(name, temp=temp))
		status_array = 'Status'
		if outside_white_mesh:

			# initialize node status
			run('copy-pointset-attributes', args=[white_mesh, white_mesh, init_mesh], opts={'celldata-as-pointdata': [region_id_array, status_array], 'unanimous': None})
			run('calculate-element-wise', args=[init_mesh], opts=[('point-data', status_array), ('label', 7), ('set', 0), ('pad', 1), ('out', init_mesh, 'binary')])

			# deform pial surface outwards a few millimeters
			opts={'normal-force': 1,
				  'curvature': 1,
				  'optimizer': 'EulerMethod',
					  'step': .1,
					  'steps': 100,
					  'max-displacement': max(get_voxel_size(t2w_image)),
					  'non-self-intersection': True,
					  'fast-collision-test': use_fast_collision,
					  'min-distance': .1,
					  'min-active': '10%',
					  'delta': .0001}
			offset_mesh = push_output(stack, deform_mesh(init_mesh, opts=opts))
			if debug == 0: try_remove(init_mesh)

			# blend between original position of non-cortical points and offset surface points
			blended_mesh = push_output(stack, nextname(offset_mesh))
			run('copy-pointset-attributes', args=[offset_mesh, offset_mesh, blended_mesh], opts={'celldata-as-pointdata': cortex_mask_array, 'unanimous': None})
			run('blend-surface', args=[white_mesh, blended_mesh, blended_mesh], opts={'where': cortex_mask_array, 'gt': 0, 'smooth-iterations': 3})
			run('calculate-element-wise', args=[blended_mesh], opts=[('cell-data', region_id_array), ('map', (1, 3), (2, 4)), ('out', blended_mesh)])
			if debug == 0: try_remove(offset_mesh)

			# merge white surface mesh with initial pial surface mesh
			cortex_mesh = push_output(stack, nextname(blended_mesh))
			append_surfaces(cortex_mesh, surfaces=[white_mesh, blended_mesh], merge=True, tol=0)
			info = evaluate_surface(cortex_mesh, cortex_mesh, mesh=True, opts=[('where', cortex_mask_array), ('gt', 0)])
			run('calculate-element-wise', args=[cortex_mesh], opts=[('cell-data', cortex_mask_array), ('mask', 'BoundaryMask'), ('set', 0), ('out', cortex_mesh)])
			run('calculate-element-wise', args=[cortex_mesh], opts=[('cell-data', region_id_array),   ('mask', 'BoundaryMask'), ('add', 4), ('out', cortex_mesh)])
			num_components = get_num_components(info)
			if num_components < 2:
				raise Exception("No. of cortex mask components: {} (expected 2)".format(num_components))
			elif num_components > 2:
				run('calculate-element-wise', args=[cortex_mesh],
					opts=[('cell-data', 'ComponentId'), ('threshold-gt', 2), ('pad', 0), ('mul', cortex_mask_array),
						  ('clamp-above', 1), ('out', cortex_mesh, 'binary', cortex_mask_array)])
			del_mesh_attr(cortex_mesh, name=['ComponentId', 'BoundaryMask'], celldata='DuplicatedMask')
			if debug == 0: try_remove(blended_mesh)

			# resolve intersections between white and pial surface if any
			init_mesh = push_output(stack, nextname(cortex_mesh))
			remove_white_pial_intersections(cortex_mesh, init_mesh, region_id_array=region_id_array)
			#run('copy-pointset-attributes', args=[cortex_mesh, cortex_mesh, init_mesh], opts={'celldata-as-pointdata': cortex_mask_array, 'unanimous': None})
			#remove_intersections(init_mesh, init_mesh, mask=cortex_mask_array)
			#del_mesh_attr(init_mesh, pointdata=cortex_mask_array)
			if debug == 0: try_remove(cortex_mesh)
		else:
			run('calculate-element-wise', args=[white_mesh], opts=[('cell-data', region_id_array), ('map', (1, 3), (2, 4)), ('out', init_mesh)])

		# optionally, append brainstem plus cerebellum surface
		if bs_cb_mesh:
			append_surfaces(init_mesh, surfaces=[init_mesh, bs_cb_mesh], merge=False)

		# initialize node status for pial surface reconstruction
		run('copy-pointset-attributes', args=[init_mesh, init_mesh], opts={'celldata-as-pointdata': [region_id_array, status_array], 'unanimous': None})
		run('calculate-element-wise', args=[init_mesh], opts=[('point-data', status_array), ('label', 3, 4), ('set', 1), ('pad', 0), ('out', init_mesh, 'binary')])

		# deform pial surface towards cGM/CSF image edges
		opts={'image': t2w_image,
				  'mask': mask,
				  'wm-mask': wm_mask,
				  'gm-mask': gm_mask,
			  'edge-distance': 1,
				  'edge-distance-type': 'Neonatal T2-w cGM/CSF',
				  'edge-distance-max-depth': 5,
				  'edge-distance-median': 1,
				  'edge-distance-smoothing': 1,
				  'edge-distance-averaging': [4, 2, 1],
			  'curvature': 2,
			  'gauss-curvature': .8,
				  'gauss-curvature-inside': 2,
				  'gauss-curvature-minimum': .1,
				  'gauss-curvature-maximum': .4,
				  'negative-gauss-curvature-action': 'inflate',
			  'repulsion': 2,
				  'repulsion-distance': .5,
				  'repulsion-width': 1,
			  'optimizer': 'EulerMethod',
				  'step': .2,
				  'steps': [25, 50, 100],
				  'epsilon': 1e-6,
				  'delta': .01,
				  'min-active': '5%',
				  'reset-status': True,
				  'non-self-intersection': True,
				  'fast-collision-test': use_fast_collision,
				  'min-distance': .1,
			  'remesh': remesh,
				  'min-edge-length': .3,
				  'max-edge-length': 2,
				  # Inversion may cause edges at the bottom of sulci to change from running
				  # along the sulcus (minimum curvature direction) to be inverted to then
				  # go across the sulcus (maximum curvature direction) which iteratively
				  # contributes to smoothing out the sulci which should actually be preserved.
				  'triangle-inversion': False}
		mesh = push_output(stack, deform_mesh(init_mesh, opts=opts))
		extract_surface(mesh, name, array=region_id_array, labels=[-1, -2, -3, 3, 4, 5, 6])

	return name

# ------------------------------------------------------------------------------
def split_cortical_surfaces(joined_mesh, right_name=None, left_name=None, internal_mesh=None, temp=None):
	"""Save cortical surfaces of right and left hemispheres as separate meshes."""
	# TODO: If no `internal_mesh` mesh is avaiable, use vtkFillHolesFilter
	#	   to close the non-closed genus-0 surfaces at the medial or brainstem cut.
	#	   Possibly ensure that the new cells/points are identical for
	#	   the closed surface meshes of right and left hemisphere.
	#
	#	   See extract-pointset-surface -fill-holes
	if debug > 0:
		assert joined_mesh, "Joined cortical surface mesh required"
		assert right_name or left_name, "At least either 'right_name' or 'left_name' required"
	if not temp:
		temp = os.path.dirname(joined_mesh)
	with ExitStack() as stack:
		if internal_mesh:
			name = os.path.join(temp, os.path.splitext(joined_mesh)[0] + '+internal.vtp')
			joined_mesh = push_output(stack, append_surfaces(name, surfaces=[joined_mesh, internal_mesh], merge=True, tol=0))
		if right_name: extract_surface(joined_mesh, right_name, labels=[-1, -3, 1, 3, 5])
		if left_name:  extract_surface(joined_mesh, left_name,  labels=[-1, -2, 2, 4, 6])
