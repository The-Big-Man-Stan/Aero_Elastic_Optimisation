'''
Created on 23 Jul 2011

@author: fzhu
Modified by Anthony Stannard 2018
'''

from scipy.io import netcdf_file
import numpy,os,shutil
import netCDF4
from scipy.io import *

mesh = 'rae2822.taumesh1'
no_of_boundaries = 5
cdf_file = netcdf_file(mesh,'r')

# Get Node IDs of points and their associated boundary markers
# =================================
# GET SURFACE MESH POINT IDs and other info from the netCDF file
# =================================
points_x = cdf_file.variables['points_xc'][:]
points_y = cdf_file.variables['points_yc'][:]
points_z = cdf_file.variables['points_zc'][:]
triangle_points = numpy.ndarray.tolist(cdf_file.variables['points_of_surfacetriangles'][:])
quad_points = numpy.ndarray.tolist(cdf_file.variables['points_of_surfacequadrilaterals'][:])
surface_points = triangle_points + quad_points
surface_marker = numpy.ndarray.tolist(cdf_file.variables['boundarymarker_of_surfaces'][:])

# ================================
# GET GLOBAL SURFACE IDS BELONGING TO EACH BOUNDARY
# ================================
# gsid stands for global surface ID, this is the index of the element within the
# surface_marker array
gsid_array = []
element_points = []
NodeIDs = []
data = []

# Make gsid array be a list containing number_of_boundaries lists
for _ in range(no_of_boundaries):
	gsid_array.append([])
	element_points.append([])
	NodeIDs.append([])
	data.append([])

# Get the list of surface IDs belonging to each boundary in separate lists
for gsid,sm in enumerate(surface_marker):
	gsid_array[sm-1].append(gsid)

# ================================
# GET Points and their Global Node IDs written in file
# ================================
# Element points contain the Node IDs belonging to a shape on a boundary
for i in range(no_of_boundaries):
	for j in gsid_array[i]:
		a = surface_points[j]
		element_points[i].append(a)

for i in range(no_of_boundaries):
	for j in element_points[i]:
		for k in j:
			NodeIDs[i].append(k)
for i in range(no_of_boundaries):
	temp = []
	NodeIDs[i] = sorted(list(set(NodeIDs[i]))) # Get rid of duplicates
	for j in NodeIDs[i]:
		temp.append([points_x[j],points_y[j],points_z[j],j])
	data[i] = numpy.asarray(temp)

data = numpy.asarray(data)
       
# User input required to write appropriate descriptive names for node points
# and to which boundary they belong
data_inner0 = data[0]
data_inner1 = data[1]
data_ff = data[2]
data_us = data[3]
data_ls = data[4]

data_airfoil_xu = data_us[:,0]
data_airfoil_zu = data_us[:,2]
data_airfoil_xl = data_ls[:,0]
data_airfoil_zl = data_ls[:,2]

with open('new_airfoil_coords.txt','a') as myfile:
	myfile.write(str(data_airfoil_xu))
	myfile.write('\n')
	myfile.write(str(data_airfoil_xl))
	myfile.write('\n')
	myfile.write(str(data_airfoil_zu))
	myfile.write('\n')
	myfile.write(str(data_airfoil_zl))
