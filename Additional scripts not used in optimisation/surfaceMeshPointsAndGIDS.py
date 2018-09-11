"""
Written by Jed Hollom (August 2015)

NOTE: This script is only compatable with netCDF
 files written using read_taumesh.py!

Reads in a given netcdf file and:
1. Prints a list of all dimensions with their number
2. Prints a list of all variables with there dimension and number
3. Prints a list of the boundary ids and the number of surface points for each
4. Plots the surface points of the specified boundary
"""

import sys
from scipy.io import netcdf_file
import numpy
import pdb

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

# =================================
# LOAD NETCDF FILE
# =================================
total_args = len(sys.argv)
if (total_args == 1):
	print ''
	sys.exit('Usage: python ncdump.py [netcdf filename] \n')

mesh = str(sys.argv[1])
print ''
print 'Found grid: ', mesh, '\n'

# Read file
cdf_file = netcdf_file(mesh, 'r')

# =================================
# OUTPUT DIMENSIONS AND VARIABLES
# =================================

# Print dimensions
print '---------Dimensions--------'
for entry in cdf_file.dimensions:
	print entry, ': ', cdf_file.dimensions[entry]
print ''

# Print variables
print '---------Variables---------'
for entry in cdf_file.variables:
	# Get number of entries
	variable_num = []
	for line in cdf_file.variables[entry]:
		variable_num.append(line)
	# Get length of entry type
	if isinstance(cdf_file.variables[entry][0], numpy.ndarray) == True:
		entry_leng = len(cdf_file.variables[entry][0])
	else:
		entry_leng = 1
	# Print
	print entry, ' ',entry_leng,' ',len(variable_num)	
print ''


# =================================
# GET SURFACE MARKER INFORMATION
# =================================
surface_marker_types =[]
surface_marker = []
key = 'boundary_ids'
key = 'boundarymarker_of_surfaces'

for entry in cdf_file.variables[key]:
	surface_marker.append(entry)
	if (entry in surface_marker_types) is False:
		surface_marker_types.append(entry)

print '-------Boundary IDs--------'
for entry in surface_marker_types:
	print entry,': ',surface_marker.count(entry)
print ''

# =================================
# GET GLOBAL IDS
# =================================
bid = raw_input('Boundary ID: ')
bid = int(bid)

gsid_list = []
for gsid, sm in enumerate(surface_marker):
	if sm == bid:
		gsid_list.append(gsid)

#gidfile = open('gidlist.txt', 'w')
#gidfile.write('Boundary ID: {}\n'.format(bid))
#for gid in gid_list:
#	gid = str(gid) + '\n'	
#	gidfile.write(gid)
#gidfile.close()


# =================================
# GET COORDINATE NAMES
# =================================
coord_entry_name = []
for entry in cdf_file.variables:
	if entry.startswith('x_'):
		coord_entry_name.append(entry)
for entry in cdf_file.variables:
	if entry.startswith('y_'):
		coord_entry_name.append(entry)
for entry in cdf_file.variables:
	if entry.startswith('z_'):
		coord_entry_name.append(entry)

print coord_entry_name

print ''

# =================================
# GET SURFACE MESH POINT IDs
# =================================
#boundary_choice = input('Choose boundary id: ')
#while boundary_choice not in surface_marker_types:
#	boundary_choice = input('Incorrect. Choose boundary id to plot: ')

#print ''

x_points = []
y_points = []
z_points = []

key = 'points_xc'
for entry in cdf_file.variables[key]:
        x_points.append(entry)

key = 'points_yc'
for entry in cdf_file.variables[key]:
        y_points.append(entry)

key = 'points_zc'
for entry in cdf_file.variables[key]:
        z_points.append(entry)



#sampled_x = []
#sampled_y = []
#sampled_z = []

#for val in range(len(surface_marker):
#	if surface_marker[val] == boundary_choice:
#		sampled_x.append(cdf_file.variables[coord_entry_name[0]][val])
#		sampled_y.append(cdf_file.variables[coord_entry_name[1]][val])
#		sampled_z.append(cdf_file.variables[coord_entry_name[2]][val])

# ================================
# GET Global Node IDs FROM SURFACE WITH SPECIFIED BID
# ================================
quadrilateral_points = []
triangle_points = []
surface_points = []
airfoil_points = []

key = 'points_of_surfacetriangles'
for entry in cdf_file.variables[key]:
        triangle_points.append(entry)

key = 'points_of_surfacequadrilaterals'
for entry in cdf_file.variables[key]:
        quadrilateral_points.append(entry)

surface_points = triangle_points + quadrilateral_points

#lnth = len(surface_points)

for i in range(len(gsid_list)):
        a = gsid_list[i]
	airfoil_points.append(surface_points[a])

#lnth = len(airfoil_points)
#pdb.set_trace()

NodeIDs = []
for i in range(len(airfoil_points)):
	a = airfoil_points[i]
	for j in range(len(a)):
		NodeIDs.append(a[j])

#pdb.set_trace()


# ================================
# GET Points and their Global Node IDs written in file
# ================================
xps = []
yps = []
zps = []

for i in range(len(NodeIDs)):
        a = NodeIDs[i]
        xps.append(x_points[a])
	yps.append(y_points[a])
	zps.append(z_points[a])

# Write file with info:
# x_coord y_coord z_coord GID
BPFile1 = open('BoundaryPointsx.txt', 'w')
BPFile2 = open('BoundaryPointsz.txt','w')
#BPFile.write('Boundary ID: {}\n'.format(bid))
#BPFile.write('x-coordinate     y-coordinate     z-coordinate     NodeID {}\n')
for i in range(len(NodeIDs)):
        #points = str(xps[i]) + ' ' + str(yps[i]) + ' ' + str(zps[i]) + ' ' + str(NodeIDs[i]) + '\n'
	points1 = str(xps[i]) + '\n'
	points2 =  str(zps[i]) + '\n'
        BPFile1.write(points1)
	BPFile2.write(points2)
BPFile1.close()
BPFile2.close()
#pdb.set_trace()

# =================================
# PLOT 3D AEROFOIL SECTIONS
# =================================
choice = raw_input('Plot aerofoil sections? (y/n): ')
if choice == 'y':
	fig = plt.figure()
	ax = fig.add_subplot(111, projection = '3d')
	ax.scatter(sampled_x,sampled_y,sampled_z, c='r')
	plt.xlabel('x')
	plt.ylabel('z')
	plt.show()

print ''

