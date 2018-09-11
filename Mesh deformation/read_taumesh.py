'''
Created on 23 Jul 2011

@author: fzhu
Modified by Anthony Stannard 2018
'''

from scipy.io import netcdf_file
import numpy,os,shutil
import netCDF4
from scipy.io import *

class readTaumesh():
    def read2(self,mesh,no_of_boundaries):
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

        return data_inner0, data_inner1, data_ff, data_us, data_ls, numpy.column_stack((points_x,points_y,points_z))

    def write_tau2D_2(self, file_in, file_out, data):
	
	if file_in != file_out:
		shutil.copy(file_in,file_out)
	
	mesh_file = netCDF4.Dataset(file_out,'r+')
	mesh_file.variables['points_xc'][:] = data[:,0]		
	mesh_file.variables['points_yc'][:] = data[:,1]
	mesh_file.variables['points_zc'][:] = data[:,2]
	mesh_file.close()
	return
	
    def write_tau2D_3(self, file_in, file_out, data):
	mesh_file = netCDF4.Dataset(file_out,'w','NETCDF3_CLASSIC')
	# Set the dimension of the netCDF file for the variables
	mesh_file.createDimension('no_of_points',len(data[:,0]))

	# Create the variables to be used for finite differences in TAU
	X = mesh_file.createVariable('x','f4','no_of_points')
	Y = mesh_file.createVariable('y','f4','no_of_points')
	Z = mesh_file.createVariable('z','f4','no_of_points')
	Global_ID = mesh_file.createVariable('global_id','f4','no_of_points')

	# Set variables to be used for finite differences in TAU
	X[:] = data[:,0]		
	Y[:] = data[:,1]
	Z[:] = data[:,2]
	Global_ID[:] = data[:,3]
	mesh_file.close()
	return

    def write_tau2D_4(self,file_out,data):
	BPFile = open('DEFORM.ASCII','w')
	for i in range(len(data[:,3])):
		points = str(data[i,0]) + ' ' + str(data[i,1]) + ' ' + str(data[i,2]) + ' ' + str(data[i,3]) + '\n'
		BPFile.write(points)
	BPFile.close()
	return
	# Version one is in the mods file as writeDef mesh, its to be used if netCDF4 isn't available
    # ----------------------
    
    
    
    






