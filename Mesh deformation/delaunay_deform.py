#=========================================================================
#=========================================================================
'''
author Benjamin Lee Hinchliffe
commented and modified by Gabriele Luigi Mura
commented and modified by Jed Hollom
commented and modified by Anthony Stannard

this script outputs the "vertex normal" for each mesh points
although normal are defined with respect to a face the definition can 
be forced and it can be possible to assigne a normal to a vertex
see TAU user guide
'''
#=========================================================================
#=========================================================================




#=================
#IMPORTING MODULES
#=================
import sys
from scipy.io import *
import read_taumesh as readm
import numpy as np
import cPickle as pkl
import lcst
from math import *
from TwoDfuns import *
import time



#====================
#USER DEFINE FUNCTION
#====================

def SUMM(sequence, start=0):
    for value in sequence:
        start = start + value
    return start


def DOT(AA, BB):
  return SUMM(((A*B) for A, B in zip(AA, BB)), start=0)

def delaunay_deform(new_a,name):
	#===================================================================================
	#READING & SETTING THE VARIABLES FROM THE RELATIVE VOLUME COEFFICIENTS (NETCDF FILES)
	#===================================================================================
	
	print ('\n')
	print('[1] reading relative volume coefficients')
	start = time.clock()
	#Read in solution of 1st step
	f = netcdf_file('relvolcoef_netcdf','r')
	eco = f.variables['Rel_volume_coeffs'][:]      #deluanay "e" coefficient
	vert = f.variables['Simplex_vert_index'][:]    #vertices of the delaunay triangulation
	po_ref = f.variables['Point_index'][:]         #point reference
	f.close()
	
	#================
	#READING THE MESH
	#================
	#Automatic mesh path finding
	#import os
	#for file in os.listdir('../'):
	#	if file.endswith('.taumesh') or file.endswith('.grid'):
	#		mesh = '../'+file
	#chord=1.0
	
	
	#bds=[[3],[4],[5],[6]] #Defining the boundaries
	
	'''
	RAE2822 2D Wing
	[3] ff
	[4] upper_surface
	[5] lower_surface
	'''
	
	mesh = 'rae2822.taumesh'
	read_in = readm.readTaumesh() #Read mesh in
	no_of_bounds = 5
	data_inner0, data_inner1, data_ff, data_us, data_ls, mesh_points = read_in.read2(mesh,no_of_bounds)
	mesh_points0,mesh_points1,na,na = separateMeshPoints(mesh_points) # Only need mesh points on y=0 plane
	
	# Store Node IDs for y = 0 plane and y = 1 plane
	NodeIDs0 = data_inner0[:,3]
	NodeIDs1 = data_inner1[:,3]
	'''
	-data- are all the points along with xyz marker and ID on the boundaries -ff,up,down and
	-mesh_points- these are all the points including the volume internal mesh point
	'''
	
	
	sys.stdout.flush() #Force Printout
	
	#==================================
	#READING FROM PICKLE THE BOUNDARIES
	#==================================
	'''
	here we are reading from a pickle file, althogh not necessary for such a small mesh it will save space later on
	when using large mesh
	if you run the command 
	print ('to compare'), data_up[0]
	you will see that tha data are the same, what changes is only the format of the file
	'''
	
	#Load the pickle file (the mesh file but in a different format)
	print('[3] reading pickle file')
	bdry_data = pkl.load(open('bdry_data.p','rb'))
	
	data_ff = bdry_data["Data_ff"]
	data_ls = bdry_data["Data_ls"]
	data_us = bdry_data["Data_us"]
	
	
	sys.stdout.flush() #Force Printout
	
	#=====================
	#ROTATING SURFACE GRID 
	#=====================
	'''
	initialisation of the new point, so the new point have the same coordinates execpt for the one we are moving. this make sense 
	only if we are going to move each mesh point at the time
	ech of this point has three column for x y z + ID but not sure about this one, it depends on how the data were extracted
	'''
	
	print('[4] deforming the surface')
	
	#Make copy of original data
	#Upper surface
	original_data_us0,original_data_us1,na,na = separateSymmetryPlanes(data_us)
	#Lower surface
	original_data_ls0,original_data_ls1,na,na = separateSymmetryPlanes(data_ls)
	
	# Original a values
	# au: [ 0.12919237  0.12013562  0.17774294  0.07507853  0.29833847 -0.02940065 0.39216305  0.04469565  0.26027172  0.17421293  0.19707076  0.20970594]
	# al: [-0.12846913 -0.13816236 -0.12472814 -0.21670836  0.0306742  -0.53413441  0.27668371 -0.45672158  0.11430043 -0.14970269  0.01428077  0.05581   ]
	
	# Get new coordinates of airfoil
	# Create instance of lcst class
	LCST = lcst.lcst()
	upperSurf,lowerSurf = LCST.get_profile(new_a,data_us[:,0],data_ls[:,0])
	
	# New coordinates of airfoil after change in the parameters
	data_ls[:,2] = lowerSurf[1]
	data_us[:,2] = upperSurf[1]
	
	sys.stdout.flush() #Force Printout
	
	#=====================
	#COMPUTING SENSITIVITY
	#=====================
	'''
	now we have made a manual deformation sign epsilon value
	as we are using the eco (e coefficients), we need to use them to find the new point
	this is done by using a inner point product between the e coeffiecient and the new coordinate of the new points
	'''
	
	
	print('[5] computing surface sensitivities')
	#Farfield
	data_ff0,data_ff1,na,na = separateSymmetryPlanes(data_ff)
	
	#Upper surface
	data_us0,data_us1,na,na = separateSymmetryPlanes(data_us)
	
	#Lower surface
	data_ls0,data_ls1,na,na = separateSymmetryPlanes(data_ls)
	
	#Orig_points0 must be in the same order points0 in Find Vol ratios script for vertices info to be accurate
	orig_points0 = np.row_stack([original_data_ls0[:,:4],original_data_us0[:,:4],data_ff0[:,:4]])
	def_points0 = np.row_stack([data_ls0[:,:4],data_us0[:,:4],data_ff0[:,:4]])
	#orig_points1 = np.row_stack([original_data_ls1[:,:4],original_data_us1[:,:4],data_ff1[:,:4]])
	#def_points1 = np.row_stack([data_ls1[:,:4],data_us1[:,:4],data_ff1[:,:4]])
	
	sys.stdout.flush() #Force Printout
	
	#======================================================
	#PROPAGATING THE DEFORMATION USING THE DELAUNAY MAPPING
	#======================================================
	
	'''
	initialising the new vector that will containt all the mesh internal points
	'''
	def_mesh_points = mesh_points0.copy() # In the correct node order
	
	'''
	do that along x,y and z using user defined function DOT as we are working with small number it is preferred to work without using python built-in functions
	'''
	print('')
	print('[6] propagating the deformation')
	for i in po_ref:
	
	    # i is the Node ID
	    def_mesh_points[i,0] = DOT(eco[i],def_points0[vert[i],0])     
	    def_mesh_points[i,1] = DOT(eco[i],def_points0[vert[i],1])
	    def_mesh_points[i,2] = DOT(eco[i],def_points0[vert[i],2])
	
	
	sys.stdout.flush() #Force Printout
	
	#=============================
	#PRINTING STATISTICS ON SCREEN
	#=============================
	
	print('')
	print ('[7] check on the original mesh')
	counter=0
	
	for j in range(3):
		for i in range(len(def_mesh_points)):
			a=mesh_points[i,j] 
			b=def_mesh_points[i,j]
		
			if a!=b and abs(a-b)>1e-12:	
				#print a, b
				#print ('difference'), abs(a-b)
				counter +=1
	print ('point different'), counter/3
	print ('point in total'), len(def_mesh_points)
	
	print('')
	
	
	
	
	sys.stdout.flush() #Force Printout
	
	#========================================================================
	#SAVING THE UPDATED PICKLE FILE & THE NEW MESH FILE CONTAINING DEFORMATION
	#========================================================================
	# Add the y=1 plane mesh points
	def_mesh_points1 = def_mesh_points.copy()
	def_mesh_points1[:,1] = 1
	def_mesh_points1[:,3] = mesh_points1[:,3]
	def_mesh_points = np.row_stack([def_mesh_points,def_mesh_points1])
	
	# Get boundary ID info for every mesh point
	NodeIDs1 = data_inner0[:,3]
	NodeIDs2 = data_inner1[:,3]
	NodeIDsff = data_ff[:,3]
	NodeIDsus = data_us[:,3]
	NodeIDsls = data_ls[:,3]
	
	print('')
	print ('[8] saving')
	read_in.write_tau2D_2(mesh,name,def_mesh_points)
	#read_in.write_tau2D_4(mesh+'_def3',def_mesh_points)
	
	sys.stdout.flush() #Force Printout
	end = time.clock()
	print ('this script took seconds'), end-start
	#===
	#EOF
	#===
