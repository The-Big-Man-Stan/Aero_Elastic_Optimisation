#=====================================================================
#=====================================================================
'''
this script defines a class called "rel_funcs"					     	
inside which there are more than 10 functions which are called in by 
the main script
@Author: Benjamin Lee Hinchliffe (The University of Sheffield)
@Modified by Gabriele Luigi Mura (The University of Sheffield)
@Modified by Michael Jed Hollom (The University of Sheffield)
'''
#=====================================================================
#=====================================================================


#=================
#IMPORTING MODULES
#=================

import sys
from scipy.spatial import Delaunay
import scipy
from scipy import linalg
import numpy as np
from itertools import product, islice
import cPickle as pkl



#===============
#DEFINING CLASS
#===============


class rel_funcs():

	def find_boundaries(self,data_original,bds):
		'''
		THIS IS SPECIFIC FOR THE THE MESH USED: NLF0416
		'''
		#Get length of data
		num_tot = len(data_original[:,0])
  
  		#Define arrays to fill
		data_ff = []
		data_te = []
		data_ls = []
		data_us = []  		
	   
	   	#Fill arrays with data
		for i in range(num_tot):
			if int(data_original[i,4]) in bds[0]:
				data_ff.append(data_original[i,:])
			if int(data_original[i,4]) in bds[1]:
				data_te.append(data_original[i,:])
			if int(data_original[i,4]) in bds[2]:
				data_ls.append(data_original[i,:])
			if int(data_original[i,4]) in bds[3]:
				data_us.append(data_original[i,:])

		#Return seperated data arrays
		return np.array(data_ff),np.array(data_us),np.array(data_ls),np.array(data_te)


	'''
	this function was defined so the ocomputed determinant would a precision stored in float 128
	this is because we have noticed some float cancellation error using the old built in function in python    
	'''
	
	
	def DET(self,m):
				
		a0003=(m[0,0]-m[0,3])
		a0103=(m[0,1]-m[0,3])	
		a0203=(m[0,2]-m[0,3])
		
		a= a0003 * ( ((m[1,1]-m[1,3]) * (m[2,2]-m[2,3])) - ((m[2,1]-m[2,3])*(m[1,2]-m[1,3])) ) 
		b= a0103 * ( ((m[1,0]-m[1,3]) * (m[2,2]-m[2,3])) - ((m[2,0]-m[2,3])*(m[1,2]-m[1,3])) ) 	
		c= a0203 * ( ((m[1,0]-m[1,3]) * (m[2,1]-m[2,3])) - ((m[2,0]-m[2,3])*(m[1,1]-m[1,3])) )
		M=a-b+c

		return M
		




		#THE FOLLOWING USE THE SCIPY MODULE CAPABILITIES


	def find_coeffs_scipy(self,tri,j,points,inner_points):
		
		simp = tri.find_simplex(inner_points) # Find simplicity
		vert = tri.vertices[simp] #find vertices of triangle each point is in
		temp_mat1 = np.ones((3,3)) # Put 'vert' into a matrix
		temp_mat1[:,:2]=points[vert] # Gives matrix of Delaunay vertices
		V = np.linalg.det(temp_mat1) #Compute determinant using own method
		V_mat = points[vert] # Get matrix for first 2 columns to build sub-elements from
		e_co= []		
		

		if V==0.0:
			sys.exit('ERROR: Volume = 0 at iteration: '+str(j))
		else:
			e_co= []
			sub_VV=[]    
			for i in range(3):  	
			
				temp_mat = np.ones((3,3))
				temp_mat[:,:2] = V_mat
				temp_mat[i,:2] = inner_points
				
				temp_mat=temp_mat.T		
				sub_V=np.linalg.det(temp_mat)
				
				c=abs(sub_V/V)
				
				e_co.append(c)
			
			
		#sum_e=abs(e_co[0])+abs(e_co[1])+abs(e_co[2])
		return e_co, vert,# sum_e #E coeffieicient, verticies and sum on E
		



		

#=========
#RESOURCES
#=========

'''
see http://www.had2know.com/academics/tetrahedron-volume-4-vertices.html
'''



#===
#EOF
#===
