#=========================================================================
#=========================================================================
'''
author Benjamin Lee Hinchliffe
commented and modified by Gabriele Luigi Mura
commented and modified by Jed Hollom

thi script output the "vertex normal" for each mesh points
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
sys.path.append('../mods/')
from scipy.io import netcdf_file
import read_taumesh as readm
import read_taumesh_normal as readm_normal
import numpy as np
import cPickle as pkl
from math import *



#====================
#USER DEFINE FUNCTION
#====================

def SUMM(sequence, start=0):
    for value in sequence:
        start = start + value
    return start


def DOT(AA, BB):
  return SUMM(((A*B) for A, B in zip(AA, BB)), start=0)


def ROTATE(bdr,alpha):
	
	theta=radians(alpha)
	
	RM=np.matrix([[cos(theta), -sin(theta)],[sin(theta), cos(theta)]])
	#print RM
	#print RM[0,:]
	final=[]
	for i in range(len(bdr)):
		for j in range(3):
			
		  	if j==1:
		  		bdr[i,j]=bdr[i,j]
		  	else:	
				dummy=[]		
  				for k in range(2):
					
					A=np.ones(shape=(2,2))
				
					B=np.array([bdr[i,0],bdr[i,2]])
					C=np.array([RM[k,0],RM[k,1]])
				
					D=DOT(C,B)
					
					dummy.append(D)

		final.append(dummy)
	return np.array(final)
  


#===================================================================================
#READING & SETTING THE VARIABLES FROM THE RELATIVE VOLUME COEFFICIENTS (NETCDF FILES)
#===================================================================================

print ('\n')
print('[1] reading relative volume coefficients')

#Read in solution of 1st step
f = netcdf_file('../relvolcoef_netcdf','r')
eco = f.variables['Rel_volume_coeffs'][:]      #deluanay "e" coefficient
vert = f.variables['Simplex_vert_index'][:]    #vertices of the delaunay triangulation
po_ref = f.variables['Point_index'][:]         #point reference
f.close()


#================
#READING THE MESH
#================
#Automatic mesh path finding
import os
for file in os.listdir('../'):
	if file.endswith('.taumesh') or file.endswith('.grid'):
		mesh = '../'+file
chord=1.0


bds=[[3],[4],[5],[6]] #Defining the boundaries

'''
NLF0416 2D Wing
[3] ff
[4] trailing_edge
[6] lower_surface
[6] upper_surface
'''


read_in = readm.readTaumesh() #Read mesh in
data, mesh_points = read_in.read(mesh,bds) #Get data and additional inner points
'''
-data- are all the points along with xyz marker and ID on the boundaries -ff,up,down and te-
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
bdry_data = pkl.load(open('../bdry_data.p','rb'))

data_ff = bdry_data["Data_ff"]
data_te = bdry_data["Data_te"]
data_ls = bdry_data["Data_ls"]
data_us=bdry_data["Data_us"]


sys.stdout.flush() #Force Printout

#=====================
#ROTATING SURFACE GRID 
#=====================
'''
initialisation of the new point, so the new point have the same coordinates execpt for the one we are moving. this make sense 
only if we are going to move each mesh point at the time
ech of this point has three column for x y z + ID but not sure about this one, it depends on how the data were extracted
'''

print('[4] rotating the surface')
#Shift and rotation angles
shift = 0.25
alpha=-10.e-0
print ('\n')
print ('rotating by an angle of'), alpha
print ('shifting by a distance of'), shift

#displacing at the 1/4 chord
#what we do is subtract 0.25 and add again however pay attention because we also need to redo this step when computing the differences

#Make copy of original data
original_data_te = data_te.copy()
original_data_ls = data_ls.copy()
original_data_us =data_us.copy()

print 'original: '+str(original_data_te[0])

#Shift
data_te[:,0]= data_te [:,0]-shift
data_ls[:,0] = data_ls[:,0]-shift
data_us[:,0]=data_us[:,0]-shift

#Make copy of shifted data
shift_data_te = data_te.copy()
shift_data_ls = data_ls.copy()
shift_data_us =data_us.copy()

print 'shift 1: '+str(shift_data_te[0])

#Rotate
rotate_te=ROTATE(data_te,alpha)
rotate_ls=ROTATE(data_ls,alpha)
rotate_us=ROTATE(data_us,alpha)

# #Make copy of rotated data
# rotate_data_te = data_te.copy()
# rotate_data_ls = data_ls.copy()
# rotate_data_us =data_us.copy()

print 'rotate: '+str(rotate_te[0])

#Rotate and Shift
data_te[:,0]=rotate_te[:,0]+shift
data_ls[:,0]=rotate_ls[:,0]+shift
data_us[:,0]=rotate_us[:,0]+shift

print ''
print data_te[0,2]
print rotate_te[0,1]
print ''

# Second Rotation
data_te[:,2]=rotate_te[:,1]
data_ls[:,2]=rotate_ls[:,1]
data_us[:,2]=rotate_us[:,1]


sys.stdout.flush() #Force Printout

#=====================
#COMPUTING SENSITIVITY
#=====================
'''
now we have made a manual deformation usign epsilon value
as we are using the eco (e coefficients), we need to use them to find the new point
this is done by using a inner point product between the e coeffiecient and the new coordinate of the new points
'''


print('[5] computing surface/rotation sensitivities')
orig_del_points = np.row_stack([data_ff[:,:4],original_data_te[:,:4],original_data_ls[:,:4],original_data_us[:,:4]])
def_del_points = np.row_stack([data_ff[:,:4],data_te[:,:4],data_ls[:,:4],data_us[:,:4]])

sys.stdout.flush() #Force Printout

#======================================================
#PROPAGATING THE DEFORMATION USING THE DELAUNAY MAPPING
#======================================================

'''
initialising the new vector that will containt all the mesh internal points
'''
def_mesh_points = mesh_points.copy()

'''
do that along x,y and z using user defined function DOT as we are working with small number it is preferred to work without using python built-in functions
'''
print('')
print('[6] propagating the deformation')
for i in xrange(len(po_ref)):

    print i,(' - '),len(po_ref)
    def_mesh_points[po_ref[i],0] = DOT(eco[i],def_del_points[vert[i],0])     
    def_mesh_points[po_ref[i],1] = DOT(eco[i],def_del_points[vert[i],1])
    def_mesh_points[po_ref[i],2] = DOT(eco[i],def_del_points[vert[i],2])


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


print('')
print ('[8] saving')
read_in.write_ver2(mesh,mesh+'_def',def_mesh_points,po_ref)
bdry_data_def = {"Data_ff":data_ff,"Data_te":data_te,"Data_ls":data_ls,"Data_us":data_us}


# Save as pickle
pkl.dump(bdry_data_def,open("../bdry_data_def.p",'wb'))

# Save as text file
# bdry_data_def_txt = open('../scripts/bdry_data_def.txt', 'w')
# for key in bdry_data_def.iterkeys():
# 	for entry in range(len(bdry_data_def[key])):
# 		line = str(bdry_data_def[key][entry][0])+' '+str(bdry_data_def[key][entry][2])
# 		bdry_data_def_txt.write(line+'\n')
# bdry_data_def_txt.close()


sys.stdout.flush() #Force Printout

#===
#EOF
#===
