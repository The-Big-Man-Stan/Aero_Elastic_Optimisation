#=====================================================================													
#=====================================================================													
'''  													    
Finds the Area/Volume ratios and outputs them into a text file.                                                                    													
To use this as a deformation method a second part will be required to 													
interpret the area/volume ratios to displace the points.     													
The main output is a file which contains the point number 													
and the ratios which are associated with it. 													
@Author: Benjamin Lee Hinchliffe (The University of Sheffield)													
@Modified by Gabriele Luigi Mura (The University of Sheffield)													
@Modified by Michael Jed Hollom (The University of Sheffield)													
@Modified by Anthony Stannard (The University of Sheffield)													
'''													
#=====================================================================													
#=====================================================================													
													
													
													
#=================													
#IMPORTING MODULES													
#=================													
from scipy.spatial import Delaunay													
from scipy.io import netcdf_file													
import numpy as np													
import read_taumesh as readm #Within the mods folder													
import time													
from rel_vol_funcs import rel_funcs #Within the 'rel_vol_funcs.py' file													
import cPickle as pkl													
from TwoDfuns import *													
import pdb													
													
#====================													
#SET SCRIPT VARIABLES													
#====================													
funs = rel_funcs()													
start = time.clock() #Starting timing													
													
#=========													
#READ MESH													
#=========													
'''													
This is the format of the information stored in the data lists													
[x coordinate, z coordinate, y coordinates, NodeIDs]													
'''													
#Enter mesh manually													
mesh = 'rae2822.taumesh'													
no_of_boundaries = 5													
#Create instance of readTaumesh class													
read_in = readm.readTaumesh() 													
#Get points with their corresponding Node IDs and boundary information													
#Inner points are points in Node ID order													
data_inner0, data_inner1, data_ff, data_us, data_ls, inner_points = read_in.read2(mesh,no_of_boundaries)													
													
#Separate arrays into y=1 and y=0 arrays due to 3D nature of Tau and put in 2D for triangulation													
#Inner nodes 0													
data_inner0,notUsed,points_inner0,notUsed = separateSymmetryPlanes(data_inner0)													
													
#Farfield													
data_ff0,data_ff1,points_ff0,points_ff1 = separateSymmetryPlanes(data_ff)													
													
#Upper surface													
data_us0,data_us1,points_us0,points_us1 = separateSymmetryPlanes(data_us)													
													
#Lower surface													
data_ls0,data_ls1,points_ls0,points_ls1 = separateSymmetryPlanes(data_ls)													
#============================													
#COMPUTING THE E-COEFFICEINTS													
#============================													
print('Computing the delaunay decomposition')													
# Create points for Delaunay domain													
points0 = points_ls0 + points_us0 + points_ff0													
													
# Make points into numpy array													
points = np.asarray(points0)													
													
#Give delaunay the points information to produce the tetrahedral/triangle													
tri = Delaunay(points) 													
													
#Create matrix to store area ratio information in													
e_coeffs = []													
vert_mat = []													
													
# Get area ratios for both planes of the mesh													
for i in range(len(points_inner0)):   													
	#Returns area ratio for each mesh point												
	e_co,vert = funs.find_coeffs_scipy(tri,i,points,points_inner0[i])												
	e_coeffs.append(e_co)												
	vert_mat.append(vert)												
	#Node IDs corresponding to the e_coeffs are stored in data_inner0[:,3] and data_inner1[:,3]												
													
#KEY:													
# tri: Triangle data which includes coordinates, triangle ID etc													
# points: The Delaunay boundary points													
# points_inner: The mesh points which will also contain surface points.													
													
# Make lists into numpy arrays													
e_coeffs = np.array(e_coeffs)													
vert_mat = np.array(vert_mat)													
													
sys.stdout.flush() #Force Printout													
													
#=========================													
#SAVING IN NETCDF FORMAT 0													
#=========================													
print('[7] saving in netcdf format')													
file_out = netcdf_file('../relvolcoef_netcdf','w')													
													
file_out.createDimension('no_points',i+1)													
file_out.createDimension('Delaunay_vert',3)													
													
temp1 = file_out.createVariable('Point_index','i',('no_points',))													
temp2 = file_out.createVariable('Rel_volume_coeffs','d',('no_points','Delaunay_vert'))													
temp3 = file_out.createVariable('Simplex_vert_index','i',('no_points','Delaunay_vert'))													
													
temp1[:] = data_inner0[:,3] # NodeIDs for symmetry plane y=0													
temp2[:,:] = e_coeffs[:,:]													
temp3[:,:] = vert_mat[:,:]													
													
sys.stdout.flush() #Force Printout													
file_out.close()													
													
#====================													
#SAVING PICKLE FILE													
#====================													
print('[8] compressing the file')													
bdry_data = {"Data_ff":data_ff0,"Data_us":data_us0,"Data_ls":data_ls0}													
pkl.dump(bdry_data,open("../bdry_data.p",'wb'))													
													
file_out = open('../eco_coef.p','wb')													
Rel_data = {"Point_index":a,"Relative_volume_coefficients":e_coeffs,"Simplex_vert_index":vert_mat}													
pkl.dump(Rel_data,file_out)													
file_out.close()													
sys.stdout.flush() #Force Printout													
													
end = time.clock() #Ending timer													
print ('this script took seconds'), end-start													
													
#=========													
#RESOURCES													
#=========													
'''													
computes and stores centroid of tetrahedron (3D) or triangle (2D)													
http://www.globalspec.com/reference/52702/203279/4-8-the-centroid-of-a-tetrahedron													
'''													
													
													
#===													
#EOF													
#===													
												