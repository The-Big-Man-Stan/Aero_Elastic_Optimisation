# Imports
import subprocess
import sys
import os
import time
import numpy as np
from get_sensitivity import *
from delaunay_deform import *
from get_drag import *
from f import *
import copy

def fprime(lcst_parameters):
	# Call the flow solver if needed
	process = subprocess.Popen(['rm',str(lcst_parameters[1])])
	exitcode = process.wait()
	if exitcode != 0: 
		drag = f(lcst_parameters)
	
	#
	# Residual part of optimisation
	#
	# Create parameter file for residual run
	process = subprocess.Popen(['cp','rae2822.para_primal','rae2822.para_residual'])
	process.wait()
	
	# Write necessary lines on bottom of residual parameter file
	appended_text = 'Solver type: Residual_only'
	with open('rae2822.para_residual','a') as myfile:
		myfile.write('\n')
		myfile.write(appended_text)
	
	# Submit residual job to Tau
	process = subprocess.Popen(['mpirun','-n','7','ptau3d.turb1eq','rae2822.para_residual','log/log_file.residual','use_mpi'])
	process.wait()
	
	#
	# Adjoint part of optimisation
	#
	# Create adjoint parameter file
	process = subprocess.Popen(['cp','rae2822.para_residual','rae2822.para_adjoint'])
	process.wait()
	
	# Write necessary lines on bottom of adjoint parameter file
	with open('rae2822.para_adjoint','a') as myfile:
		myfile.write('\n')
		myfile.write('Solver type: DAdjoint')
		myfile.write('\n')
		myfile.write('Solve dissipation error equation (0/1): 0')
		myfile.write('\n')
		myfile.write('Solve linear problem on grid level: 1')
		myfile.write('\n')
		myfile.write('Cost function part total/pressure/viscous (0/1/2): 0')
		myfile.write('\n')
		myfile.write('Point of point pressure cost function: 0')
		myfile.write('\n')
		myfile.write('Cost function: C-drag')
		myfile.write('\n')
		myfile.write('Design variable: Farfield-alpha')
		myfile.write('\n')
		myfile.write('Monitoring values: Residual_dJ/da-1_dJ/da-2_dJ/da')
		myfile.write('\n')
		myfile.write(' Monitoring significant figures: 4_8_8_8')
		myfile.write('\n')
		myfile.write('Jacobian variables: Cons')
		myfile.write('\n')
		myfile.write('Jacobian constant laminar viscosity (0/1): 0')
		myfile.write('\n')
		myfile.write('Jacobian frozen turbulence (0/1): 0')
		myfile.write('\n')
		myfile.write('Jacobian constant dissipation coeffs (0/1): 0')
		myfile.write('\n')
		myfile.write('Krylov loop: GMRes')
		myfile.write('\n')
		myfile.write('GMRes inner iterations: 40')
		myfile.write('\n')
		myfile.write('GMRes preconditioning iterations: 40')
		myfile.write('\n')
		myfile.write('Minimum residual: 1e-4')
		myfile.write('\n')
		myfile.write('Preconditioning: (none)')
	
	# Submit adjoint job
	process  = subprocess.Popen(['mpirun','-n','7','ptau3d.turb1eq','rae2822.para_adjoint','log/adjoint.C-func','use_mpi'])
	process.wait()
	
	#
	# Calculating the gradient section
	#
	# Create the Volgrad parameter file to get the sensitivity
	process = subprocess.Popen(['cp','rae2822.para_adjoint','rae2822.para_volgrad'])
	process.wait()
	
	with open('rae2822.para_volgrad','a') as myfile:
		myfile.write('\n')
		myfile.write('Primary grid filename: rae2822.taumesh_def')
		myfile.write('\n')
		myfile.write('Grid prefix: grid2/')
		myfile.write('\n')
		myfile.write('Solver type: Volgrad')
		myfile.write('\n')
		myfile.write('Automatic parameter update (0/1): 0')
	
	# Create del to use for each parameter when doing finite difference
	delta = 0.00000001*np.ones(len(lcst_parameters))
	gradient = range(len(delta))
	
	# Create deformed meshes for each parameter
	for i in range(len(delta)):
		print('Volgrad iteration: ',str(i))
		temp = copy.deepcopy(lcst_parameters)
		temp[i] = lcst_parameters[i] + delta[i]
		delaunay_deform(temp,'rae2822.taumesh_def') # This outputs the deformed mesh
		process = subprocess.Popen(['ptau3d.preprocessing','rae2822.para_volgrad'])
		process.wait()
		process = subprocess.Popen(['mpirun','-n','7','ptau3d.turb1eq','rae2822.para_volgrad','log/log_file.volgrad','use_mpi'])
		process.wait()
		sensitivity = get_sensitivity()
		gradient[i] = sensitivity/delta[i]
	
	# Clear up parameter files for next gradient call
	process = subprocess.call(['rm','rae2822.para_volgrad'])
	process = subprocess.call(['rm','rae2822.para_adjoint'])
	process = subprocess.call(['rm','rae2822.para_residual'])
	gradient = np.asarray(gradient)
	
	# Write all gradients to file
	with open('gradients.txt','a') as myfile:
		myfile.write('\n')
		myfile.write(str(gradient))
		myfile.write('\n')
	gradient = 0.01 * gradient
	return gradient
