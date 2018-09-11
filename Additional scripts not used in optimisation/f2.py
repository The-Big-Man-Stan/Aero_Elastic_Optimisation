# Imports
import subprocess
import os
import time
import numpy as np
from get_sensitivity import *
from delaunay_deform import *
from get_drag import *

def f(lcst_parameters):
	# Write all lcst iterations to file
        with open('lcst_iterations.txt','a') as myfile:
		myfile.write('f call')
                myfile.write('\n')
                myfile.write(str(lcst_parameters))
                myfile.write('\n')
		myfile.write('\n')

	#
	# Flow solver part of optimisation
	#
	# Create the mesh file and update Volume ratio cofficients
	delaunay_deform(lcst_parameters,'rae2822.taumesh1')
	#process = subprocess.Popen(['python','Find_vol_ratios.py'])
	#process.wait()
	
	# Preprocess grid for flow solvera
	process = subprocess.Popen(['ptau3d.preprocessing','rae2822.para_primal'])
	exitcode = process.wait()

	# If grid is bad,send back to path search
        if exitcode != 0:
                return 'error'

	# Submit job to Tau and wait for solution to be found
	process = subprocess.Popen(['mpirun','-n','7','ptau3d.turb1eq','rae2822.para_primal','log/log_file','use_mpi'])
	exitcode = process.wait()

        # If simulation doesn't converge
        if exitcode != 0:
                return 'error'

	# Get value of drag from current iteration
	drag = get_drag()
	with open('drag.txt','a') as myfile:
		myfile.write(str(drag))
		myfile.write('\n')
	
	# Output dummy file so fprime can know if f needs to be called
	with open(str(lcst_parameters[1]),'a') as myfile:
		myfile.write('dummy file')

	return drag
