# Optimisation run script for rae2822 Adjoint optimisation with BFGS update and Delaunay graph mapping
#from BFGS import *
from fprime import *
from f import *
import numpy as np
import time
import subprocess
#from scipy import optimize
from BFGS import *

start = time.clock()
start1 = time.time()
process = subprocess.call(['cp','Optimisation_scripts/rae2822.bmap','rae2822.bmap'])

# Initialise optimisation with airfoils lcst parameters
rae2822_lcst_parameters = [0.12919237,0.12013562,0.17774294,0.07507853,0.29833847,-0.02940065,0.39216305,0.04469565,0.26027172,0.17421293,0.19707076,0.20970594,-0.12846913,-0.13816236,-0.12472814,-0.21670836,0.0306742,-0.53413441,0.27668371,-0.45672158,0.11430043,-0.14970269,0.01428077,0.05581]

gradient_convergence = 1e-5
rae2822_lcst_parameters = np.asarray(rae2822_lcst_parameters)

# Run optimisation
xopt,gopt,hopt,gradcalls,iterations = fmin_bfgs(f,fprime,rae2822_lcst_parameters)
end = time.clock()
end1 = time.time()

# Write results to file
result_file = open('GradientResults.txt','w')
result_file.write('Optimised lcst parameters\n')
result_file.write(str(xopt))
result_file.write('\n')
result_file.write('Final gradient value\n')
result_file.write(str(gopt))
result_file.write('\n')
result_file.write('Final Hessian value\n')
result_file.write(str(hopt))
result_file.write('\n')
result_file.write('Number of gradcalls\n')
result_file.write(str(gradcalls))
result_file.write('\n')
result_file.write('Number of iterations\n')
result_file.write(str(iterations))
result_file.write('\n')
result_file.write('Optimisiation run time\n')
result_file.write(str(end-start))
result_file.write('\n')
result_file.write(str(end1-start1))
result_file.close()
