# Imports
import subprocess
import os
def get_alpha_sensitivities():
	# Use grep to get the Farfield-alpha sensitivities
	command = 'grep Farfield-alpha log/adjoint.C-func.solver.stdout.0 >farfield_alpha_line'
	os.system(command)
	#process = subprocess.Popen(['grep','Farfield-alpha','log/adjoint.C-func.solver.stdout.0','>farfield_alpha_lines'])
	
	# There are 2 instances of Fafield alpha in the log file and we are interested in the second
	# instance
	with open('farfield_alpha_lines','r') as f:
	        lines = f.readlines()
		
	
	# Want to extract the first 2 numbers in this line. The first is the partial lift
	# with respect to the angle of attack. The second is the adjoint field multiplied
	# by the residual partially differentiated with respect to the angle of attack
	line = lines[1]
	count = 0
	
	for i,e in list(enumerate(line)):
		if e.isdigit() and count == 0:
			if line[i-1] == '-':
				index1 = i-1
				count = 1
			else:
				index1 = i
				count = 1
		if count == 1 and e == ' ':
			#index2 = i-1 
			lift_sensitivity = float(line[index1:i])
			count = 2
		if count == 2 and e.isdigit():
			if line[i-1] == '-':
				index3 = i-1
				count = 3
			else:
				index3 = i
				count = 3
		if count == 3 and e == ' ':
			#index4 = i-1
			residual_sensitivity = float(line[index3:i])
			break
	
	return([lift_sensitivity,residual_sensitivity])
[a,b] = get_alpha_sensitivities()
print(a)
print(b)
