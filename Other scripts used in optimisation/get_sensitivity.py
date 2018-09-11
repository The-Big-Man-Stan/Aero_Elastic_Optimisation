# Imports
def get_sensitivity():
	file_0 = open('log/log_file.volgrad.solver.stdout.0','r')
	lines = file_0.read()
	file_0.close()
	Line_index = lines.index('Sensitivity') # Index of number is 16 away
		
	#Find the first digit in the string
	index_of_number = Line_index + 16
	
	#Get sensitivity value
	if lines[index_of_number-1] == '-':
		sensitivity = lines[index_of_number-1:index_of_number+22]
	else:
		sensitivity = lines[index_of_number:index_of_number+22]

	sensitivity = float(sensitivity)
	return sensitivity
