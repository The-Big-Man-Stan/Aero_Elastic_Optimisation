# Imports
def get_drag():
	with open('log/sim.monitor.pval.dat','r') as f:
		lines = f.readlines()
	Last_line = lines[-1] # Drag is the final value in the line
	Last_line = Last_line[:-1]
	
	for i,e in reversed(list(enumerate(Last_line))):
		if e == ' ':
			index = i
			break
	drag = float(Last_line[index:])
	return(drag)
