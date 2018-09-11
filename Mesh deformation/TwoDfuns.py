import numpy as np

# This function separates the data for the Taumesh so the points stored in an array are either on the symmetry plane at y=0 or the symmetry plane at y=1
def separateSymmetryPlanes(data):
	temp0a = []
	temp0b = []
	temp1a = []
	temp1b = []
	for i in range(len(data)):
		a = data[i][1]
		if (a == 0):
			temp0a.append([data[i][0],data[i][2]])
			temp0b.append([data[i][0],data[i][1],data[i][2],data[i][3]])
		else:
			temp1a.append([data[i][0],data[i][2]])
			temp1b.append([data[i][0],data[i][1],data[i][2],data[i][3]])
	
	# Remove superflous 1st line from arrays
	data0 = np.asarray(temp0b)
	data1 = np.asarray(temp1b)
	points0 = temp0a
	points1 = temp1a

	return data0,data1,points0,points1

def separateMeshPoints(data):
	temp0a = []
	temp0b = []
	temp1a = []
	temp1b = []
	NodeIDs = range(len(data))
	print 'NodeID length',' ',len(data)
	print len(data)
	for i in range(len(data)):
		a = data[i][1]
		if (a == 0):
			temp0a.append([data[i][0],data[i][2]])
			temp0b.append([data[i][0],data[i][1],data[i][2],NodeIDs[i]])
		else:
			temp1a.append([data[i][0],data[i][2]])
			temp1b.append([data[i][0],data[i][1],data[i][2],NodeIDs[i]])
	
	# Remove superflous 1st line from arrays
	data0 = np.asarray(temp0b)
	data1 = np.asarray(temp1b)
	points0 = temp0a
	points1 = temp1a

	return data0,data1,points0,points1
