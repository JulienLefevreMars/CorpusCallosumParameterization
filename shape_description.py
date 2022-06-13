import numpy as np


def compute_isolines(fiedler_vector,nbins=100):
	vmin = np.min(fiedler_vector)
	vmax = np.max(fiedler_vector)
	isolines = np.zeros((len(fiedler_vector),))
	intervals = np.linspace(vmin,vmax,nbins)
	for i in range(0,len(intervals)-1,2):
		isolines[np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])] = i
	return isolines,intervals
	
def compute_longitudinal_description(fiedler_vector,coords,nbins=100):
	# coords is a list of tuples
	vmin = np.min(fiedler_vector)
	vmax = np.max(fiedler_vector)
	barycenters = np.zeros((nbins-1,3))
	intervals = np.linspace(vmin,vmax,nbins-1)
	for i in range(0,len(intervals)-1):
		barycenters[i,:] = np.mean(coords[np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])],axis=0)
	return barycenters,intervals
	
	
def distance(p,points):
	return np.sqrt((points[:,0]-p[0])**2 + (points[:,1]-p[1])**2 + (points[:,2]-p[2])**2)
	
	
def compute_thickness(fiedler_vector,coords,nbins=100):
	barycenters, intervals = compute_longitudinal_description(fiedler_vector, coords, nbins)
	thickness = np.zeros((len(intervals)-1,2))
	for i in range(0,len(intervals)-1):
		indices = np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])
		thickness[i,0] = 3*np.std(distance(barycenters[i,:], coords[indices,:]))
		thickness[i,1] = np.max(distance(barycenters[i,:], coords[indices,:]))
	return thickness
	
def texture_mapping(fiedler_vector, texture, intervals):
	texture_mapped = np.zeros((len(fiedler_vector),))
	for i in range(0,len(intervals)-1):
		indices = np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])
		texture_mapped[indices] = texture[i]
	return texture_mapped
