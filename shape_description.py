import numpy as np


def compute_isolines(fiedler_vector,nbins=100):
	vmin = np.min(fiedler_vector)
	vmax = np.max(fiedler_vector)
	isolines = np.zeros((len(fiedler_vector),))
	intervals = np.linspace(vmin,vmax,nbins)
	for i in range(0,len(intervals)-1,2):
		isolines[np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])] = i
	return isolines,intervals
	
def compute_longitudinal_description(fiedler_vector,coords_nodes,nbins=100):
	# coords_nodes is a list of tuples
	coords = np.zeros((len(coords_nodes),3),dtype=int)
	for i in range(len(coords)):
		coords[i,:] = coords_nodes[i]
	vmin = np.min(fiedler_vector)
	vmax = np.max(fiedler_vector)
	barycenters = np.zeros((nbins-1,3))
	intervals = np.linspace(vmin,vmax,nbins-1)
	for i in range(0,len(intervals)-1):
		barycenters[i,:] = np.mean(coords[np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])],axis=0)
	return barycenters,intervals, coords
	
	
def distance(p,points):
	return np.sqrt((points[:,0]-p[0])**2 + (points[:,1]-p[1])**2 + (points[:,2]-p[2])**2)
	
	
def compute_thickness_stats(fiedler_vector,coords,nbins=100):
	barycenters, intervals, coords = compute_longitudinal_description(fiedler_vector, coords, nbins)
	thickness = np.zeros((len(intervals)-1,2))
	for i in range(0,len(intervals)-1):
		indices = np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])
		thickness[i,0] = 3*np.std(distance(barycenters[i,:], coords[indices,:]))
		thickness[i,1] = np.max(distance(barycenters[i,:], coords[indices,:]))
	return thickness
	
	
def compute_thickness(fiedler_vector,coords,nbins=100):
	# compute cortical thickness by using on each slice the Fiedler distance of the resulting graph
	intervals = np.linspace(vmin,vmax,nbins-1)
	thickness = np.zeros((len(intervals)-1,3))
	for i in range(0,len(intervals)-1):
		barycenters[i,:] = np.mean(coords[np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])],axis=0)
		thickness[i,0] = 3*np.std(distance(barycenters[i,:], coords[indices,:]))
		thickness[i,1] = np.max(distance(barycenters[i,:], coords[indices,:]))
		# Convert a slice in a graph
		
	return thickness
	
def texture_mapping(fiedler_vector, texture, intervals):
	texture_mapped = np.zeros((len(fiedler_vector),))
	for i in range(0,len(intervals)-1):
		indices = np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])
		texture_mapped[indices] = texture[i]
	return texture_mapped
