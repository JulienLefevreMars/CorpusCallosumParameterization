import numpy as np
import voxel_spectral_analysis as vsa
import sklearn.decomposition as sd
import matplotlib.pyplot as plt

def compute_isolines(fiedler_vector,nbins=100):
	vmin = np.min(fiedler_vector)
	vmax = np.max(fiedler_vector)
	isolines = np.zeros((len(fiedler_vector),))
	intervals = np.linspace(vmin,vmax,nbins)
	for i in range(0,len(intervals)-1,2):
		isolines[np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<intervals[i+1])] = i
	return isolines,intervals

def irregular_binning(fiedler,nbins=100):
	vmin = np.min(fiedler)
	vmax = np.max(fiedler)
	# Same as histogram equalization
	hist,bin_edges = np.histogram(fiedler,bins=nbins)
	# len(bin_edges) := nbins + 1
	cdf = np.cumsum(hist)/np.sum(hist)
	plt.plot(bin_edges[0:-1],cdf)
	plt.title("Cdf of Fiedler")
	plt.show()
	new_fiedler = np.zeros(fiedler.shape)
	print(cdf)
	for i in range(len(bin_edges)-1):
		indices = np.logical_and(fiedler>=bin_edges[i],fiedler<=bin_edges[i+1])
		new_fiedler[indices] = (vmax-vmin) * cdf[i] + vmin
	plt.hist(new_fiedler,bins=int(nbins/10))
	plt.title("Histogram of reparameterized Fiedler")
	plt.show()
	intervals = np.linspace(np.min(new_fiedler),np.max(new_fiedler),nbins)
	return intervals, new_fiedler
	
def compute_longitudinal_description(fiedler_vector,coords_nodes,nbins=100,irregular_bins=False):
	'''
	Provides barycenters/skeleton of a shape by using isolines of Fiedler vector
	'''
	coords = np.zeros((len(coords_nodes),3),dtype=int)
	for i in range(len(coords)):
		coords[i,:] = coords_nodes[i]
	vmin = np.min(fiedler_vector)
	vmax = np.max(fiedler_vector)
	print("Min, max Fiedler",vmin,vmax)
	barycenters = np.zeros((nbins-1,3))
	#intervals = np.linspace(vmin,vmax,nbins-1)
	if irregular_bins:
		bins, new_fiedler_vector = irregular_binning(fiedler_vector,nbins)
	else:
		bins = np.linspace(vmin,vmax,nbins)
		new_fiedler_vector = fiedler_vector
	#print(fiedler_vector)
	print(len(bins),nbins)
	for i in range(0,nbins-1):
		barycenters[i,:] = np.mean(coords[np.logical_and(new_fiedler_vector >=bins[i],new_fiedler_vector<bins[i+1])],axis=0)
	return barycenters,bins, coords, new_fiedler_vector
	
	
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
	barycenters, bins, coords, fiedler_vector = compute_longitudinal_description(fiedler_vector, coords, nbins)
	print(bins)
	print("Min, max Fiedler",np.min(fiedler_vector),np.max(fiedler_vector))
	thickness = np.zeros((len(bins)-1,4))
	slices = []
	for i in range(0,len(bins)-1):
		print(i,bins[i],bins[i+1], np.sum(np.logical_and(fiedler_vector >=bins[i],fiedler_vector<=bins[i+1])))

	
	for i in range(0,len(bins)-1):
		print(i)
		# 1. Skeletonization
		indices = np.logical_and(fiedler_vector >=bins[i],fiedler_vector<=bins[i+1])
		# warning, if <bins[i+1] indices is empty for i = len(bins)-1 because of histogram equalization that puts a lot of values on vmax
		slice_points = coords[indices]
		barycenters[i,:] = np.mean(slice_points,axis=0)
		thickness[i,0] = 3*np.std(distance(barycenters[i,:], coords[indices,:]))
		thickness[i,1] = np.max(distance(barycenters[i,:], coords[indices,:]))
		# 2. Convert a slice in a graph
		slice_graph = vsa.points_to_graph(slice_points,graph_type="geometry")
		#print(slice_graph)
		res = vsa.get_diameter_fiedler(slice_graph)
		#print(res)
		thickness[i,2] = res[0]
		slices.append([slice_graph,res[2]])
		# 3. PCA
		pca = sd.PCA(n_components=2)
		pca.fit(slice_points)
		projection = pca.transform(slice_points)
		thickness[i,3] = np.max(projection) - np.min(projection)
		if i==9:
			plt.plot(coords[:,1],coords[:,2],'+r')
			plt.plot(slice_points[:,1],slice_points[:,2],'+')
			plt.show()
			
		print(thickness[i,:])
	return thickness, slices, bins
	
	
def texture_mapping(fiedler_vector, texture, intervals):
	texture_mapped = np.zeros((len(fiedler_vector),))
	for i in range(0,len(texture)):
		indices = np.logical_and(fiedler_vector >=intervals[i],fiedler_vector<=intervals[i+1])
		texture_mapped[indices] = texture[i]
	return texture_mapped
