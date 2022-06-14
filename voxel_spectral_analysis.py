import numpy as np
import networkx as nx
import time

def graph_to_coords(graph):
	n = len(graph.nodes)
	coords = []# np.zeros((n,),dtype=int)
	for i,node in enumerate(graph.nodes):
		coords.append(node)
		#coords[i,:] = node
	return coords
	

def image_to_points(mask,voxel_size=[1,1,1]):
	""" 
	Create 3D points from a mask of voxels
	-mask: 3D image
	-voxel_size: by default 1mm isotropic
	""" 
	indices = np.asanyarray(np.where(mask !=0))
	for i in range(len(indices)):
		indices[i] = voxel_size[i] * indices[i]
	return indices
	
def add_valid_edge(g,inds,mask,graph_type="topology"):
	"""
	add edges from current voxel among the 26 possible neighbors if belongs to mask
	-g: graph
	-inds: 3 indices of the current voxel
	-mask: 3D image
	"""
	x,y,z=inds
	weight = 1
	sigma = 1
	#print(graph_type)
	for i in range(-1,2):
		for j in range(-1,2):
			for k in range(-1,2):
				if mask[x+i,y+j,z+k]!=0:
					if graph_type=="geometry":
						weight = np.exp(-(i**2 + j**2 + k**2)/sigma**2)
						#print(weight)
					g.add_edge(inds,(x+i,y+j,z+k),weight=weight)
	

def image_to_graph(mask,graph_type="topology"):
	ix, iy, iz = image_to_points(mask)
	g = nx.Graph()
	for i in range(len(ix)):
		add_valid_edge(g, (ix[i],iy[i],iz[i]), mask,graph_type)
	return g
	
def is_neighbour(p,q):
	# Test 26 connexity for p and q: sup norm
	return np.max(np.abs(p-q))==1
		
	
def points_to_graph(points,graph_type="topology"):
	g = nx.Graph()
	weight = 1
	sigma = 1
	for i in range(len(points)):
		for j in range(len(points)):
			if is_neighbour(points[i,:],points[j,:]):
				if graph_type=="geometry":
					weight = np.exp(-np.sum((points[i,:]-points[j,:])**2)/sigma**2)
				g.add_edge(points[i,:],points[j,:],weight = weight)
	return g
	
	
def get_diameter_fiedler(graph,compute_diameter=False):
	fiedler_vector = nx.fiedler_vector(graph)
	i_min = np.argmin(fiedler_vector)
	i_max = np.argmax(fiedler_vector)
	coords = graph_to_coords(graph) # object oriented approach would be better :-)
	node_min = coords[i_min]
	node_max = coords[i_max]
	#print(node_max,node_min,len(fiedler_vector))
	#print(graph.graph)
	t0 = time.time()
	diameter_fiedler = nx.shortest_path(graph,node_min,node_max)
	print("Fiedler distance, comp. time = ",time.time()-t0)
	diameter = -1
	if compute_diameter:
		t0 = time.time()
		diameter = nx.diameter(graph)
		print("Diameter, comp. time = ",time.time()-t0)
	return diameter_fiedler, diameter, fiedler_vector
	

