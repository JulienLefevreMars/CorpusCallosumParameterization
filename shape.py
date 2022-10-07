""""
Definition of a class to represent shapes as graphs + scalar field (Fiedler vectors) + description of this field


"""

# Author: Julien Lefevre, PhD, julien.lefevre@univ-amu.fr

import networkx as nx
import nibabel as nb
import numpy as np
import shape_description as sd
import time

def is_neighbour(p,q):
	# Test 26 connexity for p and q: sup norm
	return np.max(np.abs(p-q))==1
	

def points_to_graph(points,graph_type="topology"):
	g = nx.Graph()
	weight = 1 # TO DO: to be adapted
	sigma = 0.5
	for i in range(len(points)):
		for j in range(len(points)):
			if is_neighbour(points[i,:],points[j,:]):
				if graph_type=="geometry":
					weight = np.exp(-np.sum((points[i,:]-points[j,:])**2)/sigma**2)
					#print(weight)
				g.add_edge(tuple(points[i,:]),tuple(points[j,:]),weight = weight)
	return g

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
	
def image_to_graph(mask,graph_type="topology"):
	ix, iy, iz = image_to_points(mask)
	g = nx.Graph()
	for i in range(len(ix)):
		add_valid_edge(g, (ix[i],iy[i],iz[i]), mask,graph_type)
	return g
	
def add_valid_edge(g,inds,mask,graph_type="topology"):
	"""
	add edges from current voxel among the 26 possible neighbors if belongs to mask
	-g: graph
	-inds: 3 indices of the current voxel
	-mask: 3D image
	"""
	x,y,z=inds
	weight = 1 # TO DO: to be adapted
	sigma = 2.
	#print(graph_type)
	for i in range(-1,2):
		for j in range(-1,2):
			for k in range(-1,2):
				if mask[x+i,y+j,z+k]!=0:
					if graph_type=="geometry":
						weight = np.exp(-(i**2 + j**2 + k**2)/sigma**2)
						#print(weight)
					g.add_edge(inds,(x+i,y+j,z+k),weight=weight)


def threshold(points, axis=0, thres=0.5, inequality = "larger"):
	values = points[:,axis]
	m = np.min(values)
	L = np.max(values) - m
	if inequality == "larger":
		indices = np.where(values > thres*L + m)
	else:
		indices = np.where(values < thres*L + m)
	print(indices)
	return points[indices[0],:]


class Shape:
	def __init__(self,filename=None,graph=None,graph_type="geometry",**kwargs):
		if not(filename is None):
			g = nb.load(filename)
			mask = np.asanyarray(g.dataobj)
			self.coords_triplet = image_to_points(mask,g.header['pixdim'][1:4])
			self.graph = image_to_graph(mask,graph_type=graph_type)
			self.is_empty = False
		elif not(graph is None):
			self.graph = graph
			self.coords_triplet = None
		else:
			self.is_empty = True
		self.fiedler_vector = None
		self.description = sd.ShapeDescription(None) # given by class shape_description
		self.compute_numpy_coords()
	
	def graph_perturbation(self,coord,weight=0.1):
		# Perturbation of node ind to make as if it was "far" from the rest of the graph
		for n in self.graph.nodes:
			dist = sum([ (n[i] - coord[i])**2 for i in range(3) ])
			if dist==0:
				break
		print(n)
		print(coord)
		#for e in self.graph.edges([n]):
		#	self.graph[n][e[1]]['weight']=weight
		# Add a fictive edge
		fictive_node = (n[0],n[1]-10,n[2]-10)
		self.graph.add_edge(fictive_node,n,weight = weight)
		return fictive_node

	def get_fiedler_perturbation(self,coord,weight):
		fictive_node = self.graph_perturbation(coord,weight)
		self.get_fiedler()
		self.graph.remove_node(fictive_node)
		self.fiedler_vector = self.fiedler_vector[0:-1]
		
	
	def get_fiedler(self):
		self.fiedler_vector = nx.fiedler_vector(self.graph)
		i_min = np.argmin(self.fiedler_vector)
		i_max = np.argmax(self.fiedler_vector)
		coords = self.graph_to_coords()
		if coords[i_min][1]>coords[i_max][1]: # y coordinate is approximately aligne with Fiedler vector, but BE CAREFUL, could change depending on the way the MRI acquisition is done
			self.fiedler_vector = - self.fiedler_vector
		self.i_min = i_min
		self.i_max = i_max
		
	def graph_to_coords(self):
		n = len(self.graph.nodes)
		coords_nodes = []# np.zeros((n,),dtype=int)
		for i,node in enumerate(self.graph.nodes):
			coords_nodes.append(node)
		return coords_nodes
	
	def compute_diameter(self):
		i_min = np.argmin(self.fiedler_vector)
		i_max = np.argmax(self.fiedler_vector)
		coords = self.graph_to_coords() 
		node_min = coords[i_min]
		node_max = coords[i_max]
		diameter_fiedler = nx.dijkstra_path_length(self.graph,node_min,node_max,weight="geometry")
		diameter = nx.diameter(self.graph)
		return diameter, diameter_fiedler
		
	def add_description(self,nbins=100):
		self.description = sd.ShapeDescription(self.fiedler_vector,self,nbins=nbins)
		
	def compute_isolines(self):
		self.description.compute_isolines()
		
	def compute_numpy_coords(self):
		coord_nodes = self.graph_to_coords()
		coords = np.zeros((len(coord_nodes),3),dtype=int)
		for i in range(len(coords)):
			coords[i,:] = coord_nodes[i]
		self.coords = coords
	
	def extract_rostrum(self):
		# Heuristic:
		# - 3/4 of Length
		# - 1/3 bottom of resulting shape
		# - the most at left
		
		#print(coords)
		coords = threshold(self.coords, 1, thres=0.7, inequality = "larger")
		#print(coords)
		coords = threshold(coords, 2, thres=0.3, inequality = "smaller")
		ind = np.argmin(coords[:,1])
		return coords, ind
