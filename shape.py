""""
Definition of a shape class


"""

# Author: Julien Lefevre, PhD, julien.lefevre@univ-amu.fr

import networkx as nx
import nibabel as nb
import numpy as np
import time

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


class Shape:
	def __init__(self,filename=None,graph_type="geometry",**kwargs):
		if not(filename is None):
			g = nb.load(filename)
			mask = np.asanyarray(g.dataobj)
			self.coords_triplet = image_to_points(mask,g.header['pixdim'][1:4])
			self.graph = image_to_graph(mask,graph_type=graph_type)
			self.is_empty = False
			self.fiedler_vector = None
		self.is_empty = True
		
	def get_fiedler(self):
		self.fiedler_vector = nx.fiedler_vector(self.graph)
	
	def graph_to_coords(graph):
		n = len(graph.nodes)
		coords_nodes = []# np.zeros((n,),dtype=int)
		for i,node in enumerate(graph.nodes):
			coords_nodes.append(node)
		return coords_nodes
	
	def compute_diameter(self):
		i_min = np.argmin(self.fiedler_vector)
		i_max = np.argmax(self.fiedler_vector)
		coords = self.graph_to_coords(self.graph) 
		node_min = coords[i_min]
		node_max = coords[i_max]
		diameter_fiedler = nx.dijkstra_path_length(graph,node_min,node_max,weight="geometry")
		diameter = nx.diameter(graph)
		return diameter, diameter_fiedler