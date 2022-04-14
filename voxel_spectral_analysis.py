import numpy as np
import networkx as nx

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
					g.add_edge(inds,(x+i,y+j,z+k),weight=weight)
	

def image_to_graph(mask,graph_type="topology"):
	ix, iy, iz = image_to_points(mask)
	g = nx.Graph()
	for i in range(len(ix)):
		add_valid_edge(g, (ix[i],iy[i],iz[i]), mask,graph_type)
	return g
	

