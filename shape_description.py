""""
Definition of a class to represent shape descriptions


"""

import copy
import numpy as np
import shape as sh
#import sklearn.decomposition as sd
#import matplotlib.pyplot as plt

# Author: Julien Lefevre, PhD, julien.lefevre@univ-amu.fr

class ShapeDescription:
	def __init__(self,texture,shape=None,**kwargs):
		if texture is None:
			self.init_texture = None
			self.texture = None
		else:
			self.init_texture = texture # will be fiedler_vector of a graph
			self.texture = copy.deepcopy(texture)
		self.shape = shape
		self.isolines = None
		self.intervals = None
		
	def compute_isolines(self,nbins=100):
		vmin = np.min(self.texture)
		vmax = np.max(self.texture)
		self.isolines = np.zeros((len(self.texture),))
		self.intervals = np.linspace(vmin,vmax,nbins)
		for i in range(0,len(self.intervals)-1,2):
			self.isolines[np.logical_and(self.texture >=self.intervals[i],self.texture<self.intervals[i+1])] = i
		
	def compute_skeleton(self,add_extremity=False):	
		'''
		Provides barycenters/skeleton of a shape by using isolines of the texture
		'''
		nbins = len(self.intervals)
		coords_nodes = self.shape.graph_to_coords()
		coords = np.zeros((len(coords_nodes),3),dtype=int)
		for i in range(len(coords)):
			coords[i,:] = coords_nodes[i]
		barycenters = np.zeros((nbins-1,3))
		for i in range(0,nbins-1):
			barycenters[i,:] = np.mean(coords[np.logical_and(self.texture>=self.intervals[i],self.texture<self.intervals[i+1])],axis=0)
		if add_extremity:
			barycenters = np.vstack([coords[np.argmin(self.texture),:],barycenters])
			barycenters = np.vstack([barycenters,coords[np.argmax(self.texture),:]])
		return barycenters
