""""
Definition of a class to represent shape descriptions


"""

import copy
import matplotlib.pyplot as plt
import numpy as np
import shape as sh
#import sklearn.decomposition as sd


# Author: Julien Lefevre, PhD, julien.lefevre@univ-amu.fr

class ShapeDescription:
	def __init__(self,texture,shape=None,nbins=100,**kwargs):
		if texture is None:
			self.init_texture = None
			self.texture = None
		else:
			self.init_texture = texture # will be fiedler_vector of a graph
			self.texture = copy.deepcopy(texture) # reparameterized texture
		self.shape = shape
		self.isolines = None
		self.intervals = None
		self.barycenters = None
		self.nbins = nbins
		
	def compute_isolines(self):
		vmin = np.min(self.texture)
		vmax = np.max(self.texture)
		self.isolines = np.zeros((len(self.texture),))
		self.intervals = np.linspace(vmin,vmax,self.nbins)
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
		self.barycenters = barycenters

	def reparameterize_texture(self):
		self.compute_skeleton(add_extremity=True)
		# Length
		length = np.sqrt(np.sum((self.barycenters[1:,:] - self.barycenters[0:-1,:])**2,axis=1))
		cum_length = np.cumsum(length)
		cum_length = cum_length/cum_length[-1]
		#plt.plot(cum_length)
		#plt.show()
		# Reparameterization
		reparam_fiedler = np.zeros((len(self.texture),))
		for i in range(len(length)-1):
			indices = np.logical_and(self.texture >= self.intervals[i],self.texture <= self.intervals[i+1] )
			a = (cum_length[i+1] - cum_length[i])/(self.intervals[i+1] - self.intervals[i])
			b = cum_length[i] - a * self.intervals[i]
			reparam_fiedler[indices] = a * self.texture[indices] + b
		self.texture = reparam_fiedler
		self.compute_isolines()
	'''
		# Length
	barycenters = np.vstack([coords[np.argmin(new_fiedler_vector),:],barycenters])
	length = np.sqrt(np.sum((barycenters[1:,:] - barycenters[0:-1,:])**2,axis=1))
	cum_length = np.cumsum(length)
	cum_length = cum_length/cum_length[-1]
	plt.plot(cum_length)
	plt.show()
	
	# Reparam
	reparam_fiedler = np.zeros((len(fiedler_vector),))
	print(len(intervals),len(length), len(barycenters))
	all_a = []
	all_b = []
	for i in range(len(length)-1):
		indices = np.logical_and(new_fiedler_vector >= intervals[i],new_fiedler_vector <= intervals[i+1] )
		a = (cum_length[i+1] - cum_length[i])/(intervals[i+1] - intervals[i])
		b = cum_length[i] - a * intervals[i]
		reparam_fiedler[indices] = a * new_fiedler_vector[indices] + b
		all_a.append(a)
		all_b.append(b)
		'''
