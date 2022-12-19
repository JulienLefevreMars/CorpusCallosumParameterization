""""
Definition of a class to represent shape descriptions


"""

import copy
import matplotlib.pyplot as plt
import numpy as np
import shape as sh
import sklearn.decomposition as sd


# Author: Julien Lefevre, PhD, julien.lefevre@univ-amu.fr

def distance(p,points):
	return np.sqrt((points[:,0]-p[0])**2 + (points[:,1]-p[1])**2 + (points[:,2]-p[2])**2)


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
		self.coords = None
		self.cum_length = None
		
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
		self.coords = coords
		barycenters = np.zeros((nbins-1,3))
		for i in range(0,nbins-1):
			barycenters[i,:] = np.mean(coords[np.logical_and(self.texture>=self.intervals[i],self.texture<self.intervals[i+1])],axis=0)
		if add_extremity:
			barycenters = np.vstack([coords[np.argmin(self.texture),:],barycenters])
			barycenters = np.vstack([barycenters,coords[np.argmax(self.texture),:]])
		# Remove nan
		barycenters = barycenters[~np.isnan(barycenters).any(axis=1)]
		self.barycenters = barycenters
		self.nbins = len(self.barycenters)-1
		self.compute_isolines()

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
		#print(self.texture[1225])
		print(self.intervals)
		print("Length intervals " + str(len(self.intervals)))
		print("Length barycenters " + str(len(length)))
		for i in range(len(length)-1):
			indices = np.logical_and(self.texture >= self.intervals[i],self.texture <= self.intervals[i+1] )
			a = (cum_length[i+1] - cum_length[i])/(self.intervals[i+1] - self.intervals[i])
			b = cum_length[i] - a * self.intervals[i]
			reparam_fiedler[indices] = a * self.texture[indices] + b
		self.texture = reparam_fiedler
		self.compute_isolines()
		self.compute_skeleton(add_extremity=True)
		self.cum_length = cum_length

	def compute_thickness(self,use_fiedler=False):
		nbins = len(self.intervals)
		nb_methods = 3
		if use_fiedler:
			nb_methods = 4
		thickness = np.zeros((nbins-1,nb_methods))
		barycenters = np.zeros((nbins-1,3))
		for i in range(0,len(self.intervals)-1):
			print(i)
			print(self.intervals[i:i+2])
			# 1. Skeletonization => # Could be refactored with compute_skeleton
			indices = np.logical_and(self.texture >=self.intervals[i],self.texture<=self.intervals[i+1])
			# warning, if <bins[i+1] indices is empty for i = len(bins)-1 because of histogram equalization that puts a lot of values on vmax
			slice_points = self.coords[indices]
			barycenters[i,:] = np.mean(slice_points,axis=0)
			thickness[i,0] = 3*np.std(distance(barycenters[i,:], self.coords[indices,:]))
			thickness[i,1] = np.max(distance(barycenters[i,:], self.coords[indices,:]))
			# 2. PCA
			pca = sd.PCA(n_components=2)
			pca.fit(slice_points)
			projection = pca.transform(slice_points)
			thickness[i,2] = np.max(projection) - np.min(projection)
			# 3. Convert a slice in a graph and get Fiedler diameter
			if use_fiedler:
				slice_graph = sh.points_to_graph(slice_points,graph_type="geometry")
				#print(slice_graph)
				try:
					shape = sh.Shape(filename=None,graph = slice_graph)
					res = shape.compute_diameter()
				except :
					print("Problem with networkx")
					res = [np.nan,np.nan]
				thickness[i,3] = res[0]
		return thickness
		
	def texture_mapping(self,thickness):
		texture_mapped = np.zeros((len(self.texture),))
		for i in range(0,len(thickness)):
			indices = np.logical_and(self.texture >=self.intervals[i],self.texture<=self.intervals[i+1])
			texture_mapped[indices] = thickness[i]
		return texture_mapped
		
	def curv_abs(self):
		n_bar = len(self.barycenters)
		length = np.zeros((n_bar-1,))
		for i in range(n_bar-1):
			length[i] = distance(self.barycenters[i+1,:],self.barycenters[i:i+1,:])
		return length
