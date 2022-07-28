""""
Definition of a class to represent shape descriptions


"""

import copy
import numpy as np
#import sklearn.decomposition as sd
#import matplotlib.pyplot as plt

# Author: Julien Lefevre, PhD, julien.lefevre@univ-amu.fr

class ShapeDescription:
	def __init__(self,texture,**kwargs):
		if texture is None:
			self.init_texture = None
			self.texture = None
		else:
			self.init_texture = texture # will be fiedler_vector of a graph
			self.texture = copy.deepcopy(texture)
		self.isolines = None
		self.intervals = None
		
	def compute_isolines(self,nbins=100):
		vmin = np.min(self.texture)
		vmax = np.max(self.texture)
		self.isolines = np.zeros((len(self.texture),))
		self.intervals = np.linspace(vmin,vmax,nbins)
		for i in range(0,len(self.intervals)-1,2):
			self.isolines[np.logical_and(self.texture >=self.intervals[i],self.texture<self.intervals[i+1])] = i


