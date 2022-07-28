""""
Definition of a class to represent shape descriptions


"""

# Author: Julien Lefevre, PhD, julien.lefevre@univ-amu.fr

class ShapeDescription:
	def __init__(self,texture,**kwargs):
		self.init_texture = texture # will be fiedler_vector of a graph
		self.isolines = None
		self.intervals = None
		self.texture = copy.deepcopy(texture)
		
	def compute_isolines(texture, nbins=100):
		vmin = np.min(texture)
		vmax = np.max(texture)
		self.isolines = np.zeros((len(texture),))
		self.intervals = np.linspace(vmin,vmax,nbins)
		for i in range(0,len(intervals)-1,2):
			self.isolines[np.logical_and(texture >=intervals[i],texture<intervals[i+1])] = i


