from scipy.signal import argrelmin, argrelmax
import numpy as np

class AnalyzeCurve:
	def __init__(self,curv_abs,thickness,**kwargs):
		self.curv_abs = curv_abs
		self.thickness = thickness
		self.characteristics = {}

	def total_length(self):
		return self.curv_abs[-1]

	def extract_indices(self,a,b,L):
		'''
		a, b: percentages of total length 
		L: length
		'''
		La = a * L
		Lb = b * L
		inda = np.argmin(np.abs(self.curv_abs - La))
		indb = np.argmin(np.abs(self.curv_abs - Lb))
		return inda, indb

	def local_extrema(self,a = 0,b = 2/5,min_max='min',order = 2):
		'''
		a, b: percentages of total length where to extract a local minimum
		order : size of the window to determine the extremum
		'''
		func = [argrelmin, np.argmin]
		if min_max=='max':
			func = [argrelmax, np.argmax]
		ind_a, ind_b = self.extract_indices(a,b,self.total_length())
		#print("Indices")
		#print(ind_a,ind_b)
		indices = func[0](self.thickness[ind_a:ind_b+1], order = order)
		#print(indices)
		indices = indices[0]
		if indices[0] == 0: # remove 0 if local minimum
			indices = indices[1:]
		indices = [i + ind_a for i in indices]
		good_ind = func[1](self.thickness[indices])
		good_ind = indices[good_ind]
		#print("Good index")
		#print(good_ind)
		return self.curv_abs[good_ind], self.thickness[good_ind], good_ind
		
	def mean_thickness(self,a = 1/2,b = 3/4):
		'''
		Integrate thickness (Trapeze method)
		'''
		ind_a, ind_b = self.extract_indices(a,b,self.total_length())
		y = self.thickness[ind_a:ind_b+1]
		x = self.curv_abs[ind_a:ind_b+1] 
		integrale = np.sum((x[1:] - x[0:-1]) * (y[1:]+y[0:-1])/2)
		return integrale/(x[-1]-x[0])
		
	def characteristic_corpus_callosum(self):
		# Axe des x: Splenium -> Isthme -> Corps -> Genou/Bec
		self.characteristics['min_thickness_isthmus'] = self.local_extrema()
		self.characteristics['mean_thickness_body'] = self.mean_thickness()
		self.characteristics['max_thickness_genou'] = self.local_extrema(a=0.8 , b= 1,min_max = 'max')
		self.characteristics['max_thickness_splenium'] = self.local_extrema(a=0., b= 0.3, min_max = 'max')
		self.characteristics['total_length'] = self.total_length()

		
