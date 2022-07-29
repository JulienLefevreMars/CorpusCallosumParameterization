import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import nibabel as nb
import numpy as np
import shape as sh
import vizu as vz
import sys

if __name__ =="__main__":
	if len(sys.argv)==1:
		name = "238"
	else:
		name = sys.argv[1]
	if len(sys.argv)<=2:
		fig_to_display = "11111"
	else:
		fig_to_display = sys.argv[2]
		
	if len(fig_to_display)!= 5:
		print("second parameter should be of the form abcde where a,b,c,d,e are 0 or 1 depending wether you want to display a. Fiedler vector b. Isolines c. Skeleton d. profile e. thickness")
	
	# 1. Open the data
	data_folder = "/home/julienlefevre/ownCloud/Documents/Recherche/Data/CorpusCallosum/isthme_du_corps_calleux/"
	subject_name = "corpus_callosum_mask_26c_" + name
	filename = data_folder + subject_name  +  ".nii.gz"
	shape = sh.Shape(filename = filename,graph_type = "topolgy",nbins=50)
	
	# 2. Compute Fiedler vector	and extrema
	shape.get_fiedler()
	shape.add_description(nbins=50)
	#diameter, diameter_fiedler = shape.compute_diameter()
	#print("Diameter Fiedler = ", diameter_fiedler)
	#print("Diameter  = ",diameter)
	
	if fig_to_display[0] == "1": 
		vz.visualize_fiedler(shape.graph,shape.fiedler_vector,title=subject_name)
	plt.show()

	# 3. Isolines
	shape.compute_isolines()
	if fig_to_display[1] == "1": 
		vz.visualize_fiedler(shape.graph,shape.description.isolines,title=subject_name)
	plt.show()
	

	# 4. Skeleton 
	shape.description.compute_skeleton(add_extremity = True)
	barycenters = shape.description.barycenters
	print(barycenters)
	if fig_to_display[2] == "1": 
		fig, ax = vz.visualize_fiedler(shape.graph,None,title=subject_name)
		ax.scatter(barycenters[:,0], barycenters[:,1], barycenters[:,2],c='r',s=100,marker='o')
	plt.show()
	
	shape.description.reparameterize_texture()
	if fig_to_display[3] == "1":
		fig,ax = vz.visualize_fiedler(shape.graph,shape.description.texture,title=subject_name)
		shape.description.compute_skeleton(add_extremity = True)
		barycenters = shape.description.barycenters
		print(barycenters)
		ax.scatter(barycenters[:,0], barycenters[:,1], barycenters[:,2],c='r',s=100,marker='o')
		plt.show()
	'''	
	# 5. Re-parametrization
	
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
	
	plt.plot(all_a)
	plt.show()
	for i in range(len(cum_length)-1):
			x = intervals[i:i+2]
			plt.plot(x, all_a[i]*x + all_b[i])
			#plt.plot(i,all_a[i])
	plt.show()
	
	vz.visualize_fiedler(graph,reparam_fiedler,title=subject_name)
	plt.show()
	# zeros_ls = np.logical_and(new_fiedler_vector >=-0.001,new_fiedler_vector<0.001)
	# vz.visualize_fiedler(graph,zeros_ls,title=subject_name)
	# plt.show()
	#'''
	'''
	# 6. Thickness profile
	if irregular_bins:
		n_thickness=30
	else:
		n_thickness=50
	thickness, slices, intervals = sd.compute_thickness(reparam_fiedler,coords,nbins=n_thickness)
	if fig_to_display[3] == "1": 
		#vz.thickness_profile_isometric(thickness,cum_length,subject_name)
		vz.thickness_profile(thickness,subject_name)
	
	if fig_to_display[4] == "1":
		for i in range(len(slices)):
			if i==0:
				fig, ax = vz.visualize_fiedler(slices[i][0],slices[i][1])
			else:
				fig, ax = vz.visualize_fiedler(slices[i][0],slices[i][1],"",fig,ax)
		vz.set_axes_equal(ax,coords)
	
	# 7. Thickness remapped on the image
	#print(len(intervals))
	texture_remapped = sd.texture_mapping(reparam_fiedler, thickness[:,3], intervals)
	print(thickness, intervals)
	#texture_remapped = sd.texture_mapping(fiedler_vector, np.arange(0,19,1.), intervals)
	if fig_to_display[4] == "1": 
		vz.visualize_fiedler(graph,texture_remapped,title = subject_name)
	plt.show()
	#'''
