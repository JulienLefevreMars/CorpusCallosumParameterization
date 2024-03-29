import nibabel as nb
import numpy as np
import voxel_spectral_analysis as vsa
import shape_description as sd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import vizu as vz
import sys




if __name__ =="__main__":
	if len(sys.argv)==1:
		name = "238"
	else:
		name = sys.argv[1]
	if len(sys.argv)<=2:
		fig_to_display = "00011"
	else:
		fig_to_display = sys.argv[2]
		
	if len(sys.argv)<=3:
		irregular_bins = False
	else:
		irregular_bins = sys.argv[3]=="true"
		
	if len(fig_to_display)!= 5:
		print("second parameter should be of the form abcde where a,b,c,d,e are 0 or 1 depending wether you want to print a. Fiedler vector b. Isolines c. Skeleton d. profile e. thickness")
	
	data_folder = "/home/julienlefevre/ownCloud/Documents/Recherche/Data/CorpusCallosum/isthme_du_corps_calleux/"
	subject_name = "corpus_callosum_mask_26c_" + name
	filename = data_folder + subject_name  +  ".nii.gz"
	g = nb.load(filename)
	mask = np.asanyarray(g.dataobj)
	(x,y,z) = vsa.image_to_points(mask,g.header['pixdim'][1:4])

	# 1. Visualize the points of the mask
	#vz.visualize_mask(g)

	# 2. Compute Fiedler vector	and extrema
	graph = vsa.image_to_graph(mask,graph_type="geometry")
	# print(graph.nodes)
	diameter_fiedler, diameter, fiedler_vector = vsa.get_diameter_fiedler(graph)
	print("Diameter Fiedler = ", diameter_fiedler)
	print("Diameter  = ",diameter)
	
	if fig_to_display[0] == "1": 
		vz.visualize_fiedler(graph,fiedler_vector,title=subject_name)
	
	# 3. Isolines
	if fig_to_display[1] == "1": 
		vz.visualize_fiedler(graph,sd.compute_isolines(fiedler_vector,nbins=100)[0],title=subject_name)
	
	# 4. Skeleton 
	coords = vsa.graph_to_coords(graph)
	barycenters,intervals,coords, new_fiedler_vector = sd.compute_longitudinal_description(fiedler_vector,coords,nbins=200,
	irregular_bins=irregular_bins)
	if fig_to_display[2] == "1": 
		fig, ax = vz.visualize_fiedler(graph,None,title=subject_name)
		plt.gca().scatter(barycenters[:,0], barycenters[:,1], barycenters[:,2],c='r')
	plt.show()

	# zeros_ls = np.logical_and(new_fiedler_vector >=-0.001,new_fiedler_vector<0.001)
	# vz.visualize_fiedler(graph,zeros_ls,title=subject_name)
	# plt.show()
	#'''
	# 5. Thickness profile
	if irregular_bins:
		n_thickness=30
	else:
		n_thickness=50
	thickness, slices, intervals = sd.compute_thickness(new_fiedler_vector,coords,nbins=n_thickness)
	if fig_to_display[3] == "1": 
		vz.thickness_profile(thickness,subject_name)
	
	if fig_to_display[4] == "1":
		for i in range(len(slices)):
			if i==0:
				fig, ax = vz.visualize_fiedler(slices[i][0],slices[i][1])
			else:
				fig, ax = vz.visualize_fiedler(slices[i][0],slices[i][1],"",fig,ax)
		vz.set_axes_equal(ax,coords)
	
	# 6. Thickness remapped on the image
	#print(len(intervals))
	texture_remapped = sd.texture_mapping(new_fiedler_vector, thickness[:,3], intervals)
	print(thickness, intervals)
	#texture_remapped = sd.texture_mapping(fiedler_vector, np.arange(0,19,1.), intervals)
	if fig_to_display[4] == "1": 
		vz.visualize_fiedler(graph,texture_remapped,title = subject_name)
	plt.show()
	#'''
	


