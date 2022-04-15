import nibabel as nb
import numpy as np
import voxel_spectral_analysis as vsa
import shape_description as sd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import vizu as vz
import sys

def graph_to_coords(graph):
	n = len(graph.nodes)
	coords = np.zeros((n,3))
	for i,node in enumerate(graph.nodes):
		coords[i,:] = node
	return coords

if __name__ =="__main__":
	if len(sys.argv)==1:
		name = "238"
	else:
		name = sys.argv[1]
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
	fiedler_vector = nx.fiedler_vector(graph)

	#vz.visualize_fiedler(graph,fiedler_vector,title=subject_name)
	
	# 3. Isolines
	#vz.visualize_fiedler(graph,sd.compute_isolines(fiedler_vector,nbins=50),title=subject_name)
	
	# 4. Skeleton 
	coords = graph_to_coords(graph)
	barycenters,intervals = sd.compute_longitudinal_description(fiedler_vector,coords,nbins=50)
	print(barycenters)
	fig = vz.visualize_fiedler(graph,None,title=subject_name)
	plt.gca().scatter(barycenters[:,0], barycenters[:,1], barycenters[:,2],c='r')
	plt.show()
	
	# 5. Thickness profile
	thickness = sd.compute_thickness(fiedler_vector,coords,nbins=50)
	plt.figure
	plt.plot(thickness)
	plt.title(subject_name)
	plt.show()
	
	# 6. Thickness remapped on the image
	texture_remapped = sd.texture_mapping(fiedler_vector, thickness, intervals)
	vz.visualize_fiedler(graph,texture_remapped,title = subject_name)
	plt.show()
	
	


