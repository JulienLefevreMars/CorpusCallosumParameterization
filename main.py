import nibabel as nb
import numpy as np
import voxel_spectral_analysis as vsa
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
	data_folder = "/home/julienlefevre/ownCloud/Documents/Recherche/Data/CorpusCallosum/isthme_du_corps_calleux/"
	subject_name = "corpus_callosum_mask_26c_" + name
	filename = data_folder + subject_name  +  ".nii.gz"
	g = nb.load(filename)
	mask = np.asanyarray(g.dataobj)
	(x,y,z) = vsa.image_to_points(mask,g.header['pixdim'][1:4])
	
	# 1. Visualize the points of the mask
	#vz.visualize_mask(g)

	# 2. Compute Fiedler vector	and extrema
	graph = vsa.image_to_graph(mask)
	fiedler_vector = nx.fiedler_vector(graph)

	vz.visualize_fiedler(graph,fiedler_vector,title=subject_name)

	


