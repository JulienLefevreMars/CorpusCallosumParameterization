import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import nibabel as nb
import numpy as np
import os
import shape as sh
import vizu as vz
import scipy.ndimage as sp
import sys
import analyze_curve as ac

use_fiedler = False # use Fiedler diameter to compute thickness
nbins = 75 # number of bins to obtain slices of Fiedler vector
 # size of window to smooth the thickness curves (not in .csv files)
graph_type = "topology" # "geometry" # or 
data_folder = "/home/julienlefevre/ownCloud/Documents/Recherche/Data/CorpusCallosum/isthme_du_corps_calleux/"


def analyse_profiles(abs_curv=True,sigma_smooth = 1.):
	# Process .csv files
	dir_list = os.listdir(data_folder)
	nb_subj = 0
	data_subjects = []
	data_names = []
	for filename in dir_list:
		pos = filename.find("_thickness")
		if filename[-3::]=="csv" and not(filename[pos-3:pos]=="612"):
			data = np.loadtxt(data_folder + filename,delimiter=",")
			nb_subj +=1
			data_subjects.append(data)
			data_names.append(filename[pos-3:pos])
	all_extrema = []
	for i in range(nb_subj):
		data = data_subjects[i]
		if abs_curv:
			xvalues = data[:,0]
			xlabel = "Skeleton length"
		else:
			xvalues = np.arange(len(data[:,0]))
			xlabel = "Incremental indexing"
		curve = ac.AnalyzeCurve(xvalues,sp.gaussian_filter1d(data[:,1],sigma = sigma_smooth))
		#print(curve.total_length())
		#print(curve.local_extrema())
		#print(curve.mean_thickness())
		curve.characteristic_corpus_callosum()
		all_extrema.append(curve.characteristics)
		plt.plot(xvalues,curve.thickness)
		#plt.plot(xvalues,data[:,1])
	plt.legend(data_names)
	plt.xlabel(xlabel)
	plt.ylabel("Thickness (a.u.)")
	for i in range(nb_subj):
		print(all_extrema[i].items())
		for r in (all_extrema[i].items()):
			print(r)
			if not(type(r[1])==np.float64): # (position, value, index)
				plt.plot(r[1][0],r[1][1],'.',markersize=10,color='k')
	plt.show()
	save_features(all_extrema,data_folder + "all_characteristics.csv")

def save_features(all_extrema,filename,delimiter=","):
	nb_subj = len(all_extrema)
	nb_features = len(all_extrema[0])
	
	row1 = ""
	for key in all_extrema[0].items():
		if not(type(key[1]) == np.float64):
			row1 += key[0] + " (position) " + delimiter
			row1 += key[0] + " (value) " + delimiter
		else:
			row1 += key[0] + " " + delimiter 
	data = []
	for i in range(nb_subj):
		row = []
		for key in all_extrema[i].items():
			if not(type(key[1]) == np.float64):
				row.append(key[1][0])
				row.append(key[1][1])
			else:
				row.append(key[1])
		data.append(row)
	print(data)
	with open(filename, 'w') as f:
		f.writelines(row1 + '\n')
		np.savetxt(f, data[0:] , delimiter = delimiter, fmt ='%f')		


def save_profiles(thickness,shape,filename):
	length = shape.description.curv_abs()
	n = len(length)-1
	data = np.reshape(np.cumsum(length[0:-1]),(n,1))
	data = np.concatenate((data,np.reshape(thickness,(n,1))),axis=1)
	np.savetxt(filename, data , delimiter = ",")
	#print(length)
	#plt.plot(np.cumsum(length[0:-1]),thickness)
	#plt.show()
	
def process_subject(name,fig_to_display):
	# 1. Open the data
	subject_name = "corpus_callosum_mask_26c_" + name
	filename = data_folder + subject_name  +  ".nii.gz"
	shape = sh.Shape(filename = filename,graph_type = graph_type)
	#print(shape.graph_to_coords())
	
	# 2. Compute Fiedler vector	and extrema
	shape.get_fiedler()
	shape.add_description(nbins=nbins)
	#diameter, diameter_fiedler = shape.compute_diameter()
	#print("Diameter Fiedler = ", diameter_fiedler)
	#print("Diameter  = ",diameter)
	
	if fig_to_display[0] == "1": 
		title = subject_name + "\n \n Fiedler vector"
		vz.visualize_fiedler(shape.graph,shape.fiedler_vector,title=title,extrema=True)
	plt.show()

	# 3. Isolines
	shape.compute_isolines()
	if fig_to_display[1] == "1": 
		title = subject_name + "\n \n Isolines of Fiedler vector"
		vz.visualize_fiedler(shape.graph,shape.description.isolines,title=title)
	plt.show()
	

	# 4. Skeleton 
	shape.description.compute_skeleton(add_extremity = True)
	barycenters = shape.description.barycenters
	print(barycenters)
	if fig_to_display[2] == "1": 
		title = subject_name + "\n \n Skeleton"
		fig, ax = vz.visualize_fiedler(shape.graph,None,title=title)
		ax.scatter(barycenters[:,0], barycenters[:,1], barycenters[:,2],c='r',s=100,marker='o')
		ax.view_init(0,0)
	plt.show()
	
	for i in range(10): # HOW TO SET THAT ?
		shape.description.reparameterize_texture()
		
	if fig_to_display[1] == "1": 
		title = subject_name + "\n \n Isolines of reparameterized Fiedler vector"
		vz.visualize_fiedler(shape.graph,shape.description.isolines,title=title)
	plt.show()
	'''
	if fig_to_display[3] == "1":
		fig,ax = vz.visualize_fiedler(shape.graph,shape.description.texture,title=subject_name)
		shape.description.compute_skeleton(add_extremity = True)
		barycenters = shape.description.barycenters
		print(barycenters)
		ax.scatter(barycenters[:,0], barycenters[:,1], barycenters[:,2],c='k',s=100,marker='o')
		plt.show()
	'''


	# 6. Thickness profile
	thickness = shape.description.compute_thickness(use_fiedler=False)
	print(len(barycenters),len(thickness))
	if fig_to_display[3] == "1": 
		#vz.thickness_profile_isometric(thickness,cum_length,subject_name)
		vz.thickness_profile(thickness,subject_name)

	# 7. Thickness remapped on the image
	#print(len(intervals))
	texture_remapped = shape.description.texture_mapping(thickness[:,2])
	print(texture_remapped)
	if fig_to_display[4] == "1": 
		title = subject_name + "\n \n Thickness, by slice"
		vz.visualize_fiedler(shape.graph,texture_remapped,title = title)
	plt.show()
	
	# 8. Save the profile(s)
	# csv file
	save_profiles(thickness[:,2],shape,data_folder + subject_name + "_thickness.csv")


if __name__ =="__main__":
	if len(sys.argv)==1:
		analyse_profiles(True) # True: skeleton length in x
	else:
		# Process one subject
		name = sys.argv[1]
		if len(sys.argv)<=2:
			fig_to_display = "11111"
		else:
			fig_to_display = sys.argv[2]
		
		if len(fig_to_display)!= 5:
			print("second parameter should be of the form abcde where a,b,c,d,e are 0 or 1 depending wether you want to display a. Fiedler vector b. Isolines c. Skeleton d. Thickness profile e. thickness remapped on the CC")
		process_subject(name,fig_to_display)

