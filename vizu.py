import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def set_axes_equal(ax,coords=None):
	'''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''
	# Taken from https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
	
	
	if coords is None :
		
		x_limits = ax.get_xlim3d()
		y_limits = ax.get_ylim3d()
		z_limits = ax.get_zlim3d()
	else:
		x_limits = [np.min(coords[:,0]),np.max(coords[:,0])]
		y_limits = [np.min(coords[:,1]),np.max(coords[:,1])]
		z_limits = [np.min(coords[:,2]),np.max(coords[:,2])]
	
	x_range = abs(x_limits[1] - x_limits[0])
	x_middle = np.mean(x_limits)
	y_range = abs(y_limits[1] - y_limits[0])
	y_middle = np.mean(y_limits)
	z_range = abs(z_limits[1] - z_limits[0])
	z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
	plot_radius = 0.5*max([x_range, y_range, z_range])
	ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
	ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
	ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

'''
def visualize_mask(gifti_image):
	mask = np.asanyarray(gifti_image.dataobj)
	(x,y,z) = vsa.image_to_points(mask,gifti_image.header['pixdim'][1:4])
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x,y,z)
	set_axes_equal(ax)
	#plt.show()
'''
	
def visualize_fiedler_extrema(coords,fiedler_vector,ax):
	ind_min = np.argmin(fiedler_vector)
	ind_max = np.argmax(fiedler_vector)
	ax.scatter(coords[ind_min,0],coords[ind_min,1],coords[ind_min,2],c="y",s=100)
	ax.scatter(coords[ind_max,0],coords[ind_max,1],coords[ind_max,2],c="y",s=100)
	
	
def visualize_fiedler(graph,fiedler_vector=None,title="",fig=None,ax=None,extrema=False):
	n = len(graph.nodes)
	coords = np.zeros((n,3))
	for i,node in enumerate(graph.nodes):
		coords[i,:] = node
	if fig is None:
		fig = plt.figure(figsize=(15,9))
	#ax = fig.add_subplot(111, projection='3d')
		ax = plt.axes(projection="3d")
	#ax.set_box_aspect((1,1,1))
	ax.scatter(coords[:,0],coords[:,1],coords[:,2],c=fiedler_vector,cmap="jet",s=2)
	ax.view_init(30,0)
	if extrema:
		visualize_fiedler_extrema(coords,fiedler_vector,ax)
	plt.title(title)
	set_axes_equal(ax)
	return fig,ax
	#plt.show()
	
def thickness_profile(thickness,subject_name,legend=['3 std','max','fiedler','pca']):
	plt.figure()
	for i in range(thickness.shape[1]):
		plt.plot(range(len(thickness)),thickness[:,i])
		ind_m = np.argmin(thickness[:,i])	
		plt.scatter(ind_m,thickness[ind_m,i],c='k',s=50)
	plt.legend(legend)
	plt.title(subject_name)
	plt.show()
	
def thickness_profile_isometric(thickness,cum_length,subject_name,legend=['3 std','max','fiedler','pca']):
	plt.figure()
	for i in range(thickness.shape[1]):
		plt.plot(cum_length,thickness[:,i])
		ind_m = np.argmin(thickness[:,i])	
		plt.scatter(cum_length[ind_m],thickness[ind_m,i],c='k',s=50)
	plt.legend(legend)
	plt.title(subject_name)
	plt.show()
	

		


	
