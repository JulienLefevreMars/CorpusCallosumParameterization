# CorpusCallosumParameterization

This code allows to parameterize a (elongated) 3D shape defined by a set of voxels and to define a thickness associated to each slice of the shape.

Several steps are required:
1. Build a graph from the set of vertices (26-connectivity)
2. Compute the Fiedler vector of the graph. The weights on the edge depends on the distance $d$ between two neighboring nodes as $exp(-d^2/s^2)$ ($s = 1$ in practice).
3. Compute (pseudo-)isolines of the Fiedler vector. Regular bins are defined between the min and max value of the Fiedler vector. It provides slices of the shape, of irregular width.
4. Barycenters of each slice are computed. 
5. There are several ways to estimate a thickness for each slice, by taking the maximum distance between a barycenter and another point of the slice or by computing 3 standard deviation of the distance to the barycenter. A profile of thickness can be obtained and minima as well
6. A last step consists in backprojecting the value of the thickness on each slice, in particular to check the coherence of the metric with the visual impression.

## Dependencies

Those libraries are mandatory

- nibabel
- numpy
- networkx
- matplotlib

## How to use it ?

1. Ad hoc step: change the directory `data_folder` on line 32 of main.py
2. to run the code, several options:  
- `python main.py` : plot the profile and the thickness for subject 201
- `python main.py subject_name`: same with subject `subject_name` (201, 223, 238, 256... )
- `python main.py subject_name abcde`: where a,b,c,d,e are 0 or 1 depending wether you want to print a. the Fiedler vector b. Its isolines c. the barycenters d. the profile e. thickness


