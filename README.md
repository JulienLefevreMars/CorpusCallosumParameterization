# CorpusCallosumParameterization

This code allows to parameterize a (elongated) 3D shape defined by a set of voxels and to define a thickness associated to each slice of the shape.

Several steps are required:
1. Build a graph from the set of vertices (26-connectivity)
2. Compute the Fiedler vector of the graph. The weights on the edge depends on the distance $d$ between two neighboring nodes as $exp(-d^2/\sigma^2)$
