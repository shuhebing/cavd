from channel cimport CHANNEL
from netstorage cimport ATOM_NETWORK, VORONOI_NETWORK
from graphstorage cimport DIJKSTRA_NETWORK

cdef extern from "../basic_lib/Zeo++/feature.h":
    cdef cppclass FEATURE (CHANNEL):
        FEATURE(vector[int] nodeIds, DIJKSTRA_NETWORK* dnet,  int dim,
                int basisVecs[3][3])
        # C++ code needs modification
        int createSegments(ATOM_NETWORK*, VORONOI_NETWORK*, DIJKSTRA_NETWORK,
                char* filename, int initIndex)
