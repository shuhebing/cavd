from libcpp.vector cimport vector
from libcpp.map cimport map
from netstorage cimport VORONOI_NETWORK
from netstorage cimport ATOM_NETWORK
from graphstorage cimport DELTA_POS
from graphstorage cimport CONN
from graphstorage cimport DIJKSTRA_NODE
from graphstorage cimport DIJKSTRA_NETWORK
from geometry cimport XYZ

cdef extern from "../../zeo++/channel.h":
    cdef cppclass CHANNEL:
        map[int,int] idMappings
        map[int,int] reverseIDMappings
        vector[DIJKSTRA_NODE] nodes
        vector[CONN] connections
        
        vector[DELTA_POS] unitCells
        vector[vector[int]] ucNodes
        XYZ v_a, v_b, v_c
        int dimensionality

        CHANNEL() except +

cdef extern from "../../zeo++/channel.h" namespace "CHANNEL":
    #cdef void findChannels(DIJKSTRA_NETWORK*, vector[bint] *, vector[CHANNEL] *)
    cdef bint findChannels_new(VORONOI_NETWORK*, double, vector[CHANNEL] *)

cdef extern from "../../zeo++/channel.h":
    cdef bint c_writeToVMD "writeToVMD_new"(vector[CHANNEL] channels, char *filename)
    cdef bint c_writeToNET "writeToNET_new"(vector[CHANNEL] channels, char *filename, ATOM_NETWORK *cell)

cdef class Channel:
    cdef CHANNEL* thisptr
