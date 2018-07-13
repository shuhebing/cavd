from libcpp.vector cimport vector
from zeo.netstorage cimport VORONOI_NETWORK
from zeo.graphstorage cimport DIJKSTRA_NETWORK


cdef extern from "../../channel.h":
    cdef cppclass CHANNEL:
        CHANNEL() except +

cdef extern from "../../channel.h" namespace "CHANNEL":
    #cdef void findChannels(DIJKSTRA_NETWORK*, vector[bint] *, vector[CHANNEL] *)
    cdef void findChannels_new(VORONOI_NETWORK*, double, vector[CHANNEL] *)

cdef extern from "../../channel.h":
    cdef void c_writeToVMD "writeToVMD_new"(vector[CHANNEL] channels, char *filename)

cdef class Channel:
    cdef CHANNEL* thisptr
