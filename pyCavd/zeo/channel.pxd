from libcpp.vector cimport vector
from netstorage cimport VORONOI_NETWORK
from graphstorage cimport DIJKSTRA_NETWORK


cdef extern from "../../zeo++/channel.h":
    cdef cppclass CHANNEL:
        CHANNEL() except +

cdef extern from "../../zeo++/channel.h" namespace "CHANNEL":
    #cdef void findChannels(DIJKSTRA_NETWORK*, vector[bint] *, vector[CHANNEL] *)
    cdef bint findChannels_new(VORONOI_NETWORK*, double, vector[CHANNEL] *)

cdef extern from "../../zeo++/channel.h":
    cdef void c_writeToVMD "writeToVMD_new"(vector[CHANNEL] channels, char *filename)

cdef class Channel:
    cdef CHANNEL* thisptr
