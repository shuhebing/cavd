# distutils: language = c++
# distutils: sources = ../graphstorage.cc
"""
Cython declarations file for Zeo++ graphstorage module.
Declares Zeo++ DIJSTRA_NODE, DIJKSTRA_NETWORK classes and the associated
python wrappers, DijkstraNode, DijkstraNetwork
"""

__author__ = "Bharat Kumar Medasani"
__email__ = "mbkumar@gmail.com"
__date__ = "Jan 6, 2014"

from libcpp.vector cimport vector
from netstorage cimport VORONOI_NETWORK
from geometry cimport XYZ

cdef extern from "../../zeo++/graphstorage.h":
    cdef cppclass DELTA_POS:
        int x, y, z
        DELTA_POS() except +
        DELTA_POS(int, int, int)
    
    cdef cppclass CONN:
        int origin "from"
        int ending "to"
        double length
        double max_radius
        DELTA_POS deltaPos
        CONN() except +
        CONN(int, int, double, double, int, int, int)
    
    cdef cppclass DIJKSTRA_NODE:
        int id
        double x, y, z
        vector[CONN] connections
        double max_radius
        bint active
        DIJKSTRA_NODE() except +

    cdef cppclass DIJKSTRA_NETWORK:
        vector[DIJKSTRA_NODE] nodes
        XYZ v_a, v_b, v_c
        DIJKSTRA_NETWORK() except +

cdef extern from "../../zeo++/graphstorage.h" namespace "DIJKSTRA_NETWORK":
    cdef void buildDijkstraNetwork(VORONOI_NETWORK*, DIJKSTRA_NETWORK*)


#cdef class DijkstraNode:
#    """
#    Cython wrapper class for Zeo++ DIJKSTRA_NODE class.
#    """
#    cdef DIJKSTRA_NODE* thisptr

cdef class DeltaPos:
   cdef DELTA_POS* thisptr

cdef class Conn:
   cdef CONN* thisptr

cdef class DijkstraNode:
   cdef DIJKSTRA_NODE* thisptr

cdef class DijkstraNetwork:
    """
    Cython wrapper class for Zeo++ DIJKSTRA_NETWORK class.
    """
    cdef DIJKSTRA_NETWORK* thisptr
