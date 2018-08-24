# distutils: language = c++
# distutils: sources = ../channel.cc
#Added at 20180704
from zeo.netstorage cimport VoronoiNetwork
from zeo.netstorage cimport VORONOI_NETWORK
from channel cimport CHANNEL
from graphstorage cimport DELTA_POS
from graphstorage cimport CONN
from graphstorage cimport DIJKSTRA_NODE
from graphstorage cimport DIJKSTRA_NETWORK
from geometry cimport XYZ
from libcpp.vector cimport vector

#Added at 20180808
#Customize an exception class
class FindChannelError(Exception):
    #print("Find Channel in Voronoi Network Failed!")
    pass

cdef class Channel:
    """
    Python wrapper to Zeo++ Channel.
    """
    def __cinit__(self):
        self.thisptr = new CHANNEL()
    def __dealloc__(self):
        del self.thisptr

    property nodes:
        def __get__(self):
            nodes = []
            cdef vector[DIJKSTRA_NODE] dj_nodes = self.thisptr.nodes
            for i in range(dj_nodes.size()):
                dj_id = dj_nodes[i].id
                dj_coords = list(dj_nodes[i].x,dj_nodes[i].y,dj_nodes[i].z)
                dj_max_radius = dj_nodes[i].max_radius
                nodes.append(list(dj_id, dj_coords, dj_max_radius))
            return nodes

    property connections:
        def __get__(self):
            connections = []
            cdef vector[CONN] conns = self.thisptr.connections
            for i in range(conns.size()):
                conn_from = conns[i].origin
                conn_to = conns[i].ending
                conn_length = conns[i].length
                conn_max_radius = conns[i].max_radius
                conn_delta_pos = list(conns[i].deltaPos.x,conns[i].deltaPos.y,conns[i].deltaPos.z)
                conn = list(conn_from, conn_to, conn_length, conn_max_radius, conn_delta_pos)
                connections.append(conn)
            return connections
    
    property nodes_pos:
        def __get__(self):
            nodes_pos = []
            cdef vector[DELTA_POS] unitCells = self.thisptr.unitCells
            cdef vector[vector[int]] ucNodes = self.thisptr.ucNodes
            for i in range(unitCells.size()):
                cdef DELTA_POS pos = unitCells[i]
                cdef vector[int] ucNode = ucNodes[i]
                for j in range(ucNode.size()):
                    node_id = ucNode[j]
                    node_pos = list(pos.x, pos.y, pos.z)
                    nodes_pos.append(node_id, node_pos)
            return nodes_pos
            
    property dimensionality:
        def __get__(self):
            return self.thisptr.dimensionality
    
    property lattice:
        def __get__(self):
            cdef XYZ v_a = self.thisptr.v_a
            cdef XYZ v_b = self.thisptr.v_b
            cdef XYZ v_c = self.thisptr.v_c
            la = lsit(v_a.x,v_a.y,v_a.z)
            lb = lsit(v_b.x,v_b.y,v_b.z)
            lc = lsit(v_c.x,v_c.y,v_c.z)
            lattice = list(la, lb, lc)
            return lattice
    @classmethod
    def findChannelsInVornet(self, vornet, probe_rad, filename):
        cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')
        cdef char* c_filename = filename
        cdef vector[CHANNEL] channels
        if findChannels_new(c_vornet_ptr, probe_rad, &channels):
            if not c_writeToVMD(channels, c_filename):
                raise IOError
            if not c_writeToNET(channels, c_filename):
                raise IOError
        else:
            raise FindChannelError

#def fincChannelsinDijkstraNett(i

#Add at 20180823
cdef class Channel_Node:
    property id:
		def __get__(self):
		    id = list(self.thisptr.id)
			return id
		def __set__(self, id):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.id = id
    property coords:
        def __get__(self):
		    coords = list(self.thisptr.x,self.thisptr.y,self.thisptr.z)
			return coords
		def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.x = coords[0]
            self.thisptr.y = coords[1]
            self.thisptr.z = coords[2]
	property radius:
        def __get__(self): return self.thisptr.radius
        def __set__(self, radius): 
            print ("This value is not supposed to be modified")
            self.thisptr.radius = radius
			
cdef class Channel_Edge:
	property conn:
		def __get__(self):
		    conn = list(self.thisptr.from, self.thisptr.to)
			return conn
		def __set__(self, from, to):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.from = from
			self.thisptr.to = to

cdef class Channel_new:
    property lattice:
		def __get__(self):
		    lattice = list(self.thisptr.lattice)
			return lattice
		def __set__(self, lattice):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.lattice = lattice
    
