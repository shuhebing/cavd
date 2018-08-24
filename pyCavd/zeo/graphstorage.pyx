# distutils: language = c++
# distutils: sources = ../graphstorage.cc

from netstorage cimport VoronoiNetwork

cdef class DijkstraNetwork:
    """
    Python wrapper class to Zeo++ Djikstra Network
    """
    #cdef DIJKSTRA_NETWORK* thisptr
    def __cinit__(self):
        self.thisptr = new DIJKSTRA_NETWORK()
    @classmethod
    def from_VoronoiNetwork(vornet):
        """
        Build Dijkstra Net from input Voronoi Net
        """
        dijkstranet = DijkstraNetwork()
        c_vornet = (<VoronoiNetwork?>vornet).thisptr
        buildDijkstraNetwork(c_vornet, dijkstranet.thisptr)
        return dijkstranet 
    def __dealloc__(self):
        del self.thisptr

cdef class DeltaPos:
    def __cinit__(self, int x, int y, int z):
        self.thisptr = new DELTA_POS(x,y,z)

    def __init__(self, int x, int y, int z):
        pass

    def __dealloc__(self):
        del self.thisptr

    property pos:
        def __get__(self):
            pos = list(self.thisptr.x, self.thisptr.y, self.thisptr.z)
            return pos
        def __set__(self, pos):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.x = pos[0]
            self.thisptr.y = pos[1]
            self.thisptr.z = pos[2]

cdef class Conn:
    def __cinit__(self):
        self.thisptr = new Conn()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    property origin:
        def __get__(self):
            origin = self.thisptr.origin
            return origin
        def __set__(self, origin):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.origin = origin
    
    property ending:
        def __get__(self):
            ending = self.thisptr.ending
            return ending
        def __set__(self, ending):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.ending = ending
    
    property length:
        def __get__(self):
            length = self.thisptr.length
            return length
        def __set__(self, length):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.length = length
    
    property max_radius:
        def __get__(self):
            max_radius = self.thisptr.max_radius
            return max_radius
        def __set__(self, max_radius):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.max_radius = max_radius

    property pos:
        def __get__(self):
            pos = list(self.thisptr.deltaPos.x, self.thisptr.deltaPos.y, self.thisptr.deltaPos.z)
            return pos
        def __set__(self, pos):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.deltaPos.x = pos[0]
            self.thisptr.deltaPos.y = pos[1]
            self.thisptr.deltaPos.z = pos[2]

cdef class Conn:
    def __cinit__(self):
        self.thisptr = new Conn()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    property origin:
        def __get__(self):
            origin = self.thisptr.origin
            return origin
        def __set__(self, origin):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.origin = origin
    
    property coord:
        def __get__(self):
            ending = self.thisptr.ending
            return ending
        def __set__(self, ending):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.ending = ending
    
    property length:
        def __get__(self):
            length = self.thisptr.length
            return length
        def __set__(self, length):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.length = length
    
    property max_radius:
        def __get__(self):
            max_radius = self.thisptr.max_radius
            return max_radius
        def __set__(self, max_radius):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.max_radius = max_radius

    property pos:
        def __get__(self):
            pos = list(self.thisptr.deltaPos.x, self.thisptr.deltaPos.y, self.thisptr.deltaPos.z)
            return pos
        def __set__(self, pos):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.deltaPos.x = pos[0]
            self.thisptr.deltaPos.y = pos[1]
            self.thisptr.deltaPos.z = pos[2]

cdef class DIJKSTRA_NODE:
    def __cinit__(self):
        self.thisptr = new DIJKSTRA_NODE()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    property id:
        def __get__(self):
            id = self.thisptr.id
            return id
        def __set__(self, id):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.id = id
    
    property coords:
        def __get__(self):
            coords = list(self.thisptr.x, self.thisptr.y, self.thisptr.z)
            return coords
        def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.x = coords[0]
            self.thisptr.y = coords[1]
            self.thisptr.z = coords[2]
    
    property max_radius:
        def __get__(self):
            max_radius = self.thisptr.max_radius
            return max_radius
        def __set__(self, max_radius):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.max_radius = max_radius

    property connections:
        def __get__(self):
            cdef vector[CONN] conns = self.thisptr.connections
            connections = []
            for i in range(conns.size()):
                conn_from = conns[i].origin
                conn_to = conns[i].ending
                conn_length = conns[i].length
                conn_max_radius = conns[i].max_radius
                conn_delta_pos = list(conns[i].deltaPos.x,conns[i].deltaPos.y,conns[i].deltaPos.z)
                conn = list(conn_from, conn_to, conn_length, conn_max_radius, conn_delta_pos)
                connections.append(conn)
            return connections
        