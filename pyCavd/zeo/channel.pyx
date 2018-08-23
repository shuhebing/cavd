# distutils: language = c++
# distutils: sources = ../channel.cc
#Added at 20180704
from zeo.netstorage cimport VoronoiNetwork
from zeo.netstorage cimport VORONOI_NETWORK
from channel cimport CHANNEL
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
    
