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
            c_writeToVMD(channels, c_filename)
        else:
            raise FindChannelError

#def fincChannelsinDijkstraNett(i


	
    
