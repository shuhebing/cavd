# distutils: language = c++
# distutils: sources = ../libs/Zeo++/networkio.cc

"""
    Updated by Ye Anjiang for extending the In/Output.
    yeanjiang@shu.edu.cn
    May 18, 2019
"""

from netstorage cimport AtomNetwork, VoronoiNetwork
from netstorage cimport ATOM_NETWORK, VORONOI_NETWORK


#Added at 20180704
from libcpp.vector cimport vector
from cavd.voronoicell cimport VOR_CELL, BASIC_VCELL
from cavd.netstorage import PerformVDError
# Define the python definitions for the zeo++ functions

# Easier to implement in python
#def void parseFilename(const char* filename, char* name, char* extension):

#def bint checkInputFile(char* filename)

def readCiffile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readCIFFile(c_filename, atmnet.thisptr, radialflag):
        raise ValueError        # Find the appropriate error and return it
    return atmnet

def readArcfile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readARCFile(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def readCucfile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readCUCFile(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def readCssrfile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readCSSRFile(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def readV1file(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readV1File(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def writeCssrfile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToCSSR(c_filename, c_atmnet):
        raise IOError

def writeCiffile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToCIF(c_filename, c_atmnet):
        raise IOError

def writeV1file(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToV1(c_filename, c_atmnet):
        raise IOError

def writeNt2file(filename, vornet, minRad = None):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    if minRad:
        if not writeToNt2(c_filename, c_vornet_ptr, minRad):
            raise IOError
    else:
        if not writeToNt2(c_filename, c_vornet_ptr):
            raise IOError

def writeXyzfile(filename, atmnet, supercell_flag, is_duplicate_perimeter_atoms):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToXYZ(c_filename, c_atmnet, supercell_flag, 
            is_duplicate_perimeter_atoms):
        raise IOError

def writeVtkfile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToVTK(c_filename, c_atmnet):
        raise IOError

def writeMopacfile(filename, atmnet, supercell_flag):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToMOPAC(c_filename, c_atmnet, supercell_flag):
        raise IOError

def read_from_RemoveMigrantCif(filename,migrant,radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
       filename = (<unicode>filename).encode('utf8')
    if isinstance(migrant, unicode):
       migrant = (<unicode>migrant).encode('utf8')
    cdef char* c_filename = filename
    cdef const char* c_migrant = migrant
    if not readRemoveMigrantCif(c_filename, atmnet.thisptr, c_migrant, radialflag):
        raise ValueError
    return atmnet

# write to .vasp file. 
# Added at 20180426
def writeVaspFile(filename, atmnet, vornet, minRad = None, maxRad = None):
    if isinstance(filename, unicode):
         filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr

    if minRad and maxRad:
        if not writeToVasp(c_filename, c_atmnet, c_vornet_ptr, minRad, maxRad):
            raise IOError
    else:
        if not writeToVasp(c_filename, c_atmnet, c_vornet_ptr):
            raise IOError

# write to atomnetwork to .vasp file. 
# Added at 20180827
def writeAtomNetVaspFile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeAtmntToVasp(c_filename, c_atmnet):
        raise IOError


# write to .net file. 
# Added at 20190518
def writeNETFile(filename, atmnet, vornet, minRad = None, maxRad = None):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    if minRad or maxRad:
        if not writeToNET(c_filename, c_atmnet, c_vornet_ptr, minRad, maxRad):
            raise IOError
    else:
        if not writeToNET(c_filename, c_atmnet, c_vornet_ptr):
            raise IOError