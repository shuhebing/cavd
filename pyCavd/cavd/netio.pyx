# distutils: language = c++
# distutils: sources = ../networkio.cc

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

#class PerformVDError(Exception):
#    #print("Perform Voronoi Decompition failured!")
#    pass

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

# write to .bi file. Added at 20180408
def writeBIFile(filename, atmnet, vornet, minRad = None):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    if minRad:
        if not writeToBI(c_filename, c_atmnet, c_vornet_ptr, minRad):
            raise IOError
    else:
        if not writeToBI(c_filename, c_atmnet, c_vornet_ptr):
            raise IOError

# remove migrant ion. Added at 20180408
def getRemoveMigrantFilename(filename,migrant):
    if isinstance(filename, unicode):
       filename = (<unicode>filename).encode('utf8')
    if isinstance(migrant, unicode):
       migrant = (<unicode>migrant).encode('utf8')
    cdef char* c_filename = filename
    cdef const char* c_migrant = migrant
    pretreatedFilename = pretreatCifFilename(c_filename,c_migrant)
    if pretreatedFilename == "":
        #raise IOError("Can't Open ", filename, " or Can't Write to outputfile.")
        raise IOError
    else:
        return pretreatedFilename

# write to .vasp file. Added at 20180426
def writeVaspFile(filename, atmnet, vornet, storeRadius = False, minRad = None, maxRad = None):
    if isinstance(filename, unicode):
         filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    if storeRadius:
        if minRad and maxRad:
            if not writeToVasp(c_filename, c_atmnet, c_vornet_ptr, True, minRad, maxRad):
                raise IOError
        else:
            if not writeToVasp(c_filename, c_atmnet, c_vornet_ptr, True):
                raise IOError
    else:
        if minRad and maxRad:
            if not writeToVasp(c_filename, c_atmnet, c_vornet_ptr, False, minRad, maxRad):
                raise IOError
        else:
            if not writeToVasp(c_filename, c_atmnet, c_vornet_ptr, False):
                raise IOError

# write to atomnetwork to .vasp file. Added at 20180827
def writeAtomNetVaspFile(filename, atmnet, storeRadius = False):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeAtmntToVasp(c_filename, c_atmnet, storeRadius):
        raise IOError


# write to .zvis file. In this file, all AtomNetwork and VoronoiNetwork information contained. Added at 20180704
def writeZVisFile(filename, rad_flag, atmnet, vornet):
    if isinstance(filename, unicode):
         filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    cdef vector[VOR_CELL] vcells
    cdef vector[BASIC_VCELL] bvcells
    if not performVoronoiDecomp(rad_flag, c_atmnet, c_vornet_ptr, &vcells, True, &bvcells):
        raise PerformVDError
    if not writeZVis(c_filename, &vcells, c_atmnet, c_vornet_ptr):
        raise IOError

# write to .net file. Added at 20190518
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