# distutils: language = c++
# distutils: sources = ../libs/Zeo++/networkio.cc

from cavd.netstorage cimport ATOM_NETWORK, VORONOI_NETWORK
from libcpp.string cimport string

#Added at 20180704
from libcpp.vector cimport vector
from cavd.voronoicell cimport VOR_CELL, BASIC_VCELL
from cavd.channel cimport CHANNEL

cdef extern from '../libs/Zeo++/networkio.h':
    cdef void parseFilename(const char* filename, char* name, char* extension)

    cdef bint checkInputFile(char* filename)

    cdef bint readCIFFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readARCFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readCUCFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readCSSRFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readV1File(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint writeToCSSR(char *filename, ATOM_NETWORK *cell)

    cdef bint writeToCIF(char *filename, ATOM_NETWORK *cell)

    cdef bint writeToV1(char * filename, ATOM_NETWORK *cell)

    cdef bint writeToNt2(char *filename, VORONOI_NETWORK *vornet, double minRad)

    cdef bint writeToNt2(char *filename, VORONOI_NETWORK *vornet)

    cdef bint writeToXYZ(char *filename, ATOM_NETWORK *cell, bint is_supercell,
                         bint is_duplicate_perimeter_atoms)

    cdef bint writeToVTK(char *filename, ATOM_NETWORK *cell)

    cdef bint writeToMOPAC(char *filename, ATOM_NETWORK *cell, bint is_supercell)

    cdef bint writeVornetToXYZ "writeToXYZ"(char *filename, VORONOI_NETWORK*, double)

# remove migrant ion added at 20180408
    cdef bint readRemoveMigrantCif(char *filename, ATOM_NETWORK *cell, const char *migrant, bint radial)
    
# writeToVasp
    #cdef bint writeToVasp(char *filename, ATOM_NETWORK *cell, VORONOI_NETWORK *vornet, bint storeRadius, double minRad)
    #edited at 20180530
    cdef bint writeToVasp(char *filename, ATOM_NETWORK *cell, VORONOI_NETWORK *vornet, double minRad, double maxRad)
    cdef bint writeToVasp(char *filename, ATOM_NETWORK *cell, VORONOI_NETWORK *vornet)
    cdef bint writeAtmntToVasp(char *filename, ATOM_NETWORK *cell)

    #add at 20190518
    cdef bint writeToNET(char *filename, ATOM_NETWORK *cell, VORONOI_NETWORK *vornet, double minRad, double maxRad)
    cdef bint writeToNET(char *filename, ATOM_NETWORK *cell, VORONOI_NETWORK *vornet)

# At present  the return value of performVoronoiDecomp is void*
# Compile it after void* is changed to bool in the original source file
cdef extern from "../libs/Zeo++/network.h":
    cdef bint performVoronoiDecomp(bint, ATOM_NETWORK*, VORONOI_NETWORK*, 
            vector[VOR_CELL]*, bint, vector[BASIC_VCELL]*)