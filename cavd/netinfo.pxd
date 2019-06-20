# distutils: language = c++
# distutils: sources = ../networkinfo.cc

from libcpp.map cimport map
from libcpp.string cimport string

cdef extern from "../basic_lib/Zeo++/networkinfo.h":
    cdef void zeo_initializeRadTable "initializeRadTable"()

    cdef void zeo_initializeCovRadTable "initializeCovRadTable"()

    cdef void zeo_initializeMassTable "initializeMassTable"()

    cdef void zeo_initializeAtomCharacterTable "initializeAtomCharacterTable"()

    cdef void zeo_initializeAtomicNumberTable "initializeAtomicNumberTable"()

    cdef void zeo_readRadTable "readRadTable"(char *filename)

    cdef void zeo_readMassTable "readMassTable"(char *filename)

    cdef double zeo_lookupRadius "lookupRadius"(string atomType, bint radial)

    cdef double zeo_lookupCovRadius "lookupCovRadius"(string atomType)

    cdef double zeo_lookupMass "lookupMass"(string atomType)

    cdef int zeo_lookupAtomicNumber "lookupAtomicNumber"(string atomType)

    cdef bint zeo_isMetal "isMetal"(string atomType)
    
    # Added at 20180606
    cdef void zeo_readIonRadTable "readIonRadTable"(map[string,double] radMap)
    cdef double zeo_lookupIonRadius "lookupIonRadius"(string atomType, bint radial)
    
    #Added at 20180627
    cdef void zeo_initializeIonRadTable "initializeIonRadTable"()
    cdef void zeo_readIonRadTableFile "readIonRadTableFile"(char *filename)


