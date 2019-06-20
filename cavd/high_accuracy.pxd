from libcpp.string cimport string
from cavd.netstorage cimport ATOM_NETWORK 

cdef extern from "../basic_lib/Zeo++/sphere_approx.h":
    cdef void setupHighAccuracyAtomNetwork(ATOM_NETWORK *atmnet, 
            string AccSetting)
