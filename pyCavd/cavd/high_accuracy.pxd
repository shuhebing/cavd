from libcpp.string cimport string
from cavd.netstorage cimport ATOM_NETWORK 

cdef extern from "../../zeo++/sphere_approx.h":
    cdef void setupHighAccuracyAtomNetwork(ATOM_NETWORK *atmnet, 
            string AccSetting)
