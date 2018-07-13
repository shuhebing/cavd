from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "../../zeo++/string_additions.h":
    cdef int strCmpList(vector[string] list, string str)
 
