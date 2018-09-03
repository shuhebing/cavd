# distutils: language = c++
# distutils: sources = ../networkinfo.cc
"""
Wrapper functions to Zeo++ atomic definitons and related functions
"""
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.string cimport string

#Python definitions for the cdefinitions in .pxd file
def initializeRadTable():
    """
    Populate the atomic radius table with Zeo++ default values
    """
    zeo_initializeRadTable()

def initializeCovRadTable():
    """
    Populate the covalent tradius table with Zeo++ default values
    """
    zeo_initializeCovRadTable()

def initializeMassTable():
    """
    Populate the atomic mass table with Zeo++ default values
    """
    zeo_initializeMassTable()

def initializeAtomCharacterTable():
    """
    Populate the Atom symbol table with Zeo++ default values
    """
    zeo_initializeAtomCharacterTable()

def initializeAtomicNumberTable():
    """
    Populate the atomic number table with Zeo++ default values
    """
    zeo_initializeAtomicNumberTable()

def readRadTable(filename):
    """
    Read atomic radii values from input file and replace the default values
    """
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    zeo_readRadTable(c_filename)

def readMassTable(filename):
    """
    Read atomic mass values from input file and replace the default values
    """
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    zeo_readMassTable(c_filename)

def lookupRadius(element):
    """"
    Args:
        element:
            Element name in conventional shorthand 
            Ex: Al for aluminum 
                Si for silicon 
    Returns:
        radius of the input element
    """
    radius = zeo_lookupRadius(element, True)
    return radius
    
def lookupCovRadius(element):
    return zeo_lookupCovRadius(element)

def lookupMass(element):
    return zeo_lookupMass(element)

def lookupAtomicNumber(element):
    return zeo_lookupAtomicNumber(element)

def isMetal(element):
    return zeo_isMetal(element)

#Added at 20180420
#def initializeGoldschmidtIonRadTable():
    """
    Populate the Goldschmidt Ion radius table with Zeo++ default values
    """
#    zeo_initializeGoldschmidtIonRadTable()

#def lookupGoldschmidtIonRadius(element):
    """"
    Args:s
        element:
            Element name in conventional shorthand 
            Ex: Al for aluminum 
                Si for silicon 
    Returns:
        radius of the input element
    """
    #added at 20180604
#    if isinstance(element, unicode):
#        element = (<unicode>element).encode('utf8')
#    cdef string c_element = element
#    radius = zeo_lookupGoldschmidtIonRadius(c_element, True)
#    return radius    

#Added at 20180606
def readIonRadTable(ionicRadDic):
    """
    Read Ionic radius values from input Dictionary
    """
    cdef map[string, double] ionRadMap
    cdef string c_key
    cdef double c_value
    for key in ionicRadDic:
        c_key = (<unicode>key).encode('utf8')
        c_value = ionicRadDic[key]
        ionRadMap.insert(pair[string,double](c_key,c_value))       
    zeo_readIonRadTable(ionRadMap)

def lookupIonRadius(element):
    """"
    Args:s
        element:
            Element name in conventional shorthand 
            Ex: Al for aluminum 
                Si for silicon 
    Returns:
        radius of the input element
    """
    #added at 20180606
    if isinstance(element, unicode):
        element = (<unicode>element).encode('utf8')
    cdef string c_element = element
    radius = zeo_lookupIonRadius(c_element, True)
    return radius    
	
#Added at 20180627
def initializeIonRadTable():
    """
    Populate the Goldschmidt Ion radius table with Zeo++ default values
    """
    zeo_initializeIonRadTable()

def readIonRadTableFile(filename):
    """
    Read Ionic radii values from input file and replace the default values
    """
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    zeo_readIonRadTableFile(c_filename)
