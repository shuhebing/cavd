"""
Cython file defining methods for AtomNetwork and VoronoiNetowrk 
declared in netstorage.pxd file. 
"""

__author__ = "Bharat Medasani"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__date__ = "Dec 12, 2013"

from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as inc

#from cavd.voronoicell cimport VorCell, BasicVCell, VOR_CELL, BASIC_VCELL
cimport cavd.netinfo
from cavd.voronoicell cimport VOR_CELL, BASIC_VCELL, VOR_FACE
from cavd.geometry cimport CPoint, Point

#Added at 20180606
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.string cimport string

#STUFF='Hi'

#Added at 20180807
#Customize an exception class
class PerformVDError(Exception):
    #print("Perform Voronoi Decompition failured!")
    pass

cdef class Atom:
    """
    Class to store the information about atom (or ion) in a structure.
    """
    def __cinit__(self):
        self.thisptr = new ATOM()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    property cart_coords:
        def __get__(self):
            coords = [self.thisptr.x, self.thisptr.y, self.thisptr.z]
            return coords
        def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.x = coords[0]
            self.thisptr.y = coords[1]
            self.thisptr.z = coords[2]

    property frac_coords:
        def __get__(self):
            coords = [self.thisptr.a_coord, self.thisptr.b_coord, self.thisptr.c_coord]
            return coords
        def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.a_coord = coords[0]
            self.thisptr.b_coord = coords[1]
            self.thisptr.c_coord = coords[2]

    property radius:
        def __get__(self): return self.thisptr.radius
        def __set__(self, radius): 
            print ("This value is not supposed to be modified")
            self.thisptr.radius = radius

    property atom_type:
        def __get__(self): return self.thisptr.atom_type.decode('utf-8')
        def __set__(self, atom_type): 
            print ("This value is not supposed to be modified")
            self.thisptr.atom_type = atom_type

    property label:
        def __get__(self): return self.thisptr.label.decode('utf-8')
        def __set__(self, label): 
            print ("This value is not supposed to be modified")
            self.thisptr.label = label

    property specialID:
        def __get__(self): return self.thisptr.specialID
        def __set__(self, specialID): 
            print ("This value is not supposed to be modified")
            self.thisptr.specialID = specialID

    property mass:
        def __get__(self): return self.thisptr.mass
        def __set__(self, mass): 
            print ("This value is not supposed to be modified")
            self.thisptr.mass = mass

    property charge:
        def __get__(self): return self.thisptr.charge
        def __set__(self, charge): 
            print ("This value is not supposed to be modified")
            self.thisptr.charge = charge

cdef class AtomNetwork:
    """
    Class to store and manipulate the input atom network.
    """
    #Cython wrapper for Zeo++ ATOM_NETWORK class.
    #Contains a pointer to ATOM_NETWORK and a flag denoting whether radius
    #for each atomic species is non-zero. 
    def __cinit__(self):
        self.thisptr = new ATOM_NETWORK()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr
    
    property lattice_para:
        def __get__(self):
            return [self.thisptr.a, self.thisptr.b, self.thisptr.c]

    property lattice_angle:
        def __get__(self):
            return [self.thisptr.alpha, self.thisptr.beta, self.thisptr.gamma]
    
    property lattice:
        def __get__(self):
            la = [self.thisptr.v_a.x, self.thisptr.v_a.y, self.thisptr.v_a.z]
            lb = [self.thisptr.v_b.x, self.thisptr.v_b.y, self.thisptr.v_b.z]
            lc = [self.thisptr.v_c.x, self.thisptr.v_c.y, self.thisptr.v_c.z]
            lattice = [la, lb, lc]
            return lattice

    property atoms_num:
        def __get__(self): return self.thisptr.no_atoms
        def __set__(self, atoms_num): 
            print ("This value is not supposed to be modified")
            self.thisptr.no_atoms = atoms_num

    property atoms:
        def __get__(self):
            atoms = []
            cdef vector[ATOM] c_atoms = self.thisptr.atoms
            for i in range(c_atoms.size()):
                atom_type = c_atoms[i].atom_type.decode('utf-8')
                atom_label = c_atoms[i].label.decode('utf-8')
                radius = c_atoms[i].radius
                #atom_coords = [c_atoms[i].a_coord,c_atoms[i].b_coord,c_atoms[i].c_coord]
                atom_coords = [c_atoms[i].x,c_atoms[i].y,c_atoms[i].z]
                atoms.append([atom_label, atom_type, radius, atom_coords])
            return atoms

    def copy(self):
        """
        Create a copy of the AtomNetwork instance
        """
        newatmnet = AtomNetwork()
        self.thisptr.copy(newatmnet.thisptr)
        newatmnet.rad_flag = self.rad_flag
        return newatmnet

    def relative_to_absolute(self, a, b, c):
        return [self.thisptr.abc_to_xyz(a, b, c).vals[0], self.thisptr.abc_to_xyz(a, b, c).vals[1], self.thisptr.abc_to_xyz(a, b, c).vals[2]]

    def absolute_to_relative(self, x, y, z):
        return [self.thisptr.xyz_to_abc(x, y, z).vals[0], self.thisptr.xyz_to_abc(x, y, z).vals[1], self.thisptr.xyz_to_abc(x, y, z).vals[2]]

    @classmethod
    def read_from_CIF(cls, filename, radii, rad_flag=True, rad_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a CIF file.
        Arguments:
            filename: 
                Input CIF file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, Zeo++ default values are used.
        Returns:
            Instance of AtomNetwork
        """
        #Calls Zeo++ readCIFFile function defined in networkio.cc.
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')
        
        cdef char* c_rad_file
        if rad_flag:
            if not rad_file:
                #edited at 20180526
                #cavd.netinfo.zeo_initializeRadTable()
                cavd.netinfo.zeo_initializeIonRadTable()
            else:       # rad_file is defined
                c_rad_file = rad_file
                cavd.netinfo.zeo_readIonRadTableFile(c_rad_file)

		#Added at 20180606
        cdef map[string, double] ionRadMap
        cdef string c_key
        cdef double c_value
        if radii:
            for key in radii:
                c_key = (<unicode>key).encode('utf8')
                c_value = radii[key]
                ionRadMap.insert(pair[string,double](c_key,c_value))       
            cavd.netinfo.zeo_readIonRadTable(ionRadMap)
        
        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readCIFFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_ARC(cls, filename, rad_flag=True, rad_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a ARC file.
        Arguments:
            filename: 
                Input ARC file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ readARCFile function defined in networkio.cc.
        cdef char* c_rad_file = rad_file
        if rad_flag:
            if not rad_file:
                cavd.netinfo.zeo_initializeRadTable()
            else:       # rad_file is defined
                c_rad_file = rad_file
                cavd.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readARCFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_CSSR(cls, filename, rad_flag=True, rad_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a CSSR file.
        Arguments:
            filename: 
                Input CSSR file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ readCSSRFile function defined in networkio.cc.
        cdef char* c_rad_file
        print (rad_flag, rad_file)
        if rad_flag:
            #if not rad_file:
            cavd.netinfo.zeo_initializeRadTable()
            if rad_file:       # rad_file is defined
                c_rad_file = rad_file
                cavd.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readCSSRFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_V1(cls, filename, rad_flag=True, rad_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a V1 file.
        Arguments:
            filename: 
                Input V1 file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
        Returns:
            Instance of AtomNetwork
        """
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ readV1File function defined in networkio.cc.
        cdef char* c_rad_file = rad_file
        if rad_flag:
            if not rad_file:
                cavd.netinfo.zeo_initializeRadTable()
            else:       # rad_file is defined
                cavd.netinfo.zeo_readRadTable(c_rad_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readV1File(c_filename, atmnet.thisptr, rad_flag):
            raise IOError
        atmnet.rad_flag = rad_flag
        return atmnet

    def write_to_CSSR(self, filename):
        """
        Writes the atom data in AtomNetwork to a CSSR file.
        Arguments:
            filename: 
                Output CSSR file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToCSSR function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToCSSR(c_filename, self.thisptr):
            raise IOError

    def write_to_CIF(self, filename):
        """
        Writes the atom data in AtomNetwork to a CIF file.
        Arguments:
            filename: 
                Output CIF file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToCIF function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToCIF(c_filename, self.thisptr):
            raise IOError

    def write_to_V1(self, filename):
        """
        Writes the atom data in AtomNetwork to a V1 file.
        Arguments:
            filename: 
                Output V1 file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToV1 function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToV1(c_filename, self.thisptr):
            raise IOError

    def write_to_XYZ(self, filename, supercell_flag, 
                     is_duplicate_perimeter_atoms):
        """
        Writes the atom data in AtomNetwork to an XYZ file.
        Arguments:
            filename: 
                Output XYZ file name.
            supercell_flag:
                Flag denoting whether to write 2x2x2 supercell.
            is_duplicate_perimeter_atoms:
                Flag denoting whether perimeter atoms need to be replicated.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToXYZ function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToXYZ(c_filename, self.thisptr, supercell_flag, 
                is_duplicate_perimeter_atoms):
            raise IOError

    def write_to_VTK(self, filename):
        """
        Writes the boundary of unit cell within the AtomNetwork to a VTK file.
        Arguments:
            filename: 
                Output VTK file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToVTK function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToVTK(c_filename, self.thisptr):
            raise IOError

    def write_to_MOPAC(self, filename, supercell_flag):
        """
        Writes the atom data in AtomNetwork to a .mop file.
        Arguments:
            filename: 
                Output MOPAC file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef char* c_filename = filename
        if not writeToMOPAC(c_filename, self.thisptr, supercell_flag):
             raise IOError

# write to atomnetwork to .vasp file. Added at 20180827
    def writeAtomNetVaspFile(self, filename, storeRadius = False):
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')
        cdef char* c_filename = filename
        if not writeAtmntToVasp(c_filename, self.thisptr, storeRadius):
            raise IOError
      

    def calculate_free_sphere_parameters(self, filename):
        """
        Computes the diameters of the largest included sphere, free sphere 
        and included sphere along free sphere path. 
        Arguments:
            filename:
                Name of file where the diameters are stored.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        sucess, vornet, edge_centers, face_centers = self.perform_voronoi_decomposition(False)
        cdef char* c_fname = filename
        vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
        calculateFreeSphereParameters(vornet_ptr, c_fname, True)
		#:q:q(vornet_ptr, c_fname, False)
      
      #Added at 20180420  
    def through_VorNet(self, filename):
        """
        Computes the diameters of the largest included sphere, free sphere 
        and included sphere along free sphere path. 
        Arguments:
        filename:
            Name of file where the diameters are stored.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')
        #if isinstance(migrantRad, unicode):
        #    migrantRad = (<unicode>migrantRad).encode('utf8')
    
        vornet, edge_centers, face_centers = self.perform_voronoi_decomposition(False)
        cdef char* c_fname = filename
        #cdef double c_migrantRad = migrantRad
        #Added at 20180530
        cdef double* c_Ri_ptr
        cdef double* c_Rf_ptr
        cdef double* c_Rif_ptr
        cdef double c_Ri,c_Rf,c_Rif
        vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
        if throughVorNet(vornet_ptr, c_fname, &c_Ri, &c_Rf, &c_Rif):
            #return True
            #edited at 20180530
            return c_Ri,c_Rf,c_Rif

    def perform_voronoi_decomposition(self, saveVorCells=True):
        """
        Performs weighted voronoi decomposition of atoms in the AtomNetwork 
        to analyze void space and generate voronoi nodes, edges and faces.
        Arguments:
            saveVorCells (optional): 
                Flag to denote whether to save the VorCells.
                Reserved for future use, so ignore this.
        Returns:
            Instance of VoronoiNetwork
        """
        #Calls Zeo++ performVoronoiDecomp function defined in network.cc.
        vornet = VoronoiNetwork()  
        cdef vector[VOR_CELL] vcells
        cdef vector[BASIC_VCELL] bvcells
        #print self.rad_flag
        if not performVoronoiDecomp(self.rad_flag, self.thisptr, 
                vornet.thisptr, &vcells, saveVorCells, &bvcells):
            #edited at 20180604
            #Add a compute flag
            #success = False
            raise PerformVDError
        #else:
            #success = True
        cdef int N

        # Get the edge centers
        edge_centers = []
        cdef vector[VOR_EDGE] vedges = vornet.thisptr.edges
        cdef vector[VOR_NODE] vnodes = vornet.thisptr.nodes
        for i in range(vedges.size()):
            edge_orig =  vedges[i].origin
            edge_end =  vedges[i].ending
            o_vnode = vnodes[edge_orig]
            e_vnode = vnodes[edge_end]
            edge_center = (o_vnode.x + e_vnode.x, \
                           o_vnode.y + e_vnode.y, \
                           o_vnode.z + e_vnode.z)
            edge_center = tuple(x/2 for x in edge_center)
            if edge_center not in edge_centers:
                edge_centers.append(edge_center)



        # Get the vorcells and obtain the face centers
        face_centers = []
        cdef vector[VOR_FACE] vfaces
        cdef vector[CPoint] vertices
        cdef CPoint* cpoint_ptr 
        #cdef map[int, int] id_maps
        cdef vector[int] node_ids
        face_node_ids = set()
        for i in range(vcells.size()):
            vfaces = vcells[i].faces
            for j in range(vfaces.size()):
                node_ids = vfaces[j].node_ids
                node_id_list = []
                for k in range(node_ids.size()):
                    node_id_list.append(node_ids[k])
                node_id_set = frozenset(node_id_list)
                if not node_id_set in face_node_ids:
                    face_node_ids.add(node_id_set)
                    centroid = Point()
                    cpoint_ptr = (<Point?>centroid).thisptr
                    vertices = vfaces[j].vertices
                    for k in range(vertices.size()):
                        centroid.x = centroid.x + vertices[k].vals[0]
                        centroid.y = centroid.y + vertices[k].vals[1]
                        centroid.z = centroid.z + vertices[k].vals[2]
                    centroid.x = centroid.x/vertices.size()
                    centroid.y = centroid.y/vertices.size()
                    centroid.z = centroid.z/vertices.size()
                    face_centers.append(centroid)

        # Convert the Zeo++ Point objects in (x,y,z) tuple objects
        fcs = []
        for center in face_centers:
            cntr = (center.x,center.y,center.z)
            fcs.append(cntr)

        #bvcelllist = []
        # define copy methods for BASIC_VCELL and VOR_CELL methods
        #for i in range(vcells.size()):
        #    pass
            #vorcell = VorCell()
            #vorcell.thisptr = &(vcells[i])
            #vorcelllist.append(vcells[i])
        #for i in range(bvcells.size()):
        #    pass
            #basicvcell = BasicVCell()
            #basicvcell.thisptr = &(bvcells[i])
            #bvcelllist.append(bvcells[i])
        return vornet, edge_centers, fcs

cdef class VoronoiNode:
    """
    Class to store the voronoi nodes with coordinates and radius
    """
    def __cinit__(self):
        self.thisptr = new VOR_NODE()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    property coords:
        def __get__(self):
            coords = [self.thisptr.x, self.thisptr.y, self.thisptr.z]
            return coords
        def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.x = coords[0]
            self.thisptr.y = coords[1]
            self.thisptr.z = coords[2]
    
    property label:
        def __get__(self): return self.thisptr.label
        def __set__(self, label):
            self.thisptr.label = label

    property radius:
        def __get__(self): return self.thisptr.rad_stat_sphere
        def __set__(self, rad): 
            print ("This value is not supposed to be modified")
            self.thisptr.rad_stat_sphere = rad

cdef class VoronoiEdge:
    """
    Class to store the voronoi edges with some atrribute
    """
    def __cinit__(self):
        self.thisptr = new VOR_EDGE()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr
    
    property origin:
        def __get__(self): return self.thisptr.origin
        def __set__(self, origin): 
            print ("This value is not supposed to be modified")
            self.thisptr.origin = origin

    property ending:
        def __get__(self): return self.thisptr.ending
        def __set__(self, ending): 
            print ("This value is not supposed to be modified")
            self.thisptr.ending = ending

    property radius:
        def __get__(self): return self.thisptr.rad_moving_sphere
        def __set__(self, rad): 
            print ("This value is not supposed to be modified")
            self.thisptr.rad_moving_sphere = rad

    property leng:
        def __get__(self): return self.thisptr.length
        def __set__(self, length): 
            print ("This value is not supposed to be modified")
            self.thisptr.length = length

    property delta_uc:
        def __get__(self):
            delta_uc = [self.thisptr.delta_uc_x, self.thisptr.delta_uc_y, self.thisptr.delta_uc_z]
            return delta_uc
        def __set__(self, delta_uc):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.delta_uc_x = delta_uc[0]
            self.thisptr.delta_uc_y = delta_uc[1]
            self.thisptr.delta_uc_z = delta_uc[2]

    property bot_coords:
        def __get__(self):
            bot_coords = [self.thisptr.bottleneck_x, self.thisptr.bottleneck_y, self.thisptr.bottleneck_z]
            return bot_coords
        def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print ("This value is not supposed to be modified")
            self.thisptr.bottleneck_x = coords[0]
            self.thisptr.bottleneck_y = coords[1]
            self.thisptr.bottleneck_z = coords[2]

cdef class VoronoiNetwork:
    """
    Class to store the Voronoi network generated from Voronoi decomposition
    of atom network.
    """
    #Cython wrapper for Zeo++ VORONOI_NETWORK class.
    #Contains a pointer to ATOM_NETWORK and a flag denoting whether radisu
    #for each atomic species is non-zero. 
    def __cinit__(self):
        self.thisptr = new VORONOI_NETWORK()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    def size(self):
        return self.thisptr.nodes.size()

    property nodes:
        def __get__(self):
            nodes = []
            cdef vector[VOR_NODE] c_nodes = self.thisptr.nodes
            for i in range(c_nodes.size()):
                node_coords = [c_nodes[i].x,c_nodes[i].y,c_nodes[i].z]
                node_radius = c_nodes[i].rad_stat_sphere
                node_label = c_nodes[i].label
                nodes.append([i, node_label, node_coords, node_radius])
            return nodes

    property edges:
        def __get__(self):
            edges = []
            cdef vector[VOR_EDGE] c_edges = self.thisptr.edges
            for i in range(c_edges.size()):
                edge_origin = c_edges[i].origin
                edge_ending = c_edges[i].ending
                edge_radius = c_edges[i].rad_moving_sphere
                edge_length = c_edges[i].length
                edge_boltpos = [c_edges[i].bottleneck_x,c_edges[i].bottleneck_y,c_edges[i].bottleneck_z]
                edges.append([edge_origin, edge_ending, edge_radius, edge_length, edge_boltpos])
            return edges

    def prune(self, radius):
        """
        Removes the edges that do not allow a sphere to pass.
        Arguments:
            radius:
                Radius of the sphere
        Returns:
            Instance of VoronoiNetwork with edges pruned.
        """
        cdef VORONOI_NETWORK newcvornet = self.thisptr.prune(radius)
        newvornet = VoronoiNetwork()
        newvornet.thisptr = &newcvornet
        return newvornet

    def analyze_writeto_XYZ(self, name, double probeRad, atmnet, 
            int shift_x=0, int shift_y=0, int shift_z=0):
        """
        Create diagrams of 1) Voronoi network and 2) accessible Voronoi 
        network, and write the diagrams in VTK files and the Voronoi 
        networks in XYZ files. Useful for visualizing the Voronoi network.
        Args:
            name:
                Name to be used for output files.
            probeRad:
                Radius of the probe.
            atmnet:
                cavd.netstorage.AtomNetwork
            shift_x (default=0):
                Shift the accessible Voronoi network along x-axis
            shift_y (default=0):
                Shift the accessible Voronoi network along y-axis
            shift_z (default=0):
                Shift the accessible Voronoi network along z-axis
        """
        if isinstance(name, unicode):
            name = (<unicode>name).encode('utf8')

        cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
        cdef char* cname = name
        visVoro(name, probeRad, shift_x, shift_y, shift_z, self.thisptr, 
                c_atmnetptr)

    def write_to_XYZ(self, filename, double cutoff_radius=0):
        """
        Write only voronoi node information to XYZ file.
        Args:
            filename:
                string
                Name of file to which voronoi node info is written.
            cutoff_radius:
                float
                Threshold radius (default=0)
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef char* c_filename = filename
        if not writeVornetToXYZ(c_filename, self.thisptr, 
                cutoff_radius):
            raise ValueError

    @classmethod
    def perform_voronoi_decomposition(cls, atmnet, saveVorCells=False):
        """
        Performs weighted voronoi decomposition of atoms in the AtomNetwork 
        to analyze void space and generate voronoi nodes, edges and faces.
        Arguments:
            saveVorCells (optional): 
                Flag to denote whether to save the VorCells.
                Reserved for future use, so ignore this.
        Returns:
            Instance of VoronoiNetwork
        """
        #Calls Zeo++ performVoronoiDecomp function defined in network.cc.
        vornet = VoronoiNetwork()  
        cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
        cdef vector[VOR_CELL] vcells
        cdef vector[BASIC_VCELL] bvcells
        #print self.rad_flag
        if not performVoronoiDecomp(atmnet.rad_flag, c_atmnetptr, 
                vornet.thisptr, &vcells, saveVorCells, &bvcells):
            #raise ValueError # Change it to appropriate error
            raise PerformVDError
        #cdef int N
        #vorcelllist = []
        #bvcelllist = []
        # define copy methods for BASIC_VCELL and VOR_CELL methods
        #for i in range(vcells.size()):
        #    pass
            #vorcell = VorCell()
            #vorcell.thisptr = &(vcells[i])
            #vorcelllist.append(vcells[i])
        #for i in range(bvcells.size()):
        #    pass
            #basicvcell = BasicVCell()
            #basicvcell.thisptr = &(bvcells[i])
            #bvcelllist.append(bvcells[i])
        return vornet
    
    def parse_symmetry(self,symm_label):
        cdef vector[int] c_symm_label
        for i in range(len(symm_label)):
            c_symm_label.push_back(symm_label[i])
        parseNetworkSymmetry(c_symm_label, self.thisptr)
        return self


def substitute_atoms(atmnet, substituteSeed, radialFlag):
    """
    Attempt to substitute every other Si atom with Al atom.
    AtomNetwork may only consist of Si and O atoms, where each Si atom 
    must be bonded to exactly 4 oxygen atoms and each oxygen atom must 
    be bonded to exactly 2 Si atoms. Raises exception if the substitution
    is not successful. 
    Args:
        atmnet:
            cavd.netstorage.AtomNetwork
        substiuteSeed:
            Boolean flag to specify whether the seeded Si atom is 
            substituted or not. Since only 2 configurations are possible 
            if the structure is consistent, changing this parameter enables 
            generation of all configurations. 
        radialFlag:
            Boolean flag to specify whether atomic sizes are to be used.
    Returns:
        If successful, returns AtomNetwork instance with Si replaced with Al
        and the number of substitutions. 
    """
    cdef int substitutionNo[1]
    atmnet_copy = AtomNetwork()
    c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    if not c_substituteAtoms(c_atmnet_ptr, atmnet_copy.thisptr, substituteSeed,
            substitutionNo, radialFlag):
        raise ValueError
    subNo = substitutionNo[0]
    return atmnet_copy, subNo

def connection_values(filename, vornet):
    """
	Computes the Radius of the largest included sphere, free sphere 
	and included sphere along free sphere path. 
	Arguments:
	filename:
		Name of file where the diameters are stored.
    """
	
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')

    vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    cdef char* c_fname = filename
    cdef double c_Ri,c_Rf,c_Rif
    throughVorNet(vornet_ptr, c_fname, &c_Ri, &c_Rf, &c_Rif)
    return c_Ri,c_Rf,c_Rif	

def connection_values_list(filename, vornet):
    conn_values = []
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    cdef char* c_fname = filename
    cdef vector[double] values 
    calculateConnParameters(vornet_ptr, c_fname, &values)
    conn_values = []
    for i in range(values.size()):
        conn_values.append(values[i])
    return conn_values

