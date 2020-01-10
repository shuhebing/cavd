'''

Symmetry analysis.

Created on 2018.7.2

@author: YeAnjiang
'''
import re
import spglib
import numpy as np
import pandas as pd
import ase.spacegroup as spg
from scipy.spatial.ckdtree import cKDTree
from pymatgen.util.io_utils import clean_lines
from pymatgen.core.structure import Structure
from monty.io import zopen


class Poscar_new():
    def __init__(self, atomic_symbols, coords, lattice, comment=None, selective_dynamics=None,
                 true_names=True, velocities=None, predictor_corrector=None,
                 predictor_corrector_preamble=None):
        
        self.atomic_symbols = atomic_symbols
        self.coords = coords
        self.lattice = lattice
        self.velocities = velocities
        self.true_names = true_names
        self.selective_dynamics = selective_dynamics
        self.predictor_corrector = predictor_corrector
        self.comment = structure.formula if comment is None else comment
        self.predictor_corrector_preamble = predictor_corrector_preamble
        self.temperature = -1    

    def from_string(data, default_names=None, read_velocities=True):
        """
        Reads a Poscar from a string.

        The code will try its best to determine the elements in the POSCAR in
        the following order:
        1. If default_names are supplied and valid, it will use those. Usually,
        default names comes from an external source, such as a POTCAR in the
        same directory.
        2. If there are no valid default names but the input file is Vasp5-like
        and contains element symbols in the 6th line, the code will use that.
        3. Failing (2), the code will check if a symbol is provided at the end
        of each coordinate.

        If all else fails, the code will just assign the first n elements in
        increasing atomic number, where n is the number of species, to the
        Poscar. For example, H, He, Li, ....  This will ensure at least a
        unique element is assigned to each site and any analysis that does not
        require specific elemental properties should work fine.

        Args:
            data (str): String containing Poscar data.
            default_names ([str]): Default symbols for the POSCAR file,
                usually coming from a POTCAR in the same directory.
            read_velocities (bool): Whether to read or not velocities if they
                are present in the POSCAR. Default is True.

        Returns:
            Poscar object.
        """
        # "^\s*$" doesn't match lines with no whitespace
        chunks = re.split(r"\n\s*\n", data.rstrip(), flags=re.MULTILINE)
        try:
            if chunks[0] == "":
                chunks.pop(0)
                chunks[0] = "\n" + chunks[0]
        except IndexError:
            raise ValueError("Empty POSCAR")

        # Parse positions
        lines = tuple(clean_lines(chunks[0].split("\n"), False))
        comment = lines[0]
        scale = float(lines[1])
        lattice = np.array([[float(i) for i in line.split()]
                            for line in lines[2:5]])
        if scale < 0:
            # In vasp, a negative scale factor is treated as a volume. We need
            # to translate this to a proper lattice vector scaling.
            vol = abs(det(lattice))
            lattice *= (-scale / vol) ** (1 / 3)
        else:
            lattice *= scale

        vasp5_symbols = False
        try:
            natoms = [int(i) for i in lines[5].split()]
            ipos = 6
        except ValueError:
            vasp5_symbols = True
            symbols = lines[5].split()

            """
            Atoms and number of atoms in POSCAR written with vasp appear on 
            multiple lines when atoms of the same type are not grouped together 
            and more than 20 groups are then defined ...
            
            Example :
            
            Cr16 Fe35 Ni2
               1.00000000000000
                 8.5415010000000002   -0.0077670000000000   -0.0007960000000000
                -0.0077730000000000    8.5224019999999996    0.0105580000000000
                -0.0007970000000000    0.0105720000000000    8.5356889999999996
               Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Ni   Fe   Cr   Fe   Cr
               Fe   Ni   Fe   Cr   Fe
                 1     1     2     4     2     1     1     1     2     1     1     1     4     1     1     1     5     3     6     1
                 2     1     3     2     5
            Direct
              ...
            """
            nlines_symbols = 1
            for nlines_symbols in range(1, 11):
                try:
                    int(lines[5+nlines_symbols].split()[0])
                    break
                except ValueError:
                    pass
            for iline_symbols in range(6, 5+nlines_symbols):
                symbols.extend(lines[iline_symbols].split())
            natoms = []
            iline_natoms_start = 5+nlines_symbols
            for iline_natoms in range(iline_natoms_start,
                                      iline_natoms_start+nlines_symbols):
                natoms.extend([int(i) for i in lines[iline_natoms].split()])
            atomic_symbols = list()
            for i in range(len(natoms)):
                atomic_symbols.extend([symbols[i]] * natoms[i])
            ipos = 5+2*nlines_symbols

        postype = lines[ipos].split()[0]

        sdynamics = False
        # Selective dynamics
        if postype[0] in "sS":
            sdynamics = True
            ipos += 1
            postype = lines[ipos].split()[0]

        cart = postype[0] in "cCkK"
        nsites = sum(natoms)

        # If default_names is specified (usually coming from a POTCAR), use
        # them. This is in line with Vasp"s parsing order that the POTCAR
        # specified is the default used.
        if default_names:
            try:
                atomic_symbols = []
                for i in range(len(natoms)):
                    atomic_symbols.extend([default_names[i]] * natoms[i])
                vasp5_symbols = True
            except IndexError:
                pass

        if not vasp5_symbols:
            ind = 3 if not sdynamics else 6
            try:
                # Check if names are appended at the end of the coordinates.
                atomic_symbols = [l.split()[ind]
                                  for l in lines[ipos + 1:ipos + 1 + nsites]]
                # Ensure symbols are valid elements
                if not all([Element.is_valid_symbol(sym)
                            for sym in atomic_symbols]):
                    raise ValueError("Non-valid symbols detected.")
                vasp5_symbols = True
            except (ValueError, IndexError):
                # Defaulting to false names.
                atomic_symbols = []
                for i in range(len(natoms)):
                    sym = Element.from_Z(i + 1).symbol
                    atomic_symbols.extend([sym] * natoms[i])
                warnings.warn("Elements in POSCAR cannot be determined. "
                              "Defaulting to false names %s." %
                              " ".join(atomic_symbols))

        # read the atomic coordinates
        coords = []
        selective_dynamics = list() if sdynamics else None
        for i in range(nsites):
            toks = lines[ipos + 1 + i].split()
            crd_scale = scale if cart else 1
            coords.append([float(j) * crd_scale for j in toks[:3]])
            if sdynamics:
                selective_dynamics.append([tok.upper()[0] == "T"
                                           for tok in toks[3:6]])

        if read_velocities:
            # Parse velocities if any
            velocities = []
            if len(chunks) > 1:
                for line in chunks[1].strip().split("\n"):
                    velocities.append([float(tok) for tok in line.split()])

            # Parse the predictor-corrector data
            predictor_corrector = []
            predictor_corrector_preamble = None

            if len(chunks) > 2:
                lines = chunks[2].strip().split("\n")
                # There are 3 sets of 3xN Predictor corrector parameters
                # So can't be stored as a single set of "site_property"

                # First line in chunk is a key in CONTCAR
                # Second line is POTIM
                # Third line is the thermostat parameters
                predictor_corrector_preamble = (lines[0] + "\n" + lines[1]
                                                + "\n" + lines[2])
                # Rest is three sets of parameters, each set contains
                # x, y, z predictor-corrector parameters for every atom in orde
                lines = lines[3:]
                for st in range(nsites):
                    d1 = [float(tok) for tok in lines[st].split()]
                    d2 = [float(tok) for tok in lines[st+nsites].split()]
                    d3 = [float(tok) for tok in lines[st+2*nsites].split()]
                    predictor_corrector.append([d1,d2,d3])
        else:
            velocities = None
            predictor_corrector = None
            predictor_corrector_preamble = None

        return Poscar_new(atomic_symbols, coords, lattice, comment, selective_dynamics, vasp5_symbols,velocities=velocities, predictor_corrector=predictor_corrector,predictor_corrector_preamble=predictor_corrector_preamble)

# Using spglib to analyze the space group in .vasp file
def get_symmetry_spglib(filename, symprec=0.00001):
    with zopen(filename, "rt") as f:
        contents = f.read()
    poscar = Poscar_new.from_string(contents, False, read_velocities=False)
    positions = poscar.coords
    lattice = poscar.lattice
    atomic_symbols = poscar.atomic_symbols
    numbers = [] 
    a = ""
    j = 0
    for i in atomic_symbols:
        if i != a:
            a = i
            j = j+1
        numbers.append(j)

    cell = (lattice, positions, numbers)
    dataset = spglib.get_symmetry_dataset(cell, symprec)
    print("space group: ", dataset['international'],dataset['number'])
    print("rotations: ", dataset['rotations'])
    print("translations: ", dataset['translations'])
    print("equivalent atoms: ", dataset['equivalent_atoms'])
    print("sites wyckoffs: ", dataset['wyckoffs'])
    sym_independ = np.unique(dataset['equivalent_atoms'])
    print("independent atoms: ", sym_independ)
    for i in sym_independ:
        print("coordinates of independent atoms")
        print(positions[i])

"""
    Return the symmetry number of input positions in unitcell.
    Input:
        lattice
        positions
    Output:
        symmetry number or zero.
"""
def get_symnum_sites(lattice, positions, symprec=0.01, angle_tolerance=5):
    numbers = [1,]*len(positions)
    cell = (lattice, positions, numbers)
    dataset = spglib.get_symmetry_dataset(cell, symprec, angle_tolerance)
    if dataset:
        return dataset['number']
    else:
        return 0

"""
Function:
    Get the symmetry equivalent sites for input site.
"""
def get_equivalent_VorNodes(pos, sitesym):
    rot, trans = spg.spacegroup.parse_sitesym(sitesym)
    sympos = np.dot(rot, pos) + trans
    return sympos

"""
    Analyzing the symmetry of Voronoi Network by spglib.
"""
def get_equivalent_vornet(vornet, symprec=1e-5, angle_tolerance=5):
    positions = []
    lattice = vornet.lattice
    for i in vornet.nodes:
        positions.append(i[2])
    numbers = [1,]*len(vornet.nodes)
    cell = (lattice, positions, numbers)

    dataset = spglib.get_symmetry_dataset(cell, symprec, angle_tolerance)
    print(dataset['number'])
    voids = []
    if dataset:
        symm_label = dataset['equivalent_atoms']
        vornet_uni_symm = vornet.parse_symmetry(symm_label)
        sym_independ = np.unique(dataset['equivalent_atoms'])
        print("symm_label", symm_label)
        print("sym_independ", sym_independ)
        print("The number of symmetry distinct voids: ",len(sym_independ))
        for i in sym_independ:
            voids.append(positions[i])
        return vornet_uni_symm,voids
    else:
        return vornet,voids

"""
    Analyzing the symmetry of Voronoi Network by ourself code.
"""
def get_labeled_vornet(vornet, sitesym, symprec=1e-3):
    positions = []
    for i in vornet.nodes:
        positions.append(i[2])

    tags,tagdis = tag_sites(sitesym, positions,symprec)
    voids = []

    vornet_uni_symm = vornet.parse_symmetry(tags)
    sym_independ = np.unique(tags)
    # print("The number of symmetry distinct voids: ",len(sym_independ))
    for i in sym_independ:
         voids.append(positions[i])
    return vornet_uni_symm,voids

"""
    Get unique sites from .vasp file.
"""
def get_unique_sites(filename, sitesym, symprec=1e-3):
    with zopen(filename, "rt") as f:
        contents = f.read()
    poscar = Poscar_new.from_string(contents, False, read_velocities=False)
    positions = poscar.coords
    tags,tagdis = tag_sites(sitesym, positions, symprec)
    print(tags)
    print(tagdis)
    print(np.unique(tags))
    frame = pd.Series(tags)
    print(frame.value_counts())   

"""
    Get symmetry equivalent sites of provided scaled_positions 
    based on provided symmetry operations. This function will 
    return a mapping table of sites to symmetrically independent 
    sites.This is used to find symmetrically equivalent atoms. 
    The numbers contained are the indices of sites starting from 0, 
    i.e., the first atom is numbered as 0, and then 1, 2, 3, … 
    np.unique(equivalent_sites) gives representative symmetrically 
    independent sites.
"""
def tag_sites(sitesym, scaled_positions, symprec=1e-3):
    scaled = np.around(np.array(scaled_positions, ndmin=2),8)
    scaled %= 1.0
    scaled %= 1.0
    np.set_printoptions(suppress=True)
    tags = -np.ones((len(scaled), ), dtype=int)
    tagdis = 100*np.ones((len(scaled), ), dtype=float)
    rot, trans = spg.spacegroup.parse_sitesym(sitesym)

    siteskdTree = cKDTree(scaled)
    for i in range(len(scaled)):
        if tags[i] == -1:
            curpos = scaled[i]
            sympos = np.dot(rot, curpos) + trans
            sympos %= 1.0
            sympos %= 1.0
            sympos = np.unique(np.around(sympos,8), axis=0)
            min_dis,min_ids = siteskdTree.query(sympos,k=1)
            # print(i,len(min_dis))
            # print(min_dis)

            select = min_dis <= symprec
            select_ids = min_ids[select]
            tags[select_ids] = i
            tagdis[select_ids] = min_dis[select]
    return tags,tagdis

from cavd.local_environment import CifParser_new
