'''
Created on 2018年7月12日

@author: YeAnjiang
'''
import re
from pymatgen.util.io_utils import clean_lines
from pymatgen.core.structure import Structure
from monty.io import zopen
import numpy as np
from pymatgen.io.vasp import Poscar
import spglib
from cavd.netstorage import AtomNetwork
from cavd import LocalEnvirCom,writeVaspFile
from cavd.netio import getRemoveMigrantFilename
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


def get_Symmetry(filename):
    with zopen(filename, "rt") as f:
        contents = f.read()
    poscar = Poscar_new.from_string(contents, False, read_velocities=False)
    positions = poscar.coords
    #positions = []
    print(positions)
    lattice = poscar.lattice
    print(lattice)
    #for i in s2.sites:
    #    positions.append(i._fcoords)
    atomic_symbols = poscar.atomic_symbols
    #print(atomic_symbols)
    
    numbers = [] 
    a = ""
    j = 0
    for i in atomic_symbols:
        if i != a:
            a = i
            j= j+1
        numbers.append(j)

    #numbers = s2.atomic_numbers
    
    print(numbers)
    print(len(positions))
    print(len(numbers))
    cell = (lattice, positions, numbers)

    spacegroup = spglib.get_spacegroup(cell, symprec=0.01, angle_tolerance=5)
    print(spacegroup)  
    #symmetry = spglib.get_symmetry(cell, symprec=1e-5)
    #print(symmetry)
    dataset = spglib.get_symmetry_dataset(cell, symprec=0.01, angle_tolerance=5)
    print(len(dataset['equivalent_atoms']))
    #print(dataset['rotations'])
    #print(dataset['translations'])
    print(dataset['equivalent_atoms'])
    sym_independ = np.unique(dataset['equivalent_atoms'])
    print(len(sym_independ))
    print(sym_independ)
    for i in sym_independ:
        print(positions[i])
        #print((s2.sites[i])._fcoords)

def get_Symmetry(atmnt, vornet):
    positions = []
    #print(positions)
    lattice = atmnt.lattice
    for i in vornet.nodes:
        positions.append(atmnt.absolute_to_relative(i[1][0],i[1][1],i[1][2]))

    numbers = [1,]*len(vornet.nodes)
    
    # print(numbers)
    # print(len(positions))
    # print(len(numbers))
    cell = (lattice, positions, numbers)

    spacegroup = spglib.get_spacegroup(cell, symprec=0.01, angle_tolerance=5)
    print(spacegroup)  
    #symmetry = spglib.get_symmetry(cell, symprec=1e-5)
    #print(symmetry)
    dataset = spglib.get_symmetry_dataset(cell, symprec=0.01, angle_tolerance=5)
    print(len(dataset['equivalent_atoms']))
    #print(dataset['rotations'])
    #print(dataset['translations'])
    print(dataset['equivalent_atoms'])
    sym_independ = np.unique(dataset['equivalent_atoms'])

    for i in range(len(vornet.nodes)):
        #此处需修改定义
        #vornet.nodes[i] = [sym_independ[i],vornet.nodes[i][1], vornet.nodes[i][2]]

    # print(len(sym_independ))
    # print(sym_independ)
    for i in sym_independ:
        print(positions[i])
        #print((s2.sites[i])._fcoords)
    
    return vornet

# if __name__ == "__main__":
#     #get_Symmetry("../../examples/icsd_246817_orgin_copy.vasp")
#     # get_Symmetry("../../examples/Li2CO3-LDA_orgin.vasp")
#     #get_Symmetry("../../examples/LPS.vasp")
#     #get_Symmetry("../../examples/LPS_orgin.vasp")
#     #get_Symmetry("../../examples/icsd_20610.vasp")
#     #get_Symmetry("../../examples/icsd_20610_orgin.vasp")
#     #get_Symmetry("../../examples/icsd_29225.vasp")
#     # get_Symmetry("../../examples/icsd_29225_orgin.vasp")
#     #get_Symmetry("../../examples/LLZO-48g-180721_orgin.vasp")
#     #get_Symmetry("../../examples/custom_300001.vasp")
#     #get_Symmetry("../../examples/custom_300001_orgin.vasp")

#     radii = {}
#     remove_filename = getRemoveMigrantFilename("../../examples/Li2CO3-LDA.cif","Li")
#     radii,migrant_radius,migrant_alpha = LocalEnvirCom("../../examples/Li2CO3-LDA.cif","Li")
#     atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, True, None)
#     vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
#     writeVaspFile("../../examples/Li2CO3-LDA"+"_orgin.vasp",atmnet,vornet,True)
#     sym_vornet = get_Symmetry(atmnet, vornet)