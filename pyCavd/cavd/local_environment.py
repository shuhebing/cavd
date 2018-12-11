'''
离子局部环境分析程序。
更新日期：20180710
作者：YAJ
配位数计算原理：
    以目标离子为中心，10A为半径画球形区域，得到区域内的离子坐标。之后，按照距离目标离子中心的距离排序并存入一个列表中。
    从前往后遍历列表，判断每个离子的正负，直到出现第一个同号离子为止，此时记录下的异号离子数目即为配位数。（计算第一配位球壳层）
离子半径计算原理：
    计算给定结构文件中对应位点的配位数后，依据元素、价态和配位数三者的信息查香农1976年有效离子半径表，得到对应的半径。
alpha计算原理：
    计算出迁移离子到配位离子表面距离减去配位离子的半径，再除以迁移离子半径。max_alpha、min_alpha分别对应最远与最近配位离子。
输入：
    结构文件、迁移离子元素符号
    离子半径表：香农1976年离子半径。(ionic_radii.json)
输出：
    格式：
    自定义数据结构Coordination组成的列表。Coordination的数据结构为:label element coord frac_coord radius coord_num(配位数) coord_nei
        coord_nei为列表，列表中的每一项为(label,distance)，其中label为配位原子的label，distance为中心原子到配位原子的距离
    label: radius 组成的列表。
环境要求：
    需安装pymatgen包
    需将ionic_radii.json文件放置到同与该文件（ionic_radii.py）同一文件目录下。
使用方法：
    更改“if __name__ == "__main__":”模块中调用的“get_ionic_radii()”中的参数为需要计算的cif文件，并加上迁移离子。
程序目前的缺陷：
    无法计算具有部分占据/混占的结构文件。
    无法计算香农1976年有效离子半径表中不包含的离子半径（含H元素等）。
    仅适用离子晶体计算。
'''

import os
import json
import math
import six
import warnings
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator
from bisect import bisect_left
from pymatgen.io.cif import CifFile
from pymatgen.io.cif import CifParser
from pymatgen.io.cif import str2float
from monty.io import zopen
from six.moves import zip, cStringIO
import numpy as np
from itertools import groupby
from pymatgen.util.coord import in_coord_list_pbc, find_in_coord_list_pbc
from pymatgen.core.composition import Composition
from pymatgen.electronic_structure.core import Magmom
from pymatgen.core.operations import MagSymmOp
from collections import OrderedDict
from pymatgen.core.periodic_table import Element, Specie, get_el_sp, DummySpecie


file_dir = os.path.dirname(__file__)
#file_dir = os.path.abspath("/home/yeanjiang/yaj/CAVD/pyCavd/zeo/")
rad_file = os.path.join(file_dir, 'ionic_radii.json')
with open(rad_file, 'r') as fp:
    _ion_radii = json.load(fp)

# pymtgen实现的求离子半径
def get_ionic_radii_pymatgen(filename):
    stru = Structure.from_file(filename)
    val_eval = ValenceIonicRadiusEvaluator(stru)
    radii = val_eval.radii
    print("pymatgen radii:",radii)

class VoronoiNN_self(VoronoiNN):
    def __init__(self, tol=0, targets=None, cutoff=10.0,
                 allow_pathological=False, weight='solid_angle',
                 extra_nn_info=True):
        super(VoronoiNN_self, self).__init__()
        self.tol = tol
        self.cutoff = cutoff
        self.allow_pathological = allow_pathological
        self.targets = targets
        self.weight = weight
        self.extra_nn_info = extra_nn_info
    #根据目标离子设定范围球形区域内的异号离子数量确定配位数，考虑周期性情况
    def get_cn(self, structure, n):
        #print(structure)
        center = structure[n]
        neighbors = structure.get_sites_in_sphere(center.coords, self.cutoff)
        neighbors = [i[0] for i in sorted(neighbors, key=lambda s: s[1])]
        #print(center)
        #print(neighbors)
        # 选择异号离子
        if '+' in center.species_string:
            new_neighbors = []
            # 判断方式一：找到邻居内所有异号离子
#             for i in neighbors:
#                 if "-" in i.species_string:
#                     new_neighbors.append(i)
            # 判断方式二：找到排序后（从小到大排）的邻居内的离子，直到出现同号离子为止
            for i in range(len(neighbors)):
                if i is 0:
                    continue
                if "-" in neighbors[i].species_string:
                    new_neighbors.append(neighbors[i])
                if "+" in neighbors[i].species_string:
                    break
            neighbors = new_neighbors
               
        if "-" in center.species_string:
            new_neighbors = []
#             for i in neighbors:
#                 if "+" in i.species_string:
#                     new_neighbors.append(i)
            for i in range(len(neighbors)):
                if i is 0:
                    continue
                if "+" in neighbors[i].species_string:
                    new_neighbors.append(neighbors[i])
                if "-" in neighbors[i].species_string:
                    break
            neighbors = new_neighbors
        return len(neighbors),neighbors
    

# 自定义的Cif文件解析类
class CifParser_new(CifParser):
    """
    Parses a cif file

    Args:
        filename (str): Cif filename. bzipped or gzipped cifs are fine too.
        occupancy_tolerance (float): If total occupancy of a site is between 1
            and occupancy_tolerance, the occupancies will be scaled down to 1.
        site_tolerance (float): This tolerance is used to determine if two
            sites are sitting in the same position, in which case they will be
            combined to a single disordered site. Defaults to 1e-4.
    """

    def __init__(self, filename, occupancy_tolerance=1., site_tolerance=1e-4):
        self._occupancy_tolerance = occupancy_tolerance
        self._site_tolerance = site_tolerance
        if isinstance(filename, six.string_types):
            self._cif = CifFile.from_file(filename)
        else:
            self._cif = CifFile.from_string(filename.read())

        # store if CIF contains features from non-core CIF dictionaries
        # e.g. magCIF
        self.feature_flags = {}
        self.errors = []
        
        def is_magcif():
            """
            Checks to see if file appears to be a magCIF file (heuristic).
            """
            # Doesn't seem to be a canonical way to test if file is magCIF or
            # not, so instead check for magnetic symmetry datanames
            prefixes = ['_space_group_magn', '_atom_site_moment',
                        '_space_group_symop_magn']
            for d in self._cif.data.values():
                for k in d.data.keys():
                    for prefix in prefixes:
                        if prefix in k:
                            return True
            return False

        self.feature_flags['magcif'] = is_magcif()

        def is_magcif_incommensurate():
            """
            Checks to see if file contains an incommensurate magnetic
            structure (heuristic).
            """
            # Doesn't seem to be a canonical way to test if magCIF file
            # describes incommensurate strucure or not, so instead check
            # for common datanames
            if not self.feature_flags["magcif"]:
                return False
            prefixes = ['_cell_modulation_dimension', '_cell_wave_vector']
            for d in self._cif.data.values():
                for k in d.data.keys():
                    for prefix in prefixes:
                        if prefix in k:
                            return True
            return False

        self.feature_flags['magcif_incommensurate'] = is_magcif_incommensurate()

        for k in self._cif.data.keys():
            # pass individual CifBlocks to _sanitize_data
            self._cif.data[k] = self._sanitize_data(self._cif.data[k])

    @staticmethod
    def from_string(cif_string, occupancy_tolerance=1.):
        """
        Creates a CifParser from a string.

        Args:
            cif_string (str): String representation of a CIF.
            occupancy_tolerance (float): If total occupancy of a site is
                between 1 and occupancy_tolerance, the occupancies will be
                scaled down to 1.

        Returns:
            CifParser
        """
        stream = cStringIO(cif_string)
        return CifParser_new(stream, occupancy_tolerance)

    def _unique_coords(self, coords_in, magmoms_in=None, lattice=None):
        """
        Generate unique coordinates using coord and symmetry positions
        and also their corresponding magnetic moments, if supplied.
        """
        coords = []
        coords_num = []
        
        if magmoms_in:
            magmoms = []
            if len(magmoms_in) != len(coords_in):
                raise ValueError
            for tmp_coord, tmp_magmom in zip(coords_in, magmoms_in):
                count = 0
                for op in self.symmetry_operations:
                    coord = op.operate(tmp_coord)
                    coord = np.array([i - math.floor(i) for i in coord])
                    if isinstance(op, MagSymmOp):
                        # Up to this point, magmoms have been defined relative
                        # to crystal axis. Now convert to Cartesian and into
                        # a Magmom object.
                        magmom = Magmom.from_moment_relative_to_crystal_axes(
                            op.operate_magmom(tmp_magmom),
                            lattice=lattice
                        )
                    else:
                        magmom = Magmom(tmp_magmom)
                    if not in_coord_list_pbc(coords, coord,
                                             atol=self._site_tolerance):
                        coords.append(coord)
                        magmoms.append(magmom)
                        count = count + 1
                coords_num.append(count)   
            return coords, magmoms, coords_num
        else:
            for tmp_coord in coords_in:
                count = 0
                for op in self.symmetry_operations:
                    coord = op.operate(tmp_coord)
                    coord = np.array([i - math.floor(i) for i in coord])
                    if not in_coord_list_pbc(coords, coord,
                                             atol=self._site_tolerance):
                        coords.append(coord)
                        count = count + 1
                coords_num.append(count)
            return coords, [Magmom(0)] * len(coords), coords_num  # return dummy magmoms

    def _get_structure(self, data, primitive):
        """
        Generate structure from part of the cif.
        """

        def get_num_implicit_hydrogens(sym):
            num_h = {"Wat": 2, "wat": 2, "O-H": 1}
            return num_h.get(sym[:3], 0)

        lattice = self.get_lattice(data)

        # if magCIF, get magnetic symmetry moments and magmoms
        # else standard CIF, and use empty magmom dict
        if self.feature_flags["magcif_incommensurate"]:
            raise NotImplementedError(
                "Incommensurate structures not currently supported.")
        elif self.feature_flags["magcif"]:
            self.symmetry_operations = self.get_magsymops(data)
            magmoms = self.parse_magmoms(data, lattice=lattice)
        else:
            self.symmetry_operations = self.get_symops(data)
            magmoms = {}

        oxi_states = self.parse_oxi_states(data)
        coord_to_species = OrderedDict()
        coord_to_magmoms = OrderedDict()

        def get_matching_coord(coord):
            keys = list(coord_to_species.keys())
            coords = np.array(keys)
            for op in self.symmetry_operations:
                c = op.operate(coord)
                inds = find_in_coord_list_pbc(coords, c,
                                              atol=self._site_tolerance)
                # cant use if inds, because python is dumb and np.array([0]) evaluates
                # to False
                if len(inds):
                    return keys[inds[0]]
            return False
        
        label_el_dict = {}
        
        for i in range(len(data["_atom_site_label"])):
            try:
                # If site type symbol exists, use it. Otherwise, we use the
                # label.
                symbol = self._parse_symbol(data["_atom_site_type_symbol"][i])
                
                label = data["_atom_site_label"][i]
                
                num_h = get_num_implicit_hydrogens(
                    data["_atom_site_type_symbol"][i])
                
            except KeyError:
                symbol = self._parse_symbol(data["_atom_site_label"][i])
                
                label = data["_atom_site_label"][i]
                
                num_h = get_num_implicit_hydrogens(data["_atom_site_label"][i])
            if not symbol:
                continue

            if oxi_states is not None:
                o_s = oxi_states.get(symbol, 0)
                # use _atom_site_type_symbol if possible for oxidation state
                if "_atom_site_type_symbol" in data.data.keys():
                    oxi_symbol = data["_atom_site_type_symbol"][i]
                    o_s = oxi_states.get(oxi_symbol, o_s)
                try:
                    el = Specie(symbol, o_s)
                except:
                    el = DummySpecie(symbol, o_s)
            else:
                el = get_el_sp(symbol)
            
            x = str2float(data["_atom_site_fract_x"][i])
            y = str2float(data["_atom_site_fract_y"][i])
            z = str2float(data["_atom_site_fract_z"][i])
            magmom = magmoms.get(data["_atom_site_label"][i],
                                 np.array([0, 0, 0]))
            
            try:
                occu = str2float(data["_atom_site_occupancy"][i])
            except (KeyError, ValueError):
                occu = 1

            if occu > 0:
                coord = (x, y, z)
                match = get_matching_coord(coord)
                comp_d = {el: occu}
                
                if num_h > 0:
                    comp_d["H"] = num_h
                comp = Composition(comp_d)
                if not match:
                    coord_to_species[coord] = comp
                    coord_to_magmoms[coord] = magmom
                else:
                    coord_to_species[match] += comp
                    # disordered magnetic not currently supported
                    coord_to_magmoms[match] = None
            
            label_el_dict[coord] = label
                
        sum_occu = [sum(c.values()) for c in coord_to_species.values()
                    if not set(c.elements) == {Element("O"), Element("H")}]
        if any([o > 1 for o in sum_occu]):
            msg = "Some occupancies (%s) sum to > 1! If they are within " \
                    "the tolerance, they will be rescaled." % str(sum_occu)
            warnings.warn(msg)
            self.errors.append(msg)
        
        allspecies = []
        allcoords = []
        allmagmoms = []
        allhydrogens = []
        alllabels = []

        # check to see if magCIF file is disordered
        if self.feature_flags["magcif"]:
            for k, v in coord_to_magmoms.items():
                if v is None:
                    # Proposed solution to this is to instead store magnetic
                    # moments as Specie 'spin' property, instead of site
                    # property, but this introduces ambiguities for end user
                    # (such as unintended use of `spin` and Specie will have
                    # fictious oxidation state).
                    raise NotImplementedError(
                        'Disordered magnetic structures not currently supported.')
        if coord_to_species.items():
            for comp, group in groupby(
                    sorted(list(coord_to_species.items()), key=lambda x: x[1]),
                    key=lambda x: x[1]):
                tmp_coords = [site[0] for site in group]
                
                #print(tmp_coords)
                labels = []
                for i in tmp_coords:
                    labels.append(label_el_dict[i])
                    
                #print(labels)
                tmp_magmom = [coord_to_magmoms[tmp_coord] for tmp_coord in
                              tmp_coords]
                if self.feature_flags["magcif"]:
                    coords, magmoms, coords_num = self._unique_coords(tmp_coords,
                                                          magmoms_in=tmp_magmom,
                                                          lattice=lattice)
                else:
                    coords, magmoms, coords_num = self._unique_coords(tmp_coords)

                if set(comp.elements) == {Element("O"), Element("H")}:
                    # O with implicit hydrogens
                    im_h = comp["H"]
                    species = Composition({"O": comp["O"]})
                else:
                    im_h = 0
                    species = comp
                    
                allhydrogens.extend(len(coords) * [im_h])
                allcoords.extend(coords)
                allspecies.extend(len(coords) * [species])
                allmagmoms.extend(magmoms)
                
                for i in range(len(coords_num)):

                    alllabels.extend(coords_num[i] * [labels[i]])
                
            # rescale occupancies if necessary
            for i, species in enumerate(allspecies):
                totaloccu = sum(species.values())
                if 1 < totaloccu <= self._occupancy_tolerance:
                    allspecies[i] = species / totaloccu

        if allspecies and len(allspecies) == len(allcoords) \
                and len(allspecies) == len(allmagmoms):
            site_properties = dict()
            if any(allhydrogens):
                assert len(allhydrogens) == len(allcoords)
                site_properties["implicit_hydrogens"] = allhydrogens

            if self.feature_flags["magcif"]:
                site_properties["magmom"] = allmagmoms

            if len(site_properties) == 0:
                site_properties = None
                
            struct = Structure(lattice, allspecies, allcoords,
                               site_properties=site_properties)
            #struct = struct.get_sorted_structure()

            if primitive and self.feature_flags['magcif']:
                struct = struct.get_primitive_structure(use_site_props=True)
            elif primitive:
                struct = struct.get_primitive_structure()
                struct = struct.get_reduced_structure()  
            struct.add_site_property("_atom_site_label", alllabels)
            return struct

#自定义数据结构 
class Coordination():
    def __init__(self, label, coord, frac_coord, element, radius=None, coord_neighbors=None):
        self._label = label
        self._element = element
        self._coord = coord
        self._frac_coord = frac_coord
        self._radius = radius
        self._coord_num = len(coord_neighbors)
        self._coord_nei = coord_neighbors
    #列表labels, elements, coord_nums, radii必须具有相同的大小
    @staticmethod
    def get_coor_list(labels, sites, elements, radii, distances, coord_neighbors):
        coordination_list = []
        if len(labels) != len(sites) != len(elements) != len(coord_neighbors) != len(radii) != len(distances):
            raise ValueError("labels, sites, elements, coord_neighbors, radii, distances must be have same length!")
        for i in range(len(labels)):
            coord_nei = []
            for j in range(len(coord_neighbors[i])):
                coord_nei.append((coord_neighbors[i][j]._atom_site_label, distances[i][j]))
            coordination_list.append(Coordination(labels[i], sites[i].coords, sites[i].frac_coords, elements[i], radii[i], coord_nei).as_dict())
        return coordination_list
    
    #根据label从coordination列表中查找Coordination
    def get_coordination(coordination_list,label):
        for i in coordination_list:
            if i["label"] == label:
                return i

    def as_dict(self):
        """
        Dict representation of Radius.
        Returns:
            JSON serializable dict representation.
        """
        d = {"label": self._label,
             "element": self._element,
             "coord": self._coord,
             "frac_coord": self._frac_coord,
             "radius": self._radius,
             "coord_num": self._coord_num,
             "coord_nei": self._coord_nei
        }
        
#         # 有序字典
#         d = OrderedDict()
#         d["label"] = self._label
#         d["element"] = self._element
#         d["coord_num"] = self._coord_num
#         d["radius"] =  self._radius
        return d
       
def get_local_envir(filename):
    """
    Computes ionic radii of elements for all sites in the structure.
    If valence is zero, atomic radius is used.
    """
    radii = []
    labels = []
    cnatoms = []
    els = []
    distances = []
    sites=[]
    
    vnn = VoronoiNN_self(cutoff = 10.0)
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    def nearest_key(sorted_vals, key):
        i = bisect_left(sorted_vals, key)
        if i == len(sorted_vals):
            return sorted_vals[-1]
        if i == 0:
            return sorted_vals[0]
        before = sorted_vals[i-1]
        after = sorted_vals[i]
        if after-key < key-before:
            return after
        else:
            return before
    
    for i in range(len(stru.sites)):
        site = stru.sites[i]
        
        if isinstance(site.specie, Element):
            radius = site.specie.atomic_radius
            # Handle elements with no atomic_radius
            # by using calculated values instead.
            if radius is None:
                radius = site.specie.atomic_radius_calculated
            if radius is None:
                raise ValueError(
                        "cannot assign radius to element {}".format(
                        site.specie))
            radii.append(radius)
            continue
            
        el = site.specie.symbol
        oxi_state = int(round(site.specie.oxi_state))
        coord_no, coord_neighbors = vnn.get_cn(stru, i)
        label = site._atom_site_label
        try:
            tab_oxi_states = sorted(map(int, _ion_radii[el].keys()))
            oxi_state = nearest_key(tab_oxi_states, oxi_state)
            radius = _ion_radii[el][str(oxi_state)][str(coord_no)]
        except KeyError:
            coord_num,neig = vnn.get_cn(stru, i)
            if coord_num - coord_no > 0:
                new_coord_no = coord_no + 1
            else:
                new_coord_no = coord_no - 1
            try:
                radius = _ion_radii[el][str(oxi_state)][str(new_coord_no)]
                coord_no = new_coord_no
            except:
                tab_coords = sorted(map(int, _ion_radii[el][str(oxi_state)].keys()))
                new_coord_no = nearest_key(tab_coords, coord_no)
                i = 0
                for val in tab_coords:
                    if  val > coord_no:
                        break
                    i = i + 1
                if i == len(tab_coords):
                    key = str(tab_coords[-1])
                    radius = _ion_radii[el][str(oxi_state)][key]
                elif i == 0:
                    key = str(tab_coords[0])
                    radius = _ion_radii[el][str(oxi_state)][key]
                else:
                    key = str(tab_coords[i-1])
                    radius1 = _ion_radii[el][str(oxi_state)][key]
                    key = str(tab_coords[i])
                    radius2 = _ion_radii[el][str(oxi_state)][key]
                    radius = (radius1+radius2)/2
        #implement complex checks later

        dis = []
        for i in range(len(coord_neighbors)):
            dis.append(site.distance(coord_neighbors[i]))

        if label in labels:
            continue

        # if migrant in site.species_string:
        #     print(site)
        #     print(coord_neighbors[0])
        #     mindis = site.distance(coord_neighbors[0])
        #     maxdis = site.distance(coord_neighbors[len(coord_neighbors)-1])
        #     migrant_para_tmp.append([label,radius,coord_neighbors[0],mindis,coord_neighbors[len(coord_neighbors)-1],maxdis])
       

        labels.append(label)
        sites.append(site)
        els.append(site.species_string)
        radii.append(radius)
        distances.append(dis)
        cnatoms.append(coord_neighbors)

        label_radii_list = dict(zip(labels, radii))

        # for i in migrant_para_tmp:
        #     minnei = i[2]._atom_site_label
        #     maxnei = i[4]._atom_site_label
        #     minrad = label_radii_list[minnei]
        #     maxrad = label_radii_list[maxnei]

        #     minalpa = (i[3]-minrad)/i[1]
        #     maxalpa = (i[5]-maxrad)/i[1]
        #     migrant_para.append([i[0],i[1],minnei,i[3],minrad,i[3],minalpa,maxnei,i[5],maxrad,i[5],maxalpa])
    
    return Coordination.get_coor_list(labels, sites, els, radii, distances, cnatoms),label_radii_list

# if __name__ == "__main__":
#     coord_list = get_local_envir("../../examples/icsd_16713.cif")
#     a = Coordination.get_coordination(coord_list,"Li1")
#     print(a)
