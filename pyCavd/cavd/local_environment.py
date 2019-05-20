'''
离子局部环境分析程序。
更新日期：20190123
作者：YAJ
School of Computer Engineering and Science ShangHai University 

配位数计算原理：
    1.按距离与正负号计算配位：
        以目标原子为中心，10A为半径画球形区域，得到区域内所有原子。按照距离目标离子中心的距离排序从小到大排序，存入列表中。
        从前往后遍历列表，判断每个原子的正负，直到出现第一个同号原子为止，记录下的异号原子的数目即为配位数。（计算第一配位球壳层）
    2.按O'Keefee 提供的方法计算配位。
        原始文献为：
            [1]	M. O’Keeffe, “A proposed rigorous definition of coordination number,” Acta Crystallogr. Sect. A, vol. 35, no. 5, pp. 772–775, 1979.
        下列文献对该方法做了介绍：
            [2]	L. Yang, S. Dacek, and G. Ceder, “Proposed definition of crystal substructure and substructural similarity,” Phys. Rev. B - Condens. Matter Mater. Phys., vol. 90, no. 5, pp. 1–9, 2014.

离子半径计算原理：
    计算给定结构文件中对应位点的配位数后，依据元素、价态和配位数三者的信息查香农1976年有效离子半径表，得到对应的半径。
    ！！！----------注意----------！！！
    若香农表中无法找到结构文件中离子对应的价态或配位数，则会返回香农表中最接近文件中的价态或配位数所对应的半径，作为离子半径。
    若结构文件中为给定价态，使用原子半径作为替代。

输入：
    结构文件、迁移离子元素符号
    离子半径表：香农1976年离子半径。(ionic_radii.json)
输出：
    格式：
    Coordination组成的列表、(label,radius)组成的列表。
    Coordination定义为:
        Coordination：
            label element coord frac_coord radius coord_num(配位数) coord_nei
    分别指：
        label：cif文件中的的“_atom_site_label”字段对应的值，表征不同环境下的原子
        element：中心原子元素符号
        coord：中心原子的直角坐标
        frac_coord：中心原子的直角坐标
        radius：中心原子的离子半径（在结构文件中未给出价态时，使用原子半径替代）
        coord_num：中心原子的配位数
        coord_nei：中心原子的配位原子列表。列表中的每一项为(PeriodicSite,distance)，表示一个配位原子：
            PeriodicSite为Pymatgen中表示位点的数据结构；
            distance为中心原子到配位原子的距离，单位为A。
    label: radius 组成的列表。
环境要求：
    需安装pymatgen包
    需将ionic_radii.json文件与该文件（ionic_radii.py）放置在同一目录下。
使用方法：
    更改“if __name__ == "__main__":”模块中调用的“get_ionic_radii()”中的参数为需要计算的cif文件，并加上迁移离子。
程序目前的缺陷：
    无法计算具有部分占据/混占的结构文件。
    对于香农1976年有效离子半径表中不存在的价态（小数价态）、元素（H与H的同位素等）和配位数组合，会采用近似值的方法查表。
    仅保证针对对离子晶体计算有效。
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

class CoordEnviroComError(Exception):
    def __init__(self, message = "Coordination Environment Compute Error!"):
        self.message = message

class RadiusComError(Exception):
    def __init__(self, message = "Cannot assign radius!"):
        self.message = message

file_dir = os.path.dirname(__file__)
rad_file = os.path.join(file_dir, 'ionic_radii.json')
with open(rad_file, 'r') as fp:
    _ion_radii = json.load(fp)

# 自定义的Cif文件解析类,实现获取cif结构中的label
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
    
    
    def get_symme(self):
        for d in self._cif.data.values():
            data = d
            # Try to parse International number
            for symmetry_label in ["_space_group_IT_number",
                                    "_space_group_IT_number_",
                                    "_symmetry_Int_Tables_number",
                                    "_symmetry_Int_Tables_number_"]:
                if data.data.get(symmetry_label):
                    symm_number = int(str2float(data.data.get(symmetry_label)))
                    break
            for symmetry_label in ["_symmetry_space_group_name_H-M",
                                   "_symmetry_space_group_name_H_M",
                                   "_symmetry_space_group_name_H-M_",
                                   "_symmetry_space_group_name_H_M_",
                                   "_space_group_name_Hall",
                                   "_space_group_name_Hall_",
                                   "_space_group_name_H-M_alt",
                                   "_space_group_name_H-M_alt_",
                                   "_symmetry_space_group_name_hall",
                                   "_symmetry_space_group_name_hall_",
                                   "_symmetry_space_group_name_h-m",
                                   "_symmetry_space_group_name_h-m_"]:
                if data.data.get(symmetry_label):
                    sg = data.data.get(symmetry_label)
                    break
        return symm_number,sg

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

# pymtgen实现的求离子半径
def get_ionic_radii_pymatgen(filename):
    stru = Structure.from_file(filename)
    val_eval = ValenceIonicRadiusEvaluator(stru)
    return val_eval.radii

class VoronoiNN_self(VoronoiNN):
    def __init__(self, tol=0.5, targets=None, cutoff=10.0,
                 allow_pathological=False, weight='solid_angle',
                 extra_nn_info=True):
        super(VoronoiNN_self, self).__init__()
        self.tol = tol
        self.cutoff = cutoff
        self.allow_pathological = allow_pathological
        self.targets = targets
        self.weight = weight
        self.extra_nn_info = extra_nn_info

    #使用solid angle计算配位数
    def get_cn_solidangle(self, structure, n, use_weights=False):
        siw = self.get_nn_info(structure, n)
        return sum([e['weight'] for e in siw]) if use_weights else len(siw)

    #根据距离与异号原则计算配位数
    def get_cn_dis(self, structure, n):
        cn = 0
        center = structure[n]
        neighbors = structure.get_neighbors(center, self.cutoff)
        neighbors = [i for i in sorted(neighbors, key=lambda s: s[1])]

        if '+' in center.species_string:
            for i in range(len(neighbors)):
                if "-" in neighbors[i][0].species_string:
                   cn = cn + 1
                if "+" in neighbors[i][0].species_string:
                    break
        if "-" in center.species_string:
            for i in range(len(neighbors)):
                if "+" in neighbors[i][0].species_string:
                    cn = cn + 1
                if "-" in neighbors[i][0].species_string:
                    break
        return cn
    
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
    def get_coor_list(sites, radii, coord_neighbors):
        coordination_list = []
        if len(sites) == len(radii) == len(coord_neighbors):
            if len(coord_neighbors) == 0:
                raise CoordEnviroComError("Coordination Number is 0!")
            else:
                for i in range(len(sites)):
                    site = sites[i]
                    label = site._atom_site_label
                    coordination_list.append(Coordination(label, site.coords, site.frac_coords, site.specie.symbol, radii[label], coord_neighbors[label]).as_dict())
                return coordination_list
        else:
            raise CoordEnviroComError("labels, sites, elements, coord_neighbors, radii must be have same length!")
    
    #根据label从coordination列表中查找Coordination
    def get_coordination(coordination_list,label):
        for i in coordination_list:
            if i["label"] == label:
                return i

    def as_dict(self):
        """
        Dict representation of Coordination.
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

# 计算原子半径
# 注意：仅在结构文件中不含任何价态信息时调用
# If valence is zero, atomic radius is used.
def get_atomic_radius(site):
    radius = site.specie.atomic_radius
    # Handle elements with no atomic_radius
    # by using calculated values instead.
    if radius is None:
        radius = site.specie.atomic_radius_calculated
    return radius

#根据元素、价态、配位数查询香农表获取离子半径
def ger_ionic_radius(elem,oxi_state,coord_no):
    try:
        tab_oxi_states = sorted(map(int, _ion_radii[elem].keys()))
        oxi_state = nearest_key(tab_oxi_states, oxi_state)
        tab_coord_noes = sorted(map(int, _ion_radii[elem][str(oxi_state)].keys()))
        radius = _ion_radii[elem][str(oxi_state)][str(coord_no)]
    except KeyError:
        coord_num = coord_no
        if coord_num - coord_no > 0:
            new_coord_no = coord_no + 1
        else:
            new_coord_no = coord_no - 1
        try:
            radius = _ion_radii[elem][str(oxi_state)][str(new_coord_no)]
            coord_no = new_coord_no
        except:
            tab_coords = sorted(map(int, _ion_radii[elem][str(oxi_state)].keys()))
            new_coord_no = nearest_key(tab_coords, coord_no)
            i = 0
            for val in tab_coords:
                if  val > coord_no:
                    break
                i = i + 1
            if i == len(tab_coords):
                key = str(tab_coords[-1])
                radius = _ion_radii[elem][str(oxi_state)][key]
            elif i == 0:
                key = str(tab_coords[0])
                radius = _ion_radii[elem][str(oxi_state)][key]
            else:
                key = str(tab_coords[i-1])
                radius1 = _ion_radii[elem][str(oxi_state)][key]
                key = str(tab_coords[i])
                radius2 = _ion_radii[elem][str(oxi_state)][key]
                radius = (radius1+radius2)/2
    return radius

#获取指定位点半径值
def get_radius_value(site,coord_no):
    el = site.specie.symbol
    oxi_state = int(round(site.specie.oxi_state))
    if isinstance(site.specie, Element):
        radius = get_atomic_radius(site)
    else:
        radius = ger_ionic_radius(el,oxi_state,coord_no)
    return radius

#从结构中获取半径
def get_radii_stru(stru):       
    radii = []
    labels = [] 
    vnn = VoronoiNN_self(tol=0.5, cutoff=10.0)
    for i in range(len(stru.sites)):
        site = stru.sites[i]
        label = site._atom_site_label
        #按label获取数据
        if label in labels:
            continue
        coord_no= vnn.get_cn_solidangle(stru, i)
        labels.append(label)
        radius = get_radius_value(site,coord_no)
        if  radius is None:
            raduis = 0
        radii.append(radius)
    label_radii_dict = dict(zip(labels, radii))
    return label_radii_dict

"""
Computes ionic radii of elements for all sites in the structure.
If valence is zero, atomic radius is used.
"""
def get_radii(filename):
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    return get_radii_stru(stru)

#get coordination number from structure
def get_cns_stru(stru):
    labels = []
    cns = []
    vnn = VoronoiNN_self(tol=0.5, cutoff=10.0)
    for i in range(len(stru.sites)):
        site = stru.sites[i]
        label = site._atom_site_label
        if label in labels:
            continue
        coord_no= vnn.get_cn_solidangle(stru, i)
        labels.append(label)
        cns.append(coord_no)
    label_radii_dict = dict(zip(labels, cns))
    return label_radii_dict

#get coordination number from file
def get_cns(filename):
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    return get_cns_stru(stru)  

#获取离中心原子半径cutoff范围内最近的coordnum个原子
def get_nearest_atoms(stru, ind, coordnum, cutoff = 10.0):
    centre = stru.sites[ind]
    neighbors = stru.get_neighbors(centre, cutoff)
    neighbors = [i for i in sorted(neighbors, key=lambda s: s[1])]
    neighbors = neighbors[0:coordnum]
    return neighbors

def get_local_envir_fromstru(stru):
    labels = []
    sites = []
    radii = {}
    cnatoms = {}
    vnn = VoronoiNN_self(tol=0.5, cutoff=10.0)
    for i in range(len(stru.sites)):
        site = stru.sites[i]
        label = site._atom_site_label
        if label not in labels:
            labels.append(label)
            coord_no= vnn.get_cn_solidangle(stru, i)
            neighbors = get_nearest_atoms(stru,i,coord_no,cutoff=10.0)
            sites.append(site)
            radii[label] = get_radius_value(site,coord_no)
            cnatoms[label] = neighbors
    return Coordination.get_coor_list(sites, radii, cnatoms),radii

def get_local_envir(filename):
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    symm_number,symm_sybol = parser.get_symme()
    coord_list,radii = get_local_envir_fromstru(stru)
    return symm_number,symm_sybol, coord_list, radii

# if __name__ == "__main__":
#     a,b= get_local_envir("./Li_Na_Mg_Al_cifs/Li/icsd_281589.cif")
#     print(a)
#     print(b)