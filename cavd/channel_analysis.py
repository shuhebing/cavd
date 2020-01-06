'''
Used to analyze channel.

Created on 2019.5.8

@author: YeAnjiang
'''
import os
import spglib
import numpy as np
from pymatgen import Structure
from cavd.graphstorage import DijkstraNetwork

class Node:
    label = 0
    conn_num = 0
    def __init__(self,labe, num):
        self.label = labe
        self.conn_num = num
    def as_dict(self):
        d = {"label": self.label,
             "conn_num": self.conn_num
            }
        return d

class Edge:
    from_label = 0
    to_label = 0
    length = 0
    bt_radius = 0
    
    def __init__(self,f_label,t_label,len,rad):
        self.from_label = f_label
        self.to_label = t_label
        self.length = len
        self.bt_radius = rad
    def as_dict(self):
        d = {"from_label": self.from_label,
            "to_label": self.to_label,
            "length": self.length,
            "bt_radius": self.bt_radius
        }
        return d


def get_Symmetry_vornet(atmnt, vornet, symprec=0.01, angle_tolerance=5):
    positions = []
    lattice = atmnt.lattice
    for i in vornet.nodes:
        positions.append(i["coord"])
    numbers = [1,]*len(vornet.nodes)
   
    cell = (lattice, positions, numbers)

    #spacegroup = spglib.get_spacegroup(cell, symprec, angle_tolerance=5)
    #print(spacegroup)

    dataset = spglib.get_symmetry_dataset(cell, symprec, angle_tolerance=5)
    print(dataset['international'],dataset['number'])
    symm_label = dataset['equivalent_atoms']

    print(len(symm_label))
    print(symm_label)
    
    voidnet = DijkstraNetwork.from_VoronoiNetwork(vornet)
    print(vornet.nodes)
    print(voidnet.nodes)
    
    vornet_uni_symm = vornet.parse_symmetry(symm_label)
    voidnet_uni_sym = DijkstraNetwork.from_VoronoiNetwork(vornet_uni_symm)
    print(vornet_uni_symm.nodes)
    print(voidnet_uni_sym.nodes)
    
    sym_independ = np.unique(dataset['equivalent_atoms'])
    print(len(sym_independ))
    print(sym_independ)
    voids = []
    for i in sym_independ:
        voids.append(positions[i])
        #print(positions[i])
        #print((s2.sites[i])._fcoords)
    
    return vornet_uni_symm,voids
    
    
def get_nodes(filename):
    nodes = {}
    conns = {}
    flag_p = 0
    flag_n = 0
    file = open(filename, 'r')
    
    for line in file.readlines():
        if 'Interstitial' in line:
            flag_p = 1
            flag_n = 0
            continue
        if 'Connection' in line:
            flag_p = 0
            flag_n = 1
            continue
        if flag_p == 1:
            line = line.split('\t')
            if len(line) == 4:
                node_id = int(line[0])
                node = {
                    "void_id": int(line[0]),
                    "label": int(line[1]),
                    "fracs": list(map(float, line[2].split())),
                    "void_radius": float(line[3]),
                    "conn": []
                }
                nodes[node_id] = node
                conns[node_id] = []
                #nodes.append(node)
        if flag_n == 1:
            line = line.split('\t')
            if len(line) == 6:
                from_id = int(line[0])
                segment = {
                    "to_id": int(line[1]),
                    "to_label": -1,
                    "delta_pos": list(map(int, line[2].split())),
                    "bot_frac": list(map(float, line[3].split())),
                    "bot_radius": float(line[4]),
                    "conn_legth": float(line[5]),
                }
                conns[from_id].append(segment)
    return nodes,conns

def get_voidnet(nodes,conns):
    for key,value in nodes.items():
        void = nodes[key]
        void["conn"] = conns[key]
        for conn in void["conn"]:
            to_void = nodes[conn["to_id"]]
            #给conn["to_label"]赋值
            conn["to_label"] = to_void["label"]
    return nodes
    

def get_distict_channels(voidnet,preci=4,duplicate=False):
    distinct_voids = []
    visited_nodes = []

    for key,value in voidnet.items():
        void = voidnet[key]
        node = Node(void["label"],len(void["conn"])).as_dict()
        if node in visited_nodes:
            continue
        visited_nodes.append(node)
        
        distinct_channels = []
        visited_channels = []
        for conn in void["conn"]:
            to_void = voidnet[conn["to_id"]]
            if not duplicate and void["label"] > conn["to_label"]:
                continue
            edge = Edge(void["label"],to_void["label"],round(conn["bot_radius"],preci),round(conn["conn_legth"],preci)).as_dict()
            if edge in visited_channels:
                continue
            visited_channels.append(edge)
            distinct_channels.append(conn)
        void["conn"] = distinct_channels
        distinct_voids.append(void)
    return distinct_voids
                
def print_net(voidnet,file):
    out = open(file,"w")
    print("distinct voids:")
    out.write("distinct voids:\n")
    for void in voidnet:
        print(void["label"],void["fracs"],void["void_radius"])
        out.write(str(void["label"]) + "\t"+ str(void["fracs"]) + "\t" + str(void["void_radius"])+"\n")
    print()
    print("distinct channels:")
    out.write("distinct channels:\n")
    for void in voidnet:
        for ck in void["conn"]:
            print(void["label"],ck["to_label"],ck["conn_legth"],ck["bot_frac"],ck["bot_radius"])
            out.write(str(void["label"])+"\t"+str(ck["to_label"])+"\t"+str(ck["conn_legth"])+"\t"+str(ck["bot_frac"])+"\t"+str(ck["bot_radius"])+"\n")

def voids_classify(nodes,preci=4):
    rad_id_pairs = {}
    for key,value in nodes.items():
        radkey = round(nodes[key]["void_radius"],preci)
        if radkey not in rad_id_pairs.keys():
            rad_id_pairs[radkey] = []
        rad_id_pairs[radkey].append(nodes[key]["void_id"])
    return rad_id_pairs

def voids_relabel(nodes,rad_id_pairs):
    for value in rad_id_pairs.values():
        for id in value:
            nodes[id]["label"] = min(value)
            
    for key in nodes.keys():
        void = nodes[key]
        for conn in void["conn"]:
            to_void = nodes[conn["to_id"]]
            conn["to_label"] = to_void["label"]
    return nodes

def get_distinc(filename,preci=4,duplicate=False):
    nodes,conns = get_nodes(filename)
    voidnet = get_voidnet(nodes,conns)
    dist_voidnet = get_distict_channels(voidnet,preci,duplicate)
    
    if len(dist_voidnet) < 30:
        return dist_voidnet
    else:
        rad_label = voids_classify(voidnet)
        revoidnet = voids_relabel(voidnet,rad_label)
        dist_revoidnet = get_distict_channels(revoidnet,preci,duplicate)
        return dist_revoidnet

"""
# The code need to be updated.

#product the neb packages for one path
#posc1:file of POSCAR1POSCAR1 posc2:file of POSCAR3 posc_path:file of POSCAR_path
def path_poscar(posc1, posc2, posc_path):
    struc1 = Structure.from_file(posc1)
    struc2 = Structure.from_file(posc2)
    s = Structure.from_file(posc_path)
    
    path=[]
    path.append(struc1.sites[0].frac_coords)
    for site in s.sites:
        if site.specie.symbol == 'He':
            path.append(site.frac_coords)
    path.append(struc2.sites[0].frac_coords)
    
    nimages = len(path)-1
    images=struc1.interpolate(struc2, nimages, True)
    dir = os.path.dirname(posc_path)

    i=0
    for struc in images:
        struc.translate_sites(0,path[i]-struc.sites[0].frac_coords,frac_coords=True, to_unit_cell=True)
        num=('%02d' % i)
        if not os.path.exists(dir+'/'+num):
            os.mkdir(dir+'/'+num)
        struc.to(filename=dir+'/'+num+'/POSCAR')
        i=i+1

def poscars_after_relax(dir):
    if os.path.isdir(dir):
        for f in os.listdir(dir):
            f_dir = os.path.join(dir,f)
            if os.path.isdir(f_dir) and 'relax_POSCAR1' in os.listdir(f_dir):
                print(f_dir)
                posc_path = os.path.join(f_dir, 'POSCAR_path')
                # pos1 = os.path.join(os.path.join(f_dir, 'relax_POSCAR1'), 'CONTCAR')
                # pos1 = os.path.join(os.path.join(f_dir, 'relax_POSCAR2'), 'CONTCAR')
                pos1 = os.path.join(os.path.join(f_dir, 'relax_POSCAR1'), 'POSCAR')
                pos2 = os.path.join(os.path.join(f_dir, 'relax_POSCAR2'), 'POSCAR')
                path_poscar(pos1, pos2, posc_path)


"""
