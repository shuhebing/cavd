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
from pymatgen.core.sites import PeriodicSite
import networkx as nx

import re
from cavd.netio import *
from monty.io import zopen
from cavd.channel import Channel
from cavd.netstorage import AtomNetwork, connection_values_list
from cavd.local_environment import CifParser_new, LocalEnvirCom
from scipy.spatial.ckdtree import cKDTree

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
The above code is developed for getting the pathways between experimental mobile ion sites.

Author: Anjiang Ye & Penghui Mi
School of Computer Engineering and Science, Shanghai University
2020 01 07

"""


class MigrationPaths(object):
    def __init__(self, struc, mobileIon, channels):
        self.struc = struc
        self.mobileIon = mobileIon
        if channels:
            self.channels = channels
        else:
            raise Exception("There are no channels in this structure!")
        
        self.interstices = []
        self.channelSegs = []
        self.network = None
       
        self.keyInterstices = []
        self.keyPaths = []
        self.keyInterMobileDict = {}
        
    
    def setInterstices(self):
        for channel in self.channels:
            self.interstices.extend(channel["nodes"])

    def setChannelSegs(self):
        for channel in self.channels:
            self.channelSegs.extend(channel["conns"])
    
    def setNetwork(self):
        graph = nx.Graph()
        for node in self.interstices:
            carts = [round(cart,5) for cart in node["cart_coord"]]
            fracs = [round(frac,5) for frac in node["frac_coord"]]
            graph.add_node(node["id"], id = node["id"], label=node["label"], cart_coord=carts, frac_coord=fracs)

        for edge in self.channelSegs:
            if edge["fromId"] < edge["toId"]:
                asc = edge["toDelta"]
                des = [-1*i for i in edge["toDelta"]]
                #利用hash值区分edge
                edgeStr = str(edge["fromLabel"]) + str(edge["toLabel"]) + str(round(edge["length"],2)) + str(round(edge["bottleneckSize"], 2))
                edgeHash = str(hash(edgeStr))
                
                graph.add_edge(edge["fromId"], edge["toId"], label=edgeHash, ascPDV=asc, desPDV=des)
        
        self.network = graph
    
    # 根据cif中迁移离子label设置关键间隙
    def setKeyInterstices(self):
        stru = self.struc
        migrant = self.mobileIon
        mobileCarts = np.around(np.array([site.coords for site in stru.sites if migrant in site._atom_site_label], ndmin=2), 5)
        mobileLabels = [site._atom_site_label for site in stru.sites if migrant in site._atom_site_label]
        
        interCarts = np.around(np.array([inter["cart_coord"] for inter in self.interstices], ndmin=2), 5)
        
        intersKdTree = cKDTree(interCarts)
        minDis,minIds = intersKdTree.query(mobileCarts,k=1)
        
        for idx in range(len(minIds)):
            tmpDict = {}
            curInter = self.interstices[minIds[idx]]
            tmpDict["interId"] = curInter["id"]
            tmpDict["targetIon"] = mobileLabels[idx]
            tmpDict["dis"] = minDis[idx]
            self.keyInterstices.append(tmpDict)
        
        self.setKeyInterMobileDict()
        
    def getKeyInterstices(self):
        return self.keyInterstices
    
    def setKeyInterMobileDict(self):
        for it in self.keyInterstices:
            self.keyInterMobileDict[it["interId"]] = it["targetIon"]
    
    # 获取sourceInter与sinkInter之间所有路径
    # legalPath 定义为除起、始点外，其余点均不为keyInterstices的Path
    def getPaths(self, suorceInter, sinkInter, cutoff = 5.0):
        legalPaths = []
        paths = list(nx.all_simple_paths(self.network, suorceInter, sinkInter, cutoff))
        keyInterIds = [it["interId"] for it in self.keyInterstices]
        for path in paths:
            legalTag = True
            for idx in range(1, len(path)-1):    
                if(self.interstices[path[idx]] in keyInterIds):
                    legalTag = False
                    break
            if legalTag:
                legalPaths.append(path)
        return legalPaths
    
    def setKeyPaths(self):
        originPaths = []
        pathTags = []
        
        # 获取所有通道
        for source in self.keyInterstices:
            for sink in self.keyInterstices:
                if(source["interId"] < sink["interId"]):
                    originPaths.extend(self.getPaths(source["interId"], sink["interId"]))
                    
        for path in originPaths:
            # 去除重复通道
            uniqTag = ""
            for idx in range(len(path)-1):
                uniqTag += self.network[path[idx]][path[idx+1]]["label"]    
            if uniqTag in pathTags:
                continue
            else:
                pathTags.append(uniqTag)
                tmpDict = {}
                # 通道起点/终点间隙id
                pathStart = self.network.node[path[0]]["id"]
                pathEnd = self.network.node[path[-1]]["id"]
                # 通道起点/终点对应的迁移离子label
                pathStartMobile = self.keyInterMobileDict[pathStart]
                pathEndMobile = self.keyInterMobileDict[pathEnd]
                
                itIdsOfPath = []
                itLabelsOfPath = []
                itCartOfPath = []
                itFracOfPath = []
                itPDVofPath = []
                for idx in range(len(path)):
                    itIdsOfPath.append(self.network.node[path[idx]]["id"])
                    itLabelsOfPath.append(self.network.node[path[idx]]["label"])
                    itCartOfPath.append(self.network.node[path[idx]]["cart_coord"])
                    itFracOfPath.append(self.network.node[path[idx]]["frac_coord"])
                    if idx == 0:
                       itPDVofPath.append([0,0,0])
                    else:
                        if path[idx-1] < path[idx]:
                            nodePDV = np.array(itPDVofPath[idx-1]) + np.array(self.network[path[idx-1]][path[idx]]["ascPDV"])
                            itPDVofPath.append(nodePDV.tolist())
                        else:
                            nodePDV = np.array(itPDVofPath[idx-1]) + np.array(self.network[path[idx-1]][path[idx]]["desPDV"])
                            itPDVofPath.append(nodePDV.tolist())
                            
                tmpDict["type"] = (pathStartMobile, pathEndMobile)
                tmpDict["ids"] = itIdsOfPath
                tmpDict["labels"] = itLabelsOfPath
                tmpDict["cart_coords"] = itCartOfPath
                tmpDict["frac_coords"] = itFracOfPath
                tmpDict["pdvs"] = itPDVofPath
                
                self.keyPaths.append(tmpDict)
    
    # 外部接口
    def getKeyPaths(self):
        return self.keyPaths
     
    def comKeyPaths(self):
        self.setInterstices()
        self.setChannelSegs()
        self.setKeyInterstices()
        
        keyInterstices = self.getKeyInterstices()
        for it in keyInterstices:
            print("It",it["interId"],"<-->",it["targetIon"],"dis:",it["dis"])
            
        self.setNetwork()
        self.setKeyPaths()
        return self.getKeyPaths()
        

def outVesta(filename, migrant, ntol=0.02, rad_flag=True, lower=0.0, upper=10.0, rad_dict=None):
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    
    species = [str(sp).replace("Specie ","") for sp in stru.species]
    elements = [re.sub('[^a-zA-Z]','',sp) for sp in species]
    sitesym = parser.get_sym_opt()
    if migrant not in elements:
        raise ValueError("The input migrant ion not in the input structure! Please check it.")
    effec_radii,migrant_radius,migrant_alpha,nei_dises,coordination_list = LocalEnvirCom(stru,migrant)
    
    radii = {}
    if rad_flag:
        if rad_dict:
            radii = rad_dict
        else:
            radii = effec_radii
    
    atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, rad_flag)
	
    vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(True, ntol)
	
    prefixname = filename.replace(".cif","")
    # compute the R_T
    conn_val = connection_values_list(prefixname+".resex", vornet)
    channels = Channel.findChannels2(vornet, atmnet, lower, upper, prefixname+".net")
    
    # output vesta file for visiualization
    Channel.writeToVESTA(channels, atmnet, prefixname)
    
    migratPath = MigrationPaths(stru, migrant, channels)
    keyPaths = migratPath.comKeyPaths()
    
    return conn_val, keyPaths


if __name__ == "__main__":
    conn_val, keyPaths = outVesta("icsd_16713.cif","Li",ntol=0.02, rad_flag=True, lower=0.0, upper=10.0)
    print(keyPaths)












# # The code need to be updated.

# #product the neb packages for one path
# #posc1:file of POSCAR1POSCAR1 posc2:file of POSCAR3 posc_path:file of POSCAR_path
# def path_poscar(posc1, posc2, posc_path):
    # struc1 = Structure.from_file(posc1)
    # struc2 = Structure.from_file(posc2)
    # s = Structure.from_file(posc_path)
    
    # path=[]
    # path.append(struc1.sites[0].frac_coords)
    # for site in s.sites:
        # if site.specie.symbol == 'He':
            # path.append(site.frac_coords)
    # path.append(struc2.sites[0].frac_coords)
    
    # nimages = len(path)-1
    # images=struc1.interpolate(struc2, nimages, True)
    # dir = os.path.dirname(posc_path)

    # i=0
    # for struc in images:
        # struc.translate_sites(0,path[i]-struc.sites[0].frac_coords,frac_coords=True, to_unit_cell=True)
        # num=('%02d' % i)
        # if not os.path.exists(dir+'/'+num):
            # os.mkdir(dir+'/'+num)
        # struc.to(filename=dir+'/'+num+'/POSCAR')
        # i=i+1

# def poscars_after_relax(dir):
    # if os.path.isdir(dir):
        # for f in os.listdir(dir):
            # f_dir = os.path.join(dir,f)
            # if os.path.isdir(f_dir) and 'relax_POSCAR1' in os.listdir(f_dir):
                # print(f_dir)
                # posc_path = os.path.join(f_dir, 'POSCAR_path')
                # # pos1 = os.path.join(os.path.join(f_dir, 'relax_POSCAR1'), 'CONTCAR')
                # # pos1 = os.path.join(os.path.join(f_dir, 'relax_POSCAR2'), 'CONTCAR')
                # pos1 = os.path.join(os.path.join(f_dir, 'relax_POSCAR1'), 'POSCAR')
                # pos2 = os.path.join(os.path.join(f_dir, 'relax_POSCAR2'), 'POSCAR')
                # path_poscar(pos1, pos2, posc_path)



