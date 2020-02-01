'''
Used to analyze channel.

Created on 2019.5.8

@author: YeAnjiang
'''
import os
import spglib
import numpy as np
from cavd.graphstorage import DijkstraNetwork
from pymatgen.core.sites import PeriodicSite
import networkx as nx
import copy

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
        
        self.interstices = {}
        self.channelSegs = []
        self.network = None
       
        self.keyInterstices = []
        self.allPaths = []
        
        self.keyInterMobileDict = {}
        self.keyInterMobileIdxDict = {}
        
    
    def setInterstices(self):
        for channel in self.channels:
            for node in channel["nodes"]:
                tmpDic = {}
                tmpDic["label"] = node["label"]
                tmpDic["radius"] = node["radius"]
                tmpDic["cart_coord"] = [round(cart,5) for cart in node["cart_coord"]]
                tmpDic["frac_coord"] = [round(frac,5) for frac in node["frac_coord"]]
                self.interstices[node["id"]] = tmpDic

    def setChannelSegs(self):
        for channel in self.channels:
            self.channelSegs.extend(channel["conns"])
    
    def setNetwork(self):
        graph = nx.Graph()
        for key,value in self.interstices.items():
            carts = [round(cart,5) for cart in value["cart_coord"]]
            fracs = [round(frac,5) for frac in value["frac_coord"]]
            graph.add_node(key, label=value["label"], cart_coord=carts, frac_coord=fracs)

        for edge in self.channelSegs:
            if edge["fromId"] < edge["toId"]:
                asc = edge["toDelta"]
                des = [-1*i for i in edge["toDelta"]]
               
                # edgeStr = str(edge["fromLabel"]) + str(edge["toLabel"]) + str(round(edge["length"],2)) + str(round(edge["bottleneckSize"], 2))
                # edgeHash = str(hash(edgeStr))
                # graph.add_edge(edge["fromId"], edge["toId"], label=edgeHash, ascPDV=asc, desPDV=des)
                
                graph.add_edge(edge["fromId"], edge["toId"], length = edge["length"], ascPDV=asc, desPDV=des)
        
        self.network = graph
    
    
    def comDis(self, p1, p2):
        temp_site1 = PeriodicSite('Ar', p1, self.struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self.struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis
    
    # The key interstices are the interstices coincides the mobile ions.
    def setKeyInterstices(self):
        migrant = self.mobileIon
        
        sites = self.struc.sites
        mobiles = [(idx, sites[idx]) for idx in range(len(sites)) if migrant in sites[idx]._atom_site_label]
        iterList = list(self.interstices.items())
        
        # obtain the concidence interstice and mobile ions site
        for mobileSite in mobiles:
            minDis = 1000
            minId = -1
            for inter in iterList:
                iterId = inter[0]
                iterSite = inter[1]["frac_coord"]
                dis = self.comDis(mobileSite[1].frac_coords, iterSite)
                if dis < minDis:
                    minDis = dis
                    minId = iterId
            tmpDict = {}
            tmpDict["interId"] = minId
            tmpDict["mobileLabel"] = mobileSite[1]._atom_site_label
            tmpDict["mobileId"] = mobileSite[0]
            tmpDict["dis"] = minDis
            self.keyInterstices.append(tmpDict)
        
        self.setKeyInterMobileDict()
        self.setKeyInterMobileIdxDict()
        
    def getKeyInterstices(self):
        return self.keyInterstices
    
    def setKeyInterMobileDict(self):
        for it in self.keyInterstices:
            self.keyInterMobileDict[it["interId"]] = it["mobileLabel"]
    
    def setKeyInterMobileIdxDict(self):
        for it in self.keyInterstices:
            self.keyInterMobileIdxDict[it["interId"]] = it["mobileId"]
    

    # Obtain the paths between sourceInter and sinkInter.
    # if n>0 && n< Total(paths)ï¼return first shortest path to n-th shortest path. 
    # legalPath: midpoint of the path is not included in self.keyInterstices  
    # 
    def getPaths(self, suorceInter, sinkInter, n = 2, cutoff = 5.0):
        pathLenPairs = []
        paths = list(nx.all_simple_paths(self.network, suorceInter, sinkInter, cutoff))
        keyInterIds = [it["interId"] for it in self.keyInterstices]
        
        for path in paths:
            length = 0
            legalTag = True
            for idx in range(1,len(path)):
                length += self.network[path[idx-1]][path[idx]]["length"]
                
                if idx == len(path) - 1:
                    continue
                else:
                    midInter = self.interstices[path[idx]]
                    if(path[idx] in keyInterIds):
                        legalTag = False
                        break

            if legalTag:
                pathLenPairs.append((path, length))
        
        if pathLenPairs:
            pathLenPairs = sorted(pathLenPairs, key=lambda x: x[1],reverse=False)
            if n > 0 and n <= len(pathLenPairs):
                pathLenPairs = pathLenPairs[0:n]
            elif n > len(pathLenPairs):
                raise Exception("Given n is larger than the total number of paths between source and sink!!!")
            legalPaths = [path[0] for path in pathLenPairs]
            return legalPaths
        else:
            return []

    def setAllPaths(self, n=0):
        originPaths = []
        for source in self.keyInterstices:
            for sink in self.keyInterstices:
                if(source["interId"] < sink["interId"]):
                    originPaths.extend(self.getPaths(source["interId"], sink["interId"], n))
        
        symUniqPaths = []     
        for path in originPaths:
            itLabelsOfPath = []
            itCartOfPath = []
            itFracOfPath = []
            itPDVofPath = []
            for idx in range(len(path)):
                interId = path[idx]
                itLabelsOfPath.append(self.network.node[interId]["label"])
                itCartOfPath.append(self.network.node[interId]["cart_coord"])
                itFracOfPath.append(self.network.node[interId]["frac_coord"])
                if idx == 0:
                   itPDVofPath.append([0,0,0])
                else:
                    if path[idx-1] < path[idx]:
                        nodePDV = np.array(itPDVofPath[idx-1]) + np.array(self.network[path[idx-1]][path[idx]]["ascPDV"])
                        itPDVofPath.append(nodePDV.tolist())
                    else:
                        nodePDV = np.array(itPDVofPath[idx-1]) + np.array(self.network[path[idx-1]][path[idx]]["desPDV"])
                        itPDVofPath.append(nodePDV.tolist())
                        
            # only consider symmetry distinct channels
            if itLabelsOfPath in symUniqPaths:
                continue
            else:
                symUniqPaths.append(itLabelsOfPath)
                
                tmpDict = {}
                tmpDict["type"] = (self.keyInterMobileDict[itLabelsOfPath[0]], self.keyInterMobileDict[itLabelsOfPath[-1]])
                tmpDict["ids"] = path
                tmpDict["labels"] = itLabelsOfPath
                tmpDict["cart_coords"] = itCartOfPath
                tmpDict["frac_coords"] = itFracOfPath
                tmpDict["pdvs"] = itPDVofPath

                self.allPaths.append(tmpDict)
        
    def outPathToPOSCAR(self, fname):
        fileDir = fname + "_paths"
        if not os.path.exists(fileDir):
            os.mkdir(fileDir)
        for i in range(len(self.allPaths)):
            pathDir = fileDir + "/path_" + str(i)
            if not os.path.exists(pathDir):
                os.mkdir(pathDir)
           
            path = self.allPaths[i]
            startIdx = self.keyInterMobileIdxDict[path["ids"][0]]
            endIdx = self.keyInterMobileIdxDict[path["ids"][-1]]
    
            imageSpecies = self.struc.sites[startIdx].species
            # The image is labeled by "X-Y", where X is the label of source mobile ion site, and Y is the label of sink mobile ion site.
            imageLabel = self.keyInterMobileDict[path["ids"][0]] + "-" + self.keyInterMobileDict[path["ids"][-1]]
            imageProps = {"_atom_site_label": imageLabel}
            
            pathCoords = path["frac_coords"]
            for j in range(len(pathCoords)):
                num = ('%02d' %j)
                imageDir = pathDir + "/" +str(num)
                if not os.path.exists(imageDir):
                    os.mkdir(imageDir)
                
                stru = copy.deepcopy(self.struc)
                if j == 0:
                    imageCoords = stru.sites[startIdx].frac_coords
                elif j == len(pathCoords) - 1:
                    imageCoords = stru.sites[endIdx].frac_coords
                else: 
                    imageCoords = pathCoords[j]
                
                stru.remove_sites([startIdx])
                stru.insert(startIdx, imageSpecies, imageCoords, properties = imageProps)
                stru.remove_sites([endIdx])
                stru.to(filename=imageDir+"/POSCAR")

    # extern interface
    def comAllPaths(self, fname, n=2):
        self.setInterstices()
        self.setChannelSegs()
        self.setKeyInterstices()
        
        keyInterstices = self.getKeyInterstices()
        for it in keyInterstices:
            print("It"+str(it["interId"])+"-"+it["mobileLabel"],it["dis"])
            
        self.setNetwork()
        self.setAllPaths(n)
        self.outPathToPOSCAR(fname)
        return self.allPaths 



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

