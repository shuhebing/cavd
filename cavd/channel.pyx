# distutils: language = c++
# distutils: sources = ../channel.cc

"""
    Updated by Ye Anjiang to extend the Channel.
    yeanjiang@shu.edu.cn
    August 20, 2018
"""

from cavd.netstorage cimport VoronoiNetwork
from netstorage cimport AtomNetwork
from cavd.netstorage cimport VORONOI_NETWORK
from netstorage cimport ATOM_NETWORK
from channel cimport CHANNEL
from graphstorage cimport DELTA_POS
from graphstorage cimport CONN
from graphstorage cimport DIJKSTRA_NODE
from graphstorage cimport DIJKSTRA_NETWORK
from geometry cimport XYZ
from libcpp.vector cimport vector

#Added at 20180808
#Customize an exception class
class FindChannelError(Exception):
    #print("Find Channel in Voronoi Network Failed!")
    pass

cdef class Channel:
    """
    Python wrapper to Zeo++ Channel.
    """
    def __cinit__(self):
        self.thisptr = new CHANNEL()
    def __dealloc__(self):
        del self.thisptr

    property nodes:
        def __get__(self):
            nodes = []
            cdef vector[DIJKSTRA_NODE] dj_nodes = self.thisptr.nodes
            for i in range(dj_nodes.size()):
                dj_id = dj_nodes[i].id
                dj_label = dj_nodes[i].label
                dj_coords = [dj_nodes[i].x,dj_nodes[i].y,dj_nodes[i].z]
                dj_max_radius = dj_nodes[i].max_radius
                nodes.append([dj_id, dj_label, dj_coords, dj_max_radius])
            return nodes

    property connections:
        def __get__(self):
            connections = []
            cdef vector[CONN] conns
            cdef vector[DIJKSTRA_NODE] dj_nodes = self.thisptr.nodes
            for i in range(dj_nodes.size()):
                conns = dj_nodes[i].connections
                for j in range(conns.size()):
                    conn_from = conns[j].origin
                    conn_to = conns[i].ending
                    conn_length = conns[j].length
                    conn_bt = [conns[j].btx, conns[j].bty, conns[j].btz]
                    conn_max_radius = conns[j].max_radius
                    conn = [conn_from, conn_to, conn_length, conn_bt, conn_max_radius]
                    connections.append(conn)
            return connections
    
    property nodes_deltapos:
        def __get__(self):
            nodes_deltapos = []
            cdef vector[DELTA_POS] unitCells = self.thisptr.unitCells
            cdef vector[vector[int]] ucNodes = self.thisptr.ucNodes
            cdef DELTA_POS pos
            cdef vector[int] ucNode
            for i in range(unitCells.size()):
                pos = unitCells[i]
                ucNode = ucNodes[i]
                for j in range(ucNode.size()):
                    node_id = ucNode[j]
                    node_pos = [pos.x, pos.y, pos.z]
                    nodes_deltapos.append([node_id, node_pos])
            return nodes_deltapos
            
    property dimensionality:
        def __get__(self):
            return self.thisptr.dimensionality
    
    property lattice:
        def __get__(self):
            la = [self.thisptr.v_a.x, self.thisptr.v_a.y, self.thisptr.v_a.z]
            lb = [self.thisptr.v_b.x, self.thisptr.v_b.y, self.thisptr.v_b.z]
            lc = [self.thisptr.v_c.x, self.thisptr.v_c.y, self.thisptr.v_c.z]
            lattice = [la, lb, lc]
            return lattice
    @classmethod
    def findChannelsInVornet(cls, vornet, probe_rad, filename):
        cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')
        cdef char* c_filename = filename
        cdef vector[CHANNEL] channels
        if findChannels_new(c_vornet_ptr, probe_rad, &channels):
            if not c_writeToVMD(channels, c_filename):
                raise IOError
        else:
            raise FindChannelError
    
    @classmethod
    def findChannels(cls, vornet, atmnet, probe_rad, filename=None):
        cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
        cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
        cdef char* c_filename

        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef vector[CHANNEL] c_channels
        if not findChannels_new(c_vornet_ptr, probe_rad, &c_channels):
            raise FindChannelError
        if filename:
            c_filename = filename
            if not c_writeToNET(c_channels, c_filename, c_atmnet_ptr):
                raise IOError

        channels = []
        cdef vector[int] nodeIDs
        cdef DELTA_POS disp
        cdef DIJKSTRA_NODE curNode
        cdef DIJKSTRA_NODE otherNode
        cdef CONN curConn
        for i in range(c_channels.size()):
            nodes = {}
            conns = []
            channel = {}
            for j in range((c_channels[i].unitCells).size()):
                nodeIDs = (c_channels[i].ucNodes).at(j)
                disp = (c_channels[i].unitCells).at(j)
                for k in range(nodeIDs.size()):
                    curNode = (c_channels[i].nodes).at(nodeIDs.at(k))

                    for l in range((curNode.connections).size()):
                        curConn = curNode.connections.at(l)
                        otherNode = (c_channels[i].nodes).at(curConn.ending)
                        frac_coord = atmnet.absolute_to_relative(curConn.btx, curConn.bty, curConn.btz)
                        conns.append([curNode.id, [0, 0, 0], otherNode.id, [curConn.deltaPos.x, curConn.deltaPos.y, curConn.deltaPos.z], frac_coord , curConn.length, curConn.max_radius])
            channel["id"] = i
            channel["dim"] = c_channels[i].dimensionality
            channel["conns"] = conns
            channels.append(channel)
        return channels
    
    @classmethod
    def findChannels2(cls, vornet, atmnet, lower, upper, filename=None):
        cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
        cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
        cdef char* c_filename

        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef vector[CHANNEL] c_channels
        if not findChannels_new2(c_vornet_ptr, lower, upper, &c_channels):
            raise FindChannelError
        if filename:
            c_filename = filename
            if not c_writeToNET(c_channels, c_filename, c_atmnet_ptr):
                raise IOError

        channels = []
        cdef vector[int] nodeIDs
        cdef DELTA_POS disp
        cdef DIJKSTRA_NODE curNode
        cdef DIJKSTRA_NODE otherNode
        cdef CONN curConn
        for i in range(c_channels.size()):
            nodes = []
            conns = []
            channel = {}
           
            for j in range((c_channels[i].unitCells).size()):
                nodeIDs = (c_channels[i].ucNodes).at(j)
                disp = (c_channels[i].unitCells).at(j)
                for k in range(nodeIDs.size()):
                    curNode = (c_channels[i].nodes).at(nodeIDs.at(k))
                    node = {}
                    node["id"] = curNode.id
                    node["label"] = curNode.label
                    node["radius"] = curNode.max_radius
                    node["frac_coord"] = atmnet.absolute_to_relative(curNode.x, curNode.y, curNode.z)
                    nodes.append(node)
                    for l in range((curNode.connections).size()):
                        conn = {}
                        curConn = curNode.connections.at(l)
                        otherNode = (c_channels[i].nodes).at(curConn.ending)
                        frac_coord = atmnet.absolute_to_relative(curConn.btx, curConn.bty, curConn.btz)
                        conn["fromId"] = curNode.id
                        conn["fromLabel"] = curNode.label
                        conn["fromId"] = curNode.id
                        conn["fromDelta"] = [0, 0, 0]
                        conn["toId"] = otherNode.id
                        conn["toLabel"] = otherNode.label
                        conn["toDelta"] = [curConn.deltaPos.x, curConn.deltaPos.y, curConn.deltaPos.z]
                        conn["length"] = curConn.length
                        conn["bottleneck"] = frac_coord
                        conn["bottleneckSize"] = curConn.max_radius
                        conns.append(conn)
            channel["id"] = i
            channel["dim"] = c_channels[i].dimensionality
            channel["nodes"] = nodes
            channel["conns"] = conns
            channels.append(channel)
        return channels

    @classmethod
    def writeToVESTA(cls, channels, atmnet, filename):
        out = open(filename+".vesta","w")
        out.write("#VESTA_FORMAT_VERSION 3.3.0\n")
        out.write("\n")
        out.write("\n")
        out.write("CRYSTAL\n")
        out.write("\n")
        out.write("TITLE\n")
        out.write(filename + "\n")
        out.write("\n")
        out.write("GROUP\n")
        out.write("1 1 P 1\n")
        out.write("SYMOP\n")
        out.write(" 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1   1\n")
        out.write(" -1.0 -1.0 -1.0  0 0 0  0 0 0  0 0 0\n")
        out.write("TRANM 0\n")
        out.write(" 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1\n")
        out.write("LTRANSL\n")
        out.write(" -1\n")
        out.write(" 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000\n")
        out.write("LORIENT\n")
        out.write(" -1   0   0   0   0\n")
        out.write(" 1.000000  0.000000  0.000000  1.000000  0.000000  0.000000\n")
        out.write(" 0.000000  0.000000  1.000000  0.000000  0.000000  1.000000\n")
        out.write("LMATRIX\n")
        out.write(" 1.000000  0.000000  0.000000  0.000000\n")
        out.write(" 0.000000  1.000000  0.000000  0.000000\n")
        out.write(" 0.000000  0.000000  1.000000  0.000000\n")
        out.write(" 0.000000  0.000000  0.000000  1.000000\n")
        out.write(" 0.000000  0.000000  0.000000\n")
        out.write("CELLP\n")

        # write lattice parameters (a, b, c) and lattice angle (alpha, beta, gama).
        lattice_para = atmnet.lattice_para
        lattice_angle = atmnet.lattice_angle
        out.write(" " + str(round(lattice_para[0], 6)) + " " + str(round(lattice_para[1], 6)) + " " + str(round(lattice_para[2], 6)) + " " + 
            str(round(lattice_angle[0], 6)) + " " + str(round(lattice_angle[1], 6)) + " " + str(round(lattice_angle[2], 6)) + "\n")
        out.write(" 0.000000   0.000000   0.000000   0.000000   0.000000   0.000000\n")
        out.write("STRUC\n")

        # write Interstice parameters (a, b, c) and lattice angle (alpha, beta, gama).
        idx = 1
        for channel in channels:
            for node in channel["nodes"]:
                out.write(" " + str(idx) + " " + "He " + "He" + str(node["id"]) + " " + "1.0 " + 
                    str(round(node["frac_coord"][0], 6)) + " " +  str(round(node["frac_coord"][1], 6)) + " " + str(round(node["frac_coord"][2], 6)) + 
                    " 1a 1\n")
                out.write("                0.000000   0.000000   0.000000   0.00\n")
                idx = idx + 1
        bdx = 0
        for channel in channels:
            for conn in channel["conns"]:
                out.write(" " + str(idx) + " " + "Ne " + "Ne" + str(bdx) + " " + "1.0 " +
                   str(round(conn["bottleneck"][0], 6)) + " " + str(round(conn["bottleneck"][1], 6)) + " " + str(round(conn["bottleneck"][2], 6)) + 
                    " 1a 1\n")
                out.write("                0.000000   0.000000   0.000000   0.00\n")
                bdx = bdx + 1
                idx = idx + 1
        out.write("  0 0 0 0 0 0 0\n")
        out.write("THERI 0\n")

        count = 1
        for channel in channels:
            for node in channel["nodes"]:
                out.write(" " + str(count) + " " + "He " + "He" + str(node["id"]) + " 1.000000\n")
                count = count + 1
        bdx = 0
        for channel in channels:
            for conn in channel["conns"]:
                out.write(" " + str(count) + " " + "Ne " + "Ne" + str(bdx) + " 1.000000\n")
                count = count + 1
                bdx = bdx + 1
        out.write("  0 0 0\n")

        out.write("SHAPE\n")
        out.write("  0       0       0       0   0.000000  0   192   192   192   192\n")
        out.write("BOUND\n")
        out.write("       0        1         0        1         0        1\n")
        out.write("  0   0   0   0  0\n")
        out.write("SBOND\n")
        bond_count = 1
        for channel in channels:
            for conn in channel["conns"]:
                out.write(" " + str(bond_count) + " He" + str(conn["fromId"]) + " He" + str(conn["toId"]) + " " +
                    str(round(conn["length"] - 0.01, 6)) + " " + str(round(conn["length"] + 0.01, 6)) + 
                    " 0  0  0  1  1  0.200  1.000 8 0 148\n")
                bond_count = bond_count + 1
        out.write("  0 0 0 0\n")   
        out.write("SITET\n")
        idx = 1
        for channel in channels:
            for node in channel["nodes"]:
                out.write(" " + str(idx) + " " + "He" + str(node["id"]) + " " + str(round(node["radius"], 6)) + 
                    " 252 232 206 252 232 206 204  0\n")
                idx = idx + 1
        bdx = 0
        for channel in channels:
            for conn in channel["conns"]:
                out.write(" " + str(idx) + " " + "Ne" + str(bdx) + " " + str(round(conn["bottleneckSize"], 6)) +
                    " 254  55 181 254  55 181 204  0\n")
                idx = idx + 1
                bdx = bdx + 1
        out.write("  0 0 0 0 0 0\n")
        out.write("VECTR\n")
        out.write(" 0 0 0 0 0\n")
        out.write("VECTT\n")
        out.write(" 0 0 0 0 0\n")
        out.write("SPLAN\n")
        out.write("  0   0   0   0\n")
        out.write("LBLAT\n")
        out.write("-1\n")
        out.write("LBLSP\n")
        out.write("-1\n")
        out.write("DLATM\n")
        out.write("-1\n")
        out.write("DLBND\n")
        out.write("-1\n")
        out.write("DLPLY\n")
        out.write("-1\n")
        out.write("PLN2D\n")
        out.write("0   0   0   0\n")


# #Add at 20180823
# cdef class Channel_Node:
#     property id:
#         def __get__(self):
#             id = list(self.thisptr.id)
#             return id
#         def __set__(self, id):      # Don't set this
#             """
#             This variable is not supposed to be modified manually
#             """
#             print ("This value is not supposed to be modified")
#             self.thisptr.id = id
#     property coords:
#         def __get__(self):
#             coords = list(self.thisptr.x,self.thisptr.y,self.thisptr.z)
#             return coords
#         def __set__(self, coords):      # Don't set this
#             """
#             This variable is not supposed to be modified manually
#             """
#             print ("This value is not supposed to be modified")
#             self.thisptr.x = coords[0]
#             self.thisptr.y = coords[1]
#             self.thisptr.z = coords[2]
#     property radius:
#         def __get__(self): return self.thisptr.radius
#         def __set__(self, radius): 
#             print ("This value is not supposed to be modified")
#             self.thisptr.radius = radius
#             
# cdef class Channel_Edge:
#     property conn:
#         def __get__(self):
#             conn = [self.thisptr.from, self.thisptr.to]
#             return conn
#         def __set__(self, from, to):      # Don't set this
#             """
#             This variable is not supposed to be modified manually
#             """
#             print ("This value is not supposed to be modified")
#             self.thisptr.from = from
#             self.thisptr.to = to
# 
# cdef class Channel_new:
#     property lattice:
#         def __get__(self):
#             lattice = list(self.thisptr.lattice)
#             return lattice
#         def __set__(self, lattice):      # Don't set this
#             """
#             This variable is not supposed to be modified manually
#             """
#             print ("This value is not supposed to be modified")
#             self.thisptr.lattice = lattice
#     
