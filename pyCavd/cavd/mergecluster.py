import numpy as np
from scipy.spatial import cKDTree
import networkx as nx
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure


class Void(object):
    # 把间隙抽象为Void类，包含id、类别、分数坐标、半径、能量等属性
    def __init__(self):
        self.id = None
        self.label = None
        self.coord = None
        self.radii = None
        self.energy = None


class Channel(object):
    # 把通道抽象为Channel类，包含开始间隙点id、结束间隙点id、晶相、瓶颈分数坐标、瓶颈尺寸等属性
    def __init__(self):
        self.start = None
        self.end = None
        self.phase = None
        self.coord = None
        self.radii = None
        self.dist = None
        self.label = None


class MergeCluster(object):
    def __init__(self):
        # threshold为簇之间的阈值
        self.__struc = None
        self.voids = None
        self.channels = None
        self.mergedvoids = []
        self.mergedchannels = []
        self.threshold = 0.2
        self.clusters = []

    def setstruc(self, filename_cif):
        self.__struc = Structure.from_file(filename_cif)

    def setvoids(self, voids_list):
        self.voids = voids_list

    def setchannels(self, channels_list):
        self.channels = channels_list

    @property
    def voids_mc(self):
        return self.voids

    @property
    def channels_mc(self):
        return self.channels

    def get_absolute_dis(self, p1, p2):
        """
        在不考虑周期性的情况下，计算两点之间的距离
        :param p1: 分数坐标，例如[0.5，0.5，0.5]
        :param p2: 分数坐标
        :return: 两点之间的欧式距离，类型为float
        """
        coord1 = np.array(self.fac2cart(p1))
        coord2 = np.array(self.fac2cart(p2))
        diff = np.subtract(coord1, coord2)
        return np.linalg.norm(diff)

    def get_period_dis(self, p1, p2):
        """
        在考虑周期性的情况下，计算两点之间的距离
        :param p1: 分数坐标，例如[0.5，0.5，0.5]
        :param p2: 分数坐标
        :return:  两点之间的距离，考虑周期性
        """
        temp_site1 = PeriodicSite('Ar', p1, self.__struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self.__struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis

    @staticmethod
    def readdata(filename_cavd):
        # 读取CAVD方法计算出的连通文件bi
        voids_list = []
        channels_list = []
        flag_p = 0
        flag_n = 0
        file = open(filename_cavd, 'r')
        for line in file.readlines():
            if 'Interstitial' in line:
                flag_p = 1
                flag_n = 0
                continue
            if 'Channel' in line:
                flag_p = 0
                flag_n = 1
                continue
            if flag_p == 1:
                line = line.split()
                if len(line) > 3:
                    void = Void()
                    void.id = int(line[0])
                    void.label = int(line[1])
                    void.coord = [np.float64(line[2]), np.float64(line[3]), np.float64(line[4])]
                    void.radii = np.float64(line[5])
                    voids_list.append(void)
            if flag_n == 1:
                line = line.split()
                if len(line) > 4:
                    channel = Channel()
                    channel.start = int(line[0])
                    channel.end = int(line[1])
                    channel.phase = [int(line[2]), int(line[3]), int(line[4])]
                    channel.coord = [np.float64(line[5]), np.float64(line[6]), np.float64(line[7])]
                    channel.radii = np.float64(line[8])
                    channels_list.append(channel)
        return voids_list, channels_list

    def fac2cart(self, coord):
        """
        分数坐标转换成笛卡尔坐标
        """
        return np.dot(coord, self.__struc.lattice.matrix)

    def cart2fac(self, coord):
        """
        笛卡尔坐标转换成分数坐标
        """
        return np.dot(coord, np.linalg.inv(self.__struc.lattice.matrix))

    def cal_clusters(self):
        # 找到每个簇，结果保存为字典
        coords = []
        displacevetors = np.array([[0, 0, 0], [1, 1, 0], [1, 1, -1], [1, 0, 1], [1, 0, 0], [1, 0, -1], [1, -1, 1],
                                  [1, -1, 0], [1, -1, -1], [0, 1, 1], [0, 1, 0], [0, 1, -1], [0, 0, 1], [1, 1, 1]])
        for dv in displacevetors:
            for void in self.voids:
                coords.append(self.fac2cart(void.coord+dv))
        print(len(self.voids))
        coord_tree = cKDTree(coords)
        pair_void1 = [i for i in coord_tree.query_pairs(r=self.threshold)
                      if i[0] < len(self.voids) or i[1] < len(self.voids)]
        pair_void = []
        for i in pair_void1:
            if i[0] < len(self.voids) and i[1] < len(self.voids):
                pair_void.append([i[0], i[1]])
            elif i[0] < len(self.voids):
                pair_void.append([i[0], i[1] % len(self.voids)])
            else:
                pair_void.append([i[1], i[0] % len(self.voids)])
        if len(pair_void) > 0:
            id_void = {}
            for i in range(len(self.voids)):
                id_void[i] = self.voids[i]
            graph_clusters = nx.Graph()
            for e in pair_void:
                graph_clusters.add_edge(e[0], e[1])
            queue_clusters = []
            for sc in nx.connected_component_subgraphs(graph_clusters):
                queue_clusters.append(list(sc.nodes))
            while queue_clusters:
                subv_in = queue_clusters.pop(0)
                subv_out = []
                centre_coord = self.fac2cart(id_void[subv_in[0]].coord)
                for i in range(1, len(subv_in)):
                    disvector = np.around(np.array(self.cart2fac(centre_coord)) - np.array(id_void[subv_in[i]].coord))
                    tempcoord = self.fac2cart([id_void[subv_in[i]].coord[0] + disvector[0],
                                               id_void[subv_in[i]].coord[1] + disvector[1],
                                               id_void[subv_in[i]].coord[2] + disvector[2], ])
                    centre_coord[0] = (tempcoord[0] + centre_coord[0]) / 2.0
                    centre_coord[1] = (tempcoord[1] + centre_coord[1]) / 2.0
                    centre_coord[2] = (tempcoord[2] + centre_coord[2]) / 2.0
                centre_coord = self.cart2fac(centre_coord)
                for i in centre_coord:
                    if i < 0:
                        i += 1
                    if i > 1:
                        i -= 1
                for i in range(len(subv_in) - 1, -1, -1):
                    if self.get_period_dis(id_void[subv_in[i]].coord, centre_coord) > self.threshold:
                        subv_out.append(subv_in[i])
                        subv_in.remove(subv_in[i])
                if len(subv_out) == 0:
                    subid_in = []
                    for si in subv_in:
                        subid_in.append(id_void[si].id)
                    self.clusters.append({"voidsid": subid_in, "centre_coord": centre_coord})
                else:
                    centre_coord = self.fac2cart(id_void[subv_in[0]].coord)
                    for i in range(1, len(subv_in)):
                        disvector = np.around(
                            np.array(self.cart2fac(centre_coord)) - np.array(id_void[subv_in[i]].coord))
                        tempcoord = self.fac2cart([id_void[subv_in[i]].coord[0] + disvector[0],
                                                   id_void[subv_in[i]].coord[1] + disvector[1],
                                                   id_void[subv_in[i]].coord[2] + disvector[2], ])
                        centre_coord[0] = (self.fac2cart(tempcoord)[0] + centre_coord[0]) / 2.0
                        centre_coord[1] = (self.fac2cart(tempcoord)[1] + centre_coord[1]) / 2.0
                        centre_coord[2] = (self.fac2cart(tempcoord)[2] + centre_coord[2]) / 2.0
                    centre_coord = self.cart2fac(centre_coord)
                    for i in centre_coord:
                        if i < 0:
                            i += 1
                        if i > 1:
                            i -= 1
                    subid_in = []
                    for si in subv_in:
                        subid_in.append(id_void[si].id)
                    self.clusters.append({"voidsid": subid_in, "centre_coord": centre_coord})
                    pair_subvout = []
                    for i in range(len(subv_out)):
                        for j in range(i+1, len(subv_out)):
                            if self.get_period_dis(id_void[subv_out[i]].coord,
                                                   id_void[subv_out[j]].coord) < self.threshold:
                                pair_subvout.append([subv_out[i], subv_out[j]])
                    if len(pair_subvout) > 0:
                        graph_subvout = nx.Graph()
                        for e in pair_subvout:
                            graph_subvout.add_edge(e[0], e[1])
                        for sc in nx.connected_component_subgraphs(graph_subvout):
                            queue_clusters.append(list(sc.nodes))
        print(self.clusters)


    def process_voidsandchannels(self):
        mignet = nx.Graph()
        for n in self.voids:  # 添加点
            mignet.add_node(n.id, label=n.label, coord=n.coord, radii=n.radii)
        for e in self.channels:  # 添加边
            if e.start < e.end:
                l1 = e.phase
                l2 = [-1 * i for i in e.phase]
            else:
                l1 = [-1 * i for i in e.phase]
                l2 = e.phase
            mignet.add_edge(e.start, e.end, phase1=l1, phase2=l2, coord1=e.coord, radii1=e.radii)

        if len(self.clusters) > 0:
            id_void = {}
            for i in range(len(self.voids)):
                id_void[self.voids[i].id] = self.voids[i]

            for i in range(len(self.clusters)):
                tempvoid = Void()   # tempvoid为添加的新点
                tempvoid.label = 1000
                tempvoid.coord = self.clusters[i]['centre_coord']
                tempvoid.radii = id_void[self.clusters[i]['voidsid'][0]].radii
                for void in self.clusters[i]['voidsid']:
                    if id_void[void].label < tempvoid.label:
                        tempvoid.label = id_void[void].label
                        tempvoid.id = id_void[void].id
                    tempvoid.radii = (tempvoid.radii + id_void[void].radii) / 2.0
                tempedges = []
                nearedges = []
                for id in self.clusters[i]['voidsid']:
                    for nearvoid in list(mignet.adj[id].keys()):
                        if nearvoid not in self.clusters[i]['voidsid']:
                            nearedges.append([id, nearvoid])
                def takeSecond(elem):
                    return elem[1]
                nearedges.sort(key=takeSecond)
                for j in range(len(nearedges)):
                    if nearedges[j][1] != nearedges[j-1][1] or j == 0:
                        if nearedges[j][0] <nearedges[j][1]:
                            disvector = np.around(
                                np.array(tempvoid.coord) - np.array(id_void[nearedges[j][0]].coord)) \
                                        +  np.array(mignet[nearedges[j][0]][nearedges[j][1]]['phase1'])
                            if tempvoid.id < nearedges[j][1]:
                                p1 = disvector
                                p2 = [-1 * i for i in disvector]
                            else:
                                p2 = disvector
                                p1 = [-1 * i for i in disvector]
                        else:
                            disvector = np.around(
                                np.array(tempvoid.coord) - np.array(id_void[nearedges[j][0]].coord)) \
                                        + np.array(mignet[nearedges[j][0]][nearedges[j][1]]['phase2'])
                            if tempvoid.id < nearedges[j][1]:
                                p1 = disvector
                                p2 = [-1 * i for i in disvector]
                            else:
                                p2 = disvector
                                p1 = [-1 * i for i in disvector]


                        tempedges.append({"from": tempvoid.id, "to": nearedges[j][1],
                                          "phase1": p1,
                                          "phase2": p2,
                                          "coord1": mignet[nearedges[j][0]][nearedges[j][1]]['coord1'],
                                          "radii1": mignet[nearedges[j][0]][nearedges[j][1]]['radii1']})
                    else:
                        if mignet[nearedges[j][0]][nearedges[j][1]]['radii1'] > tempedges[-1]['radii1']:
                            tempedges.pop()
                            if nearedges[j][0] < nearedges[j][1]:
                                disvector = np.around(
                                    np.array(tempvoid.coord) - np.array(id_void[nearedges[j][0]].coord)) \
                                            + np.array(mignet[nearedges[j][0]][nearedges[j][1]]['phase1'])
                                if tempvoid.id < nearedges[j][1]:
                                    p1 = disvector
                                    p2 = [-1 * i for i in disvector]
                                else:
                                    p2 = disvector
                                    p1 = [-1 * i for i in disvector]
                            else:
                                disvector = np.around(
                                    np.array(tempvoid.coord) - np.array(id_void[nearedges[j][0]].coord)) \
                                            + np.array(mignet[nearedges[j][0]][nearedges[j][1]]['phase2'])
                                if tempvoid.id < nearedges[j][1]:
                                    p1 = disvector
                                    p2 = [-1 * i for i in disvector]
                                else:
                                    p2 = disvector
                                    p1 = [-1 * i for i in disvector]
                            tempedges.append({"from": tempvoid.id, "to": nearedges[j][1],
                                              "phase1": p1,
                                              "phase2": p2,
                                              "coord1": mignet[nearedges[j][0]][nearedges[j][1]]['coord1'],
                                              "radii1": mignet[nearedges[j][0]][nearedges[j][1]]['radii1']})
                for void in self.clusters[i]['voidsid']:
                    mignet.remove_node(void)
                mignet.add_node(tempvoid.id, label=tempvoid.label, coord=tempvoid.coord, radii=tempvoid.radii)
                for e in tempedges:
                    mignet.add_edge(e['from'], e['to'], label=e['from'], phase1=e['phase1'], phase2=e['phase2'],
                                    coord1=e['coord1'], radii1=e['radii1'])

        for node in mignet.nodes():
            tempvoid = Void()
            tempvoid.id = node
            tempvoid.label = mignet.node[node]['label']
            tempvoid.coord = mignet.node[node]['coord']
            tempvoid.radii = mignet.node[node]['radii']
            self.mergedvoids.append(tempvoid)
        for edge in mignet.edges():
            tempchannel1 = Channel()
            tempchannel2 = Channel()
            tempchannel1.start = edge[0]
            tempchannel1.end = edge[1]
            tempchannel2.end = edge[0]
            tempchannel2.start = edge[1]
            if edge[0] < edge[1]:
                tempchannel1.phase = mignet[edge[0]][edge[1]]["phase1"]
                tempchannel2.phase = mignet[edge[0]][edge[1]]["phase2"]
            else:
                tempchannel1.phase = mignet[edge[0]][edge[1]]["phase2"]
                tempchannel2.phase = mignet[edge[0]][edge[1]]["phase1"]
            tempchannel1.coord = mignet[edge[0]][edge[1]]["coord1"]
            tempchannel2.coord = mignet[edge[0]][edge[1]]["coord1"]
            tempchannel1.radii = mignet[edge[0]][edge[1]]["radii1"]
            tempchannel2.radii = mignet[edge[0]][edge[1]]["radii1"]
            self.mergedchannels.append(tempchannel1)
            self.mergedchannels.append(tempchannel2)

    def savedata(self, filename_cif):
        with open(filename_cif.split(".")[0]+'_mergecluster_network.net', 'w') as f:
            f.write('Interstitial table:\n')
            for void in self.mergedvoids:
                f.write(str(void.id)+"\t"+str(void.label)+"\t "
                        + str(void.coord[0]) + " " + str(void.coord[1]) + " "+str(void.coord[2]) + "\t "
                        + str(void.radii)
                        + "\n")
            f.write('Connection table::\n')
            for channel in self.mergedchannels:
                f.write(str(channel.start) + "\t " + str(channel.end) + "\t " + str(channel.phase[0]) + " "
                        + str(channel.phase[1]) + " " + str(channel.phase[2]) + "\t "
                        + str(channel.coord[0]) + " "
                        + str(channel.coord[1]) + " " + str(channel.coord[2]) + "\t "
                        + str(channel.radii) + "\n")


class CalChannelLabel(object):
    def __init__(self):
        self.__voids = []
        self.__channels = []
        self.nonequalchannels = {}
        self.voidid_label = {}

    def setvoids(self, voids_list):
        self.__voids = voids_list

    def setchannels(self, channels_list):
        self.__channels = channels_list

    @property
    def voids(self):
        return self.__voids

    @property
    def channels(self):
        return self.__channels

    @staticmethod
    def readdata(filename_cavd):
        # 读取CAVD方法计算出的连通文件bi
        voids_list = []
        channels_list = []
        flag_p = 0
        flag_n = 0
        file = open(filename_cavd, 'r')
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
                line = line.split()
                if len(line) > 3:
                    void = Void()
                    void.id = int(line[0])
                    void.label = int(line[1])
                    void.coord = [np.float64(line[2]), np.float64(line[3]), np.float64(line[4])]
                    void.radii = np.float64(line[5])
                    void.energy = np.float64(line[6])
                    voids_list.append(void)
            if flag_n == 1:
                line = line.split()
                if len(line) > 4:
                    channel = Channel()
                    channel.start = int(line[0])
                    channel.end = int(line[1])
                    channel.phase = [int(line[2]), int(line[3]), int(line[4])]
                    channel.coord = [np.float64(line[5]), np.float64(line[6]), np.float64(line[7])]
                    channel.radii = np.float64(line[8])
                    channel.dist = np.float64(line[9])
                    channels_list.append(channel)
        return voids_list, channels_list

    def calnonequalchannels(self, filename_cavd, preci=3):

        for void in self.__voids:
            self.voidid_label[void.id] = void.label
        i = 0
        for channel in self.__channels:
            key1 = (self.voidid_label[channel.start],  self.voidid_label[channel.end],
                    round(channel.dist, preci),round(channel.radii, preci))
            key2 = (self.voidid_label[channel.end],  self.voidid_label[channel.start],
                    round(channel.dist, preci), round(channel.radii, preci))
            if key1 not in self.nonequalchannels.keys() and key2 not in self.nonequalchannels.keys():
                self.nonequalchannels[key1] = i
                i += 1
        print(self.nonequalchannels)

    def calchannellabel(self):
        for channel in self.__channels:
            if (self.voidid_label[channel.start], self.voidid_label[channel.end], round(channel.dist, 3),
                                        round(channel.radii, 3)) in self.nonequalchannels.keys():
                key_nonequalchannel = (self.voidid_label[channel.start], self.voidid_label[channel.end],
                                        round(channel.dist, 3), round(channel.radii, 3))
            elif (self.voidid_label[channel.end], self.voidid_label[channel.start], round(channel.dist, 3),
                                        round(channel.radii, 3)) in self.nonequalchannels.keys():
                key_nonequalchannel = (self.voidid_label[channel.end], self.voidid_label[channel.start],
                                        round(channel.dist, 3), round(channel.radii, 3))
            else:
                raise ValueError("非等价路径片段计算错误")
            channel.label = self.nonequalchannels[key_nonequalchannel]

    def savedata(self, filename_cif):
        with open(filename_cif.split(".")[0]+'_calchannellabel_network.net', 'w') as f:
            f.write('Interstitial table:\n')
            for void in self.__voids:
                f.write(str(void.id)+"\t"+str(void.label)+"\t "
                        + str(void.coord[0]) + " " + str(void.coord[1]) + " "+str(void.coord[2]) + "\t "
                        + str(void.radii) + "\t " + str(void.energy) + "\n")
            f.write('Connection table::\n')
            for channel in self.__channels:
                f.write(str(channel.start) + "\t " + str(channel.end) + "\t " + str(channel.phase[0]) + " "
                        + str(channel.phase[1]) + " " + str(channel.phase[2]) + "\t "
                        + str(channel.coord[0]) + " "
                        + str(channel.coord[1]) + " " + str(channel.coord[2]) + "\t "
                        + str(channel.radii) + "\t "
                        + str(channel.dist) + "\t "
                        + str(channel.label) + "\n")


if __name__ == "__main__":
    filename_CIF = 'stucture/icsd_000467/icsd_000467.cif'
    filename_CAVD = "stucture/icsd_000467/icsd_000467_origin.net"
    mc = MergeCluster()
    mc.setstruc(filename_CIF)
    voids, channels = mc.readdata(filename_CAVD)
    mc.setvoids(voids)
    mc.setchannels(channels)
    mc.cal_clusters()
    mc.process_voidsandchannels()
    mc.savedata(filename_CIF)




