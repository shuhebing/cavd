# -*- coding: UTF-8 -*-
import os
import sys
import re
from cavd.netstorage import AtomNetwork
from cavd.netstorage import connection_values_list
from cavd.channel import Channel
from cavd.area_volume import asa_new
from cavd.netio import *
from cavd.local_environment import get_local_envir
from cavd.get_Symmetry import get_Symmetry
from cavd.high_accuracy import high_accuracy_atmnet
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator

#from local_environment import get_local_envir,Coordination


# 获取特定结构中
# 所有离子的有效半径，
# 迁移离子到最邻近离子表面距离与迁移离子半径的比值alpha，
# 迁移离子半径
def LocalEnvirCom(filename, migrant):
    # stru = Structure.from_file(filename)
    # val_eval = ValenceIonicRadiusEvaluator(stru)
    # radii = val_eval.radii
    coordination_list, radii = get_local_envir(filename)
    
    # 为了防止保存的半径信息无法匹配此处对半径信息做特殊处理，如Ag+的半径会保存为Ag、Ag+、Ag1+
    # radii_keys = list(radii.keys())
    # for key in radii_keys:
    #     radii[re.sub('[^a-zA-Z]','',key)] = radii[key]
    #     if re.search('[A-Z][a-z]*\+|[A-Z][a-z]*\-', key) != None:
    #         s1 = re.sub('\+','1+',key)
    #         radii[re.sub('\-','1-',s1)] = radii[key]

    # minRad = 0
    # for label in radii:
    #     if migrant in label:
    #         # 取最小值作为迁移离子半径
    #         if minRad == 0:
    #             minRad = radii[label]
    #         elif radii[label] < minRad:
    #             minRad = radii[label]
    
    migrant_labels = []
    coord_tmp = []
    nei_dis_tmp = []
    min_nei_dis_tmp = []
    migrant_radii_tmp = []
    migrant_paras_tmp = []
    for i in coordination_list:
        if migrant in i["label"]:
            if len(i["coord_nei"]) == 0:
                nei = i["coord_nei"][0]
                nei_dis = nei[1]

                nei_label = nei[0]
                nei_radius = radii[nei_label]

                migrant_radius = i["radius"]
                migrant_label = i["label"]

                alpha = (nei_dis - nei_radius)/migrant_radius

                migrant_labels.append(migrant_label)
                migrant_radii_tmp.append(migrant_radius)
                migrant_paras_tmp.append(alpha)
                coord_tmp.append(i["coord_num"])
                nei_dis_tmp.append(nei_dis)
                min_nei_dis_tmp.append(nei_dis - nei_radius)

    if len(migrant_labels) != 0 and len(migrant_paras) != 0 and len(nei_dises) != 0:
        migrant_radii = dict(zip(migrant_labels, migrant_radii_tmp))
        migrant_paras = dict(zip(migrant_labels, migrant_paras_tmp))
        nei_dises = dict(zip(coord_tmp, zip(nei_dis_tmp, min_nei_dis_tmp)))

        rad_sum = 0
        alpha_sum = 0
        for value in migrant_radii.values():
            rad_sum += value
        for value in migrant_paras.values():
            alpha_sum += value
        migrant_radius = rad_sum/len(migrant_radii)
        migrant_alpha = alpha_sum/len(migrant_paras)
        # print(migrant_radii)
        # print(migrant_paras)
        print(migrant_radius)
        print(migrant_alpha)
        print(migrant+" radii: ",radii)
        return radii,migrant_radius,migrant_alpha,nei_dises
    else:
        return radii,migrant_radius,migrant_alpha,nei_dises

# 获取特定结构中
# 所有离子的有效半径
def getIonicRadii(filename):
    coordination_list, radii = get_local_envir(filename)
    print(radii)
    return radii

def AllCom(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True, minRad=0.0, maxRad=0.0):
    radii = {}
    if rad_flag and effective_rad:
        #考虑如何利用migrant_radius与migrant_alpha
        radii,migrant_radius,migrant_alpha, nei_dises = LocalEnvirCom(filename,migrant)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    # high_accur_atmnet = atmnet.copy()
    # high_accuracy_atmnet(high_accur_atmnet, "S50")

    if migrant:
        os.remove(remove_filename)

    prefixname = filename.replace(".cif","")
    # vornet,edge_centers,fcs = high_accur_atmnet.perform_voronoi_decomposition(False)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)

    writeBIFile(prefixname+".bi",atmnet,vornet)

    sym_vornet,voids = get_Symmetry(atmnet, vornet)

    writeBIFile(prefixname+"_orgin.bi",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_orgin.vasp",atmnet,sym_vornet,rad_store_in_vasp)

    probe_rad = migrant_radius*migrant_alpha
    minRad = migrant_radius*migrant_alpha*0.85
    maxRad = migrant_radius*migrant_alpha*1.15
    print(minRad)
    print(maxRad)

    writeVaspFile(prefixname+"_selected.vasp",atmnet,sym_vornet,rad_store_in_vasp,minRad,maxRad)
    conn = connection_values_list(prefixname+".resex", sym_vornet)

    channels = Channel.findChannels(sym_vornet,atmnet,0.60,prefixname+".net")
    oneD,twoD,threeD = ConnStatus(minRad, conn)
    
    dims = []
    for i in channels:
        dims.append(i["dim"])
    return conn,oneD,twoD,threeD,nei_dises,dims,voids

# 使用带半径的公式进行计算
def AllCom5(filename, standard, migrant=None, rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True):
    radii = {}
    if rad_flag and effective_rad:
        #考虑如何利用migrant_radius与migrant_alpha
        radii, migrant_radius, migrant_alpha, nei_dises = LocalEnvirCom(filename,migrant)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    # high_accur_atmnet = atmnet.copy()
    # high_accuracy_atmnet(high_accur_atmnet, "S50")
    if migrant:
        os.remove(remove_filename)

    prefixname = filename.replace(".cif","")
    # vornet,edge_centers,fcs = high_accur_atmnet.perform_voronoi_decomposition(False)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    sym_vornet,voids = get_Symmetry(atmnet, vornet)

    writeBIFile(prefixname+"_orgin.bi",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_orgin.vasp",atmnet,sym_vornet,rad_store_in_vasp)

    minRad = standard*migrant_alpha*0.85
    print(minRad)

    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+".net")
    
    dims = []
    for i in channels:
        dims.append(i["dim"])
    return nei_dises,dims,voids

# 使用簇替换的方法进行计算
def AllCom4(filename, standard, migrant=None, rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True):
    radii = {}
    if rad_flag and effective_rad:
        #考虑如何利用migrant_radius与migrant_alpha
        radii, migrant_radius, migrant_alpha, nei_dises = LocalEnvirCom(filename,migrant)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    high_accur_atmnet = atmnet.copy()
    high_accuracy_atmnet(high_accur_atmnet, "S50")
    if migrant:
        os.remove(remove_filename)

    prefixname = filename.replace(".cif","")
    vornet,edge_centers,fcs = high_accur_atmnet.perform_voronoi_decomposition(False)
    sym_vornet,voids = get_Symmetry(high_accur_atmnet, vornet)

    writeBIFile(prefixname+"_orgin.bi",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_orgin.vasp",atmnet,sym_vornet,rad_store_in_vasp)

    minRad = standard*migrant_alpha*0.85
    print(minRad)

    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+".net")
    dims = []
    for i in channels:
        dims.append(i["dim"])
    return nei_dises,dims,voids

# 使用不带半径的公式进行计算
def AllCom3(filename, standard, migrant=None, rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True):
    radii = {}
    if rad_flag and effective_rad:
        #考虑如何利用migrant_radius与migrant_alpha
        radii, migrant_radius, migrant_alpha, nei_dises = LocalEnvirCom(filename,migrant)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    rad_flag = False
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    # high_accur_atmnet = atmnet.copy()
    # high_accuracy_atmnet(high_accur_atmnet, "S50")
    if migrant:
        os.remove(remove_filename)

    prefixname = filename.replace(".cif","")
    # vornet,edge_centers,fcs = high_accur_atmnet.perform_voronoi_decomposition(False)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    sym_vornet,voids = get_Symmetry(atmnet, vornet)

    writeBIFile(prefixname+"_orgin.bi",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_orgin.vasp",atmnet,sym_vornet,rad_store_in_vasp)

    minRad = standard*migrant_alpha*0.85
    print(minRad)

    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+".net")
    
    dims = []
    for i in channels:
        dims.append(i["dim"])
    return nei_dises,dims,voids

#AllCom
def AllCom2(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True, minRad=0.0, maxRad=0.0):
    radii = {}
    if rad_flag and effective_rad:
        radii = getIonicRadii(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    # high_accur_atmnet = atmnet.copy()
    # high_accuracy_atmnet(high_accur_atmnet, "S50")

    if migrant:
        os.remove(remove_filename)

    prefixname = filename.replace(".cif","")
    # vornet,edge_centers,fcs = high_accur_atmnet.perform_voronoi_decomposition(False)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)

    sym_vornet,voids = get_Symmetry(atmnet, vornet)

    writeBIFile(prefixname+"_orgin.bi",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_orgin.vasp",atmnet,sym_vornet,rad_store_in_vasp)

    writeVaspFile(prefixname+"_selected.vasp",atmnet,sym_vornet,rad_store_in_vasp,minRad,maxRad)

    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+".net")
    
    dims = []
    for i in channels:
        dims.append(i["dim"])
    return dims,voids

#计算指定结构的瓶颈和间隙
def BIComputation(filename, migrant=None, rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True,  minRad=0.0, maxRad=0.0):
    radii = {}
    if rad_flag and effective_rad:
        radii,migrant_radius,migrant_alpha, nei_dises = LocalEnvirCom(filename,migrant)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    #delete temp file
    if migrant:
        os.remove(remove_filename)
    prefixname = filename.replace(".cif","")
    writeBIFile(prefixname+"_orgin.bi",atmnet,vornet)
    writeBIFile(prefixname+"_selected.bi",atmnet,vornet)
    writeVaspFile(prefixname+"_orgin.vasp",atmnet,vornet,rad_store_in_vasp)
    writeVaspFile(prefixname+"_selected.vasp",atmnet,vornet,rad_store_in_vasp,minRad,maxRad)

#计算指定结构最大自由球体半径，最大包含球体半径和沿着最大自由球体路径上的最大包含球体半径：Rf Ri Rif
def ConnValCom(filename, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = LocalEnvirCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    if migrant:
        os.remove(remove_filename)
    prefixname = filename.replace(".cif","")
    Ri,Rf,Rif = atmnet.through_VorNet(prefixname+".res")
    return Ri,Rf,Rif
    
#计算某个结构的连通数值列表，存放a，b，c方向上的Rf
def ConnValListCom(filename, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = LocalEnvirCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    if migrant:
        os.remove(remove_filename)
    prefixname = filename.replace(".cif","")
    conn = connection_values_list(prefixname+".resex",vornet)
    return conn

#判断某个结构的连通性,给定目标离子的半径，判断它是否是1D，2D，3D导通
def ConnStatusCom(filename, radius, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    connlist = ConnValListCom(filename, migrant, rad_flag, effective_rad, rad_file)
    oneD = False
    twoD = False
    threeD = False

    af = connlist[0]
    bf = connlist[1]
    cf = connlist[2]
    
    if radius <= af:
        aconn = True
    else:
        aconn = False
    
    if radius <= bf:
        bconn = True
    else:
        bconn = False
    
    if radius <= cf:
        cconn = True
    else:
        cconn = False
    
    if aconn and bconn and cconn:
        threeD = True
    if (aconn and bconn) or (aconn and cconn) or (bconn and cconn):
        twoD = True
    if aconn or bconn or cconn:
        oneD = True
    return oneD,twoD,threeD

#根据连通数值列表，判断某个结构的连通性。给定一个原子的半径，判断它是否是1D，2D，3D导通
def ConnStatus(radius,connlist):
    oneD = False
    twoD = False
    threeD = False

    af = connlist[0]
    bf = connlist[1]
    cf = connlist[2]
    
    if radius <= af:
        aconn = True
    else:
        aconn = False
    
    if radius <= bf:
        bconn = True
    else:
        bconn = False
    
    if radius <= cf:
        cconn = True
    else:
        cconn = False
    
    if aconn and bconn and cconn:
        threeD = True
    if (aconn and bconn) or (aconn and cconn) or (bconn and cconn):
        twoD = True
    if aconn or bconn or cconn:
        oneD = True
    return oneD,twoD,threeD
    
#计算通道
def ChannelCom(filename, probe_rad, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = LocalEnvirCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    if migrant:
        os.remove(remove_filename)
    prefixname = filename.replace(".cif","")
    Channel.findChannels(vornet,atmnet,probe_rad,prefixname+".net")


#计算ASA
def ASACom(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = LocalEnvirCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    if migrant:
        os.remove(remove_filename)
    prefixname = filename.replace(".cif","")
    asa_new(prefixname+".zsa",False,atmnet,probe_rad,probe_rad,num_sample)

#计算空隙网络
def VoidNetCom(filename, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = LocalEnvirCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    if migrant:
        os.remove(remove_filename)
    prefixname = filename.replace(".cif","")
    writeZVisFile(prefixname+".zvis", False, atmnet, vornet)


# if __name__ == "__main__":
#     # radii = LocalEnvirCom("../../examples/icsd_16713.cif")
#     conn,oneD,twoD,threeD = AllCom("../../examples/Li2CO3-LDA.cif",0.5,1000,"Li", True,True,None,True,0.5,0.7)