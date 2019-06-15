# -*- coding: UTF-8 -*-
import os
import sys
import re
from cavd.netstorage import AtomNetwork
from cavd.netstorage import connection_values_list
from cavd.channel import Channel
from cavd.area_volume import asa_new
from cavd.netio import *
from cavd.local_environment import get_local_envir_fromstru
from cavd.get_Symmetry import get_symnum_sites, get_equivalent_vornet,get_labeled_vornet
from cavd.high_accuracy import high_accuracy_atmnet
from cavd.local_environment import CifParser_new
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator
from pymatgen.core.sites import PeriodicSite
from monty.io import zopen
import numpy as np
from scipy.spatial.ckdtree import cKDTree

#from local_environment import get_local_envir,Coordination


# 计算迁移离子晶格位的复原率
def rediscovery(migrate,vorosites,stru):
    labels = []
    recover_labels = []
    recover_state = {}
    true_recover_dis = {}
    #点类型，分别表示间隙、瓶颈、面心
    points_type = ["It","Bn","Fc"]


    for k in range(len(stru.sites)):
        site = stru.sites[k]
        label = site._atom_site_label
        if migrate not in label:
            continue
        #labels记录所有的label
        if label not in labels:
            labels.append(label)
        #recover_labels记录已恢复的label
        if label in recover_labels:
            continue

        for pts_idx, pts in enumerate (vorosites):
            cp_tag = np.ones((len(pts), ), dtype=int)
            for pt_idx, pt in enumerate (pts):
                if cp_tag[pt_idx] != -1:
                    print("mobile:",site,"label",label)
                    print("void:",pt)
                    #以Ar作为临时当前空隙的表示符号
                    tmp_site = PeriodicSite("Ar",pt,stru.lattice)
                    print(site.distance(tmp_site))

                    if site.distance(tmp_site) < 0.5:
                        #当某空隙位已与结构中的晶格位配对时，将该空隙位以及与该空隙位0.25A半径范围内的所有空隙位移除，后续不再判断。
                        recover_labels.append(label)

                        true_recover_dis[str(label)] = (points_type[pts_idx]+str(pt_idx),site.distance(tmp_site))
                        cp_tag[pt_idx] = -1
                        
                        for pt_idx2, pt2 in enumerate (pts):
                            tmp_site2 = PeriodicSite("Ar",list(pt2),stru.lattice)
                            if tmp_site.distance(tmp_site2) < 0.25:
                                cp_tag[pt_idx2] = -1
                        break

    #统计当前结构的恢复率
    recover_rate = len(recover_labels)/len(labels)
    for la in labels:
        if la in recover_labels:
            recover_state[str(la)] = True
        else:
            recover_state[str(la)] = False

    return recover_rate, recover_state, true_recover_dis

# 返回点的tag
def get_point_tag(id, pts_len):
    vexs_len = pts_len[0]
    bts_len = pts_len[1]
    fcs_len = pts_len[2]
    #点类型，分别表示间隙、瓶颈、面心
    if id < vexs_len:
        return "It" + str(id)
    elif id < vexs_len + bts_len:
        return "Bn" + str(id - vexs_len)
    elif id < vexs_len + bts_len + fcs_len:
        return "Fc" + str(id - vexs_len - bts_len)
    else:
        raise IndexError

# 找晶格位
def rediscovery_kdTree(migrate,vorosites,stru):
    recover_labels = []
    recover_state = {}
    migrate_mindis = {}
    
    migrate_pos_frac = np.around(np.array([site.frac_coords for site in stru.sites if migrate in site._atom_site_label], ndmin=2), 3)
    # print(migrate_pos_frac)
    migrate_pos_frac %= 1.0
    migrate_pos_frac %= 1.0
    # print(migrate_pos_frac)
    migrate_pos = [site.coords for site in stru.sites if migrate in site._atom_site_label]
    # print(migrate_pos)
    labels = [site._atom_site_label for site in stru.sites if migrate in site._atom_site_label]

    points = np.around(np.array(vorosites[0] + vorosites[1] + vorosites[2], ndmin=2), 3)
    # print(points)
    points %= 1.0
    points %= 1.0
    # print(points)
    # print(len(points))
    vorositesKdTree = cKDTree(points)
    min_dis,min_ids = vorositesKdTree.query(migrate_pos_frac,k=1)
    # print(labels)
    # print(min_dis)

    for idx in range(len(min_ids)):
        if labels[idx] in recover_labels:
            continue
        tmp_site1 = PeriodicSite("Ar", migrate_pos_frac[idx], stru.lattice)
        tmp_site2 = PeriodicSite("Ar", points[min_ids[idx]], stru.lattice)
        pts_len = [len(vorosites[0]), len(vorosites[1]), len(vorosites[2])]
        pt_tag = get_point_tag(min_ids[idx], pts_len)
        migrate_mindis[str(labels[idx])] = (pt_tag, tmp_site1.distance(tmp_site2))
        if tmp_site1.distance(tmp_site2) <= 0.5:
            recover_state[str(labels[idx])] = pt_tag
            recover_labels.append(labels[idx])
        else:
            recover_state[str(labels[idx])] = None

    recover_rate = len(recover_labels) / len(np.unique(labels))
    return recover_rate, recover_state, migrate_mindis

# 获取特定结构中
# 所有离子的有效半径，
# 迁移离子到最邻近离子表面距离与迁移离子半径的比值alpha，
# 迁移离子半径
def LocalEnvirCom(stru, migrant):
    # stru = Structure.from_file(filename)
    # val_eval = ValenceIonicRadiusEvaluator(stru)
    # radii = val_eval.radii
    coordination_list, radii = get_local_envir_fromstru(stru)
    
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
    
    coord_tmp = []
    nei_dis_tmp = []
    min_nei_dis_tmp = []
    migrant_paras = []
    migrant_radii = []
    for i in coordination_list:
        if migrant in i["label"]:
            #获取最邻近配位
            nearest_atom = i["coord_nei"][0]
            #最邻近配位原子的label
            nei_label = nearest_atom[0]._atom_site_label
            #到最邻近配位原子中心的距离
            nei_dis = nearest_atom[1]
            #最邻近配位原子的半径
            nei_radius = radii[nei_label]
            alpha_tmp = (nei_dis - nei_radius)/radii[i["label"]]
            
            #获取所有label迁移离子配位数
            coord_tmp.append(i["coord_num"])
            #获取不同label迁移离子到最邻近配位原子中心的距离与表面距离
            nei_dis_tmp.append(nei_dis)
            #获取所有label迁移离子到最邻近配位原子表面的距离
            min_nei_dis_tmp.append(nei_dis - nei_radius)
            migrant_paras.append(alpha_tmp)
            migrant_radii.append(radii[i["label"]])
            
    nei_dises = list(zip(coord_tmp, zip(nei_dis_tmp, min_nei_dis_tmp)))
    migrant_alpha = float(sum(migrant_paras))/len(migrant_paras)
    if migrant_alpha > 1.0:
        migrant_alpha = 1.0
    migrant_radius = float(sum(migrant_radii))/len(migrant_radii)
    return radii,migrant_radius,migrant_alpha,nei_dises,coordination_list

# 获取特定结构中
# 所有离子的有效半径
def getIonicRadii(filename):
    coordination_list, radii = get_local_envir(filename)
    print(radii)
    return radii

def AllCom8(filename, standard, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]

    # print(stru)
    #获取空间群号与符号
    symm_number,symm_sybol = parser.get_symme()
    #获取icsd cif文件中的对称操作
    sitesym = parser.get_sym_opt()
    radii,migrant_radius,migrant_alpha,nei_dises,coordination_list = LocalEnvirCom(stru,migrant)
    atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, True, None)
    
    prefixname = filename.replace(".cif","")
    vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(True)
    writeVaspFile(prefixname+"_origin_nofcs.vasp",atmnet,vornet)
    spg_vornet,uq_voids = get_equivalent_vornet(vornet, 0.01)

    symprec = 0.01
    sym_opt_num = len(sitesym)
    voids_num = len(vornet.nodes)

    writeNETFile(prefixname+"_origin_nofcs.net",atmnet,vornet)
    add_fcs_vornet = vornet.add_facecenters(faces)
    writeNETFile(prefixname+"_origin_addfcs.net",atmnet,add_fcs_vornet)
    writeVaspFile(prefixname+"_origin_addfcs.vasp",atmnet,add_fcs_vornet)

    spg_vornet,uq_voids = get_equivalent_vornet(add_fcs_vornet, 0.01)

    sym_vornet,voids =  get_labeled_vornet(add_fcs_vornet, sitesym, symprec)
    uni_voids_num = len(voids)

    voids_abs = []
    for void in sym_vornet.nodes:
        voids_abs.append(void[2])
    # print("voids")
    # print(voids_abs)

    bottlenecks = []
    for bt in sym_vornet.edges:
        bottlenecks.append(bt[2])
    # print("bottlenecks")
    # print(bottlenecks)
    
    # print("fcs",fcs)
    # facecenters = []
    # for fc in fcs:
    #     facecenters.append(atmnet.absolute_to_relative(fc[0], fc[1], fc[2]))
    # print("facecenters")
    # print(facecenters)

    vorosites = [voids_abs, bottlenecks, fcs]
    recover_rate, recover_state, migrate_mindis = rediscovery_kdTree(migrant,vorosites,stru)
    
    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet)

    conn_val = connection_values_list(prefixname+".resex", sym_vornet)
    minRad = standard*migrant_alpha*0.85
    dim_network,connect = ConnStatus(minRad, conn_val)
    writeVaspFile(prefixname+"_"+str(round(minRad,4))+".vasp",atmnet,sym_vornet,minRad,5.0)
    channels = Channel.findChannels(sym_vornet,atmnet,0,prefixname+"_0.net")
    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+"_"+str(round(minRad,4))+".net")
    
    return symm_sybol,symm_number,symprec,voids_num,sym_opt_num,uni_voids_num,recover_rate,recover_state,migrate_mindis

def AllCom7(filename, standard, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]  
    #获取空间群号与符号
    symm_number,symm_sybol = parser.get_symme()
    #获取icsd cif文件中的对称操作
    sitesym = parser.get_sym_opt()

    if rad_flag and effective_rad:
        #考虑如何利用migrant_radius与migrant_alpha
        radii,migrant_radius,migrant_alpha,nei_dises,coordination_list = LocalEnvirCom(stru,migrant)
    if migrant:
        atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, rad_flag, rad_file)
    else:
        atmnet = AtomNetwork.read_from_CIF(filename, radii, rad_flag, rad_file)
    
    prefixname = filename.replace(".cif","")
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(True)

    symprec = 0.01
    sym_opt_num = len(sitesym)
    voids_num = len(vornet.nodes)

    sym_vornet,voids =  get_labeled_vornet(vornet, sitesym, symprec)
    uni_voids_num = len(voids)

    # 衡量计算得到的对称性独立的间隙数与 对称性独立间隙最小数 之间的差异
    dif = abs(uni_voids_num - voids_num/sym_opt_num) / (voids_num/sym_opt_num)

    # positions = []
    # for i in vornet.nodes:
    #     positions.append(i[2])

    # print("Extend the scale to voids, bottlenecks, and face centers!")
    # print()
    bottlenecks = []
    for bt in sym_vornet.edges:
        bottlenecks.append(bt[2])
    facecenters = []
    for fidx in range(len(fcs)):
        facecenters.append(atmnet.absolute_to_relative(fcs[fidx][0], fcs[fidx][1], fcs[fidx][2]))
    
    vorosites = [voids, bottlenecks, facecenters]
    recover_rate, recover_state, true_recover_dis = rediscovery(migrant,vorosites,stru)
        
    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet)

    conn_val = connection_values_list(prefixname+".resex", sym_vornet)
    minRad = standard*migrant_alpha*0.85
    dim_network,connect = ConnStatus(minRad, conn_val)
    writeVaspFile(prefixname+"_"+str(round(minRad,4))+".vasp",atmnet,sym_vornet,minRad,5.0)

    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+"_"+str(round(minRad,4))+".net")
    dims_channel = []
    if len(channels)==0:
        dims_channel.append(0)
    else:
        for i in channels:
            dims_channel.append(i["dim"])
    return symm_sybol,symm_number,symprec,voids_num,sym_opt_num,uni_voids_num,recover_rate,recover_state,true_recover_dis


# 使用带半径的公式进行计算
def AllCom6(filename, standard, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}    
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]  
    #获取空间群号与符号
    symm_number,symm_sybol = parser.get_symme()
    #获取icsd cif文件中的对称操作
    sitesym = parser.get_sym_opt()

    if rad_flag and effective_rad:
        #考虑如何利用migrant_radius与migrant_alpha
        radii,migrant_radius,migrant_alpha,nei_dises,coordination_list = LocalEnvirCom(stru,migrant)
    if migrant:
        atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, rad_flag, rad_file)
    else:
        atmnet = AtomNetwork.read_from_CIF(filename, radii, rad_flag, rad_file)
    
    prefixname = filename.replace(".cif","")
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)

    symprec = 0.01
    sym_opt_num = len(sitesym)
    voids_num = len(vornet.nodes)

    sym_vornet,voids =  get_labeled_vornet(vornet, sitesym, symprec)
    uni_voids_num = len(voids)

    # 衡量计算得到的对称性独立的间隙数与 对称性独立间隙最小数 之间的差异
    dif = abs(uni_voids_num - voids_num/sym_opt_num) / (voids_num/sym_opt_num)

    recover_rate, recover_state, true_recover_dis = rediscovery(migrant,voids,stru)

    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet)

    conn_val = connection_values_list(prefixname+".resex", sym_vornet)
    minRad = standard*migrant_alpha*0.85
    dim_network,connect = ConnStatus(minRad, conn_val)
    writeVaspFile(prefixname+"_"+str(round(minRad,4))+".vasp",atmnet,sym_vornet,minRad,5.0)

    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+"_"+str(round(minRad,4))+".net")
    dims_channel = []
    if len(channels)==0:
        dims_channel.append(0)
    else:
        for i in channels:
            dims_channel.append(i["dim"])
    # return symm_sybol,symm_number,symprec,uni_voids_num,conn_val,connect,dim_network,dims_channel,migrant_alpha,radii,minRad,nei_dises,recover_rate, recover_state, true_recover_dis,coordination_list
    return symm_sybol, symm_number, sym_opt_num, voids_num, symprec, uni_voids_num, dis, recover_rate, recover_state, true_recover_dis
    
# 使用带半径的公式进行计算
def AllCom5(filename, standard, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    
    #获取空间群号与符号
    symm_number,symm_sybol = parser.get_symme()

    if rad_flag and effective_rad:
        #考虑如何利用migrant_radius与migrant_alpha
        radii,migrant_radius,migrant_alpha,nei_dises,coordination_list = LocalEnvirCom(stru,migrant)
    if migrant:
        atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, rad_flag, rad_file)
    else:
        atmnet = AtomNetwork.read_from_CIF(filename, radii, rad_flag, rad_file)
    
    prefixname = filename.replace(".cif","")
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)

    print("\nSymmetry number from cif: ", symm_number)
    max_symm = 0
    positions = []
    lattice = vornet.lattice
    for i in vornet.nodes:
        positions.append(i[2])
    for j in range(10):
        symprec = 0.01 + j*0.01
        symm_num_vornet =  get_symnum_sites(lattice, positions, symprec)

        if max_symm < symm_num_vornet:
            max_symm = symm_num_vornet 
            max_symm_info = (max_symm, symprec)

        if symm_num_vornet == symm_number:
            print("Distance tolerance in Cartesian coordinates to find crystal symmetry: ",symprec)
            print("Symmetry number from Voronoi network: ", symm_num_vornet)
            print("\n")
            sym_vornet, voids = get_equivalent_vornet(vornet,symprec)
            break

        #在0.01-0.10范围内均无法得到与cif文件中一致的空间群号，使用在此过程中出现的最大值代替
        elif j == 9:
            print("The Symmetry calculated from Vornet (with symprec 0.01-0.1) is different from that obtained from cif files.")
            print("Using the lagest value of symm_num_vornet instead!")
            symm_num_vornet = max_symm_info[0]
            symprec = max_symm_info[1]
            print("Distance tolerance in Cartesian coordinates to find crystal symmetry: ",symprec)
            print("Symmetry number of Voronoi network: ", symm_num_vornet)
            print("\n")
            sym_vornet, voids = get_equivalent_vornet(vornet,symprec)

    recover_rate, recover_state, true_recover_dis = rediscovery(migrant,voids,stru)

    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet)

    conn_val = connection_values_list(prefixname+".resex", sym_vornet)
    minRad = standard*migrant_alpha*0.85
    dim_network,connect = ConnStatus(minRad, conn_val)
    writeVaspFile(prefixname+"_"+str(round(minRad,4))+".vasp",atmnet,sym_vornet,minRad,5.0)

    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+"_"+str(round(minRad,4))+".net")
    dims_channel = []
    if len(channels)==0:
        dims_channel.append(0)
    else:
        for i in channels:
            dims_channel.append(i["dim"])
    return symm_sybol,symm_number,symm_num_vornet,symprec,conn_val,connect,dim_network,dims_channel,migrant_alpha,radii,minRad,nei_dises,recover_rate, recover_state, true_recover_dis,coordination_list

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

    writeBIFile(prefixname+"_origin.bi",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet,rad_store_in_vasp)

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
    if migrant:
        os.remove(remove_filename)

    prefixname = filename.replace(".cif","")
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    sym_vornet,voids = get_Symmetry(atmnet, vornet)

    writeBIFile(prefixname+"_origin.bi",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet,rad_store_in_vasp)

    minRad = standard*migrant_alpha*0.85
    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+".net")
    
    dims = []
    if len(channels)==0:
        dims.append(0)
    else:
        for i in channels:
            dims.append(i["dim"])
    return migrant_alpha,radii,minRad,nei_dises,dims,voids

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

    writeBIFile(prefixname+"_origin.bi",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet,rad_store_in_vasp)

    writeVaspFile(prefixname+"_selected.vasp",atmnet,sym_vornet,rad_store_in_vasp,minRad,maxRad)

    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+"_channel.net")
    
    dims = []
    for i in channels:
        dims.append(i["dim"])
    return dims,voids

def AllCom(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=True, rad_file=None, minRad=0.0, maxRad=0.0):
    radii = {}
    if rad_flag and effective_rad:
        #考虑如何利用migrant_radius与migrant_alpha
        symm_number, radii,migrant_radius,migrant_alpha, nei_dises,coordination_list = LocalEnvirCom(filename,migrant)
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

    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet)

    probe_rad = migrant_radius*migrant_alpha
    minRad = migrant_radius*migrant_alpha*0.85
    maxRad = migrant_radius*migrant_alpha*1.15
    print(minRad)
    print(maxRad)

    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet,minRad,maxRad)
    writeVaspFile(prefixname+"_selected.vasp",atmnet,sym_vornet,minRad,maxRad)
    
    channels = Channel.findChannels(sym_vornet,atmnet,0.60,prefixname+".net")
    conn = connection_values_list(prefixname+".resex", sym_vornet)
    oneD,twoD,threeD = ConnStatus(minRad, conn)
    
    dims = []
    for i in channels:
        dims.append(i["dim"])
    return conn,oneD,twoD,threeD,nei_dises,dims,voids,coordination_list


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
    writeBIFile(prefixname+"_origin.bi",atmnet,vornet)
    writeBIFile(prefixname+"_selected.bi",atmnet,vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,vornet,rad_store_in_vasp)
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
    return [oneD,twoD,threeD]

#根据连通数值列表，判断某个结构的连通性。给定一个原子的半径，判断它是否是1D，2D，3D导通
def ConnStatus(radius,connlist):
    connects = []
    for i in connlist:
        if radius > i:
            connects.append(False)
        else:
            connects.append(True)
    dim_net = connects.count(True)
    return dim_net,connects
    
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
