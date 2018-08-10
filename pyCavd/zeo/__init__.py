import os
import sys
import re
from zeo.netstorage import AtomNetwork
from zeo.channel import Channel
from zeo.area_volume import asa_new
from zeo.netio import *
from zeo.ionic_radii import get_ionic_radii
from pymatgen.core.structure import Structure
#from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator

#获取输入结构中的离子半径
def EffectiveRadCom(filename):
    # stru = Structure.from_file(filename)
    # val_eval = ValenceIonicRadiusEvaluator(stru)
    # radii = val_eval.radii
    radii = get_ionic_radii(filename)
	radii_keys = list(radii.keys())
    
    #为了防止保存的半径信息无法匹配此处对半径信息做特殊处理，如Ag+的半径会保存为Ag、Ag+、Ag1+ 1
    for key in radii_keys:
		radii[re.sub('[^a-zA-Z]','',key)] = radii[key]
        if re.search('[A-Z][a-z]*\+|[A-Z][a-z]*\-', key) != None:
            s1 = re.sub('\+','1+',key)
            radii[re.sub('\-','1-',s1)] = radii[key]
    print(radii)
    return radii

def AllCom(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=False, rad_file=None, rad_store_in_vasp=True, minRad=0.0, maxRad=0.0):
    radii = {}
    if rad_flag and effective_rad:
        radii = EffectiveRadCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
	vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    os.remove(remove_filename)
	writeBIFile(prefixname+"_orgin.bi",atmnet,vornet)
	writeVaspFile(prefixname+"_orgin.vasp",atmnet,vornet,rad_store_in_vasp)
	writeVaspFile(prefixname+"_selected.vasp",atmnet,vornet,rad_store_in_vasp,minRad,maxRad)
	writeZVisFile(prefixname+".zvis", rad_flag, atmnet, vornet)
    Channel.findChannelsInVornet(vornet,probe_rad,prefixname+".zchan")
    asa_new(prefixname+".zsa",False,atmnet,probe_rad,probe_rad,num_sample)
	
	conn = connection_values_list(prefixname+".resex", vornet)
    oneD,twoD,threeD = ConnValue(probe_rad, conn)
    return conn,oneD,twoD,threeD

#计算某个结构的瓶颈和间隙
def BIComputation(filename, migrant=None, rad_flag=True, effective_rad=False, rad_file=None, rad_store_in_vasp=True, minRad=0.0):
    radii = {}
    if rad_flag and effective_rad:
        radii = EffectiveRadCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
        atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
        vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
        #delete temp file
        os.remove(remove_filename)
        prefixname = filename.replace(".cif","")
        writeBIFile(prefixname+"_orgin.bi",atmnet,vornet)
        writeVaspFile(prefixname+"_orgin.vasp",atmnet,vornet,rad_store_in_vasp)
        writeVaspFile(prefixname+"_selected.vasp",atmnet,vornet,rad_store_in_vasp,minRad,maxRad)

#计算某个结构最大自由球体半径，最大包含球体半径和沿着最大自由球体路径上的最大包含球体半径：Rf Ri Rif
def ConnValCom(filename, migrant=None, rad_flag=True, effective_rad=False, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = EffectiveRadCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    prefixname = filename.replace(".cif","")
    Ri,Rf,Rif = atmnet.through_VorNet(prefixname+".res",0)
    return Ri,Rf,Rif
    
#计算某个结构的连通性状态列表，存放1D，2D，3D连通信息元素，这些元素为一个字典，字典的键为Rf、Ri、Rif，值为对应的数值
def ConnValListCom(filename, migrant=None, rad_flag=True, effective_rad=False, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = EffectiveRadCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    prefixname = filename.replace(".cif","")
    #连通性状态列表，存放1D，2D，3D连通信息元素，这些元素为一个字典，字典的键为Rf、Ri、Rif，值为对应的数值
    #需重写该函数，需返回该列表
    conn = atmnet.calculate_free_sphere_parameters(prefixname+".resex")
    return conn

#判断某个结构的连通性,给定一个原子的半径，判断它是否是1D，2D，3D导通
def ConnValueCom(filename, radius, migrant=None, rad_flag=True, effective_rad=False, rad_file=None):
    connlist = ConnValListCom(filename, migrant, rad_flag, effective_rad, rad_file)
    connection = []
    for i in connlist:
        if radius <= connlist[i]["Rf"]:
            connection[i] = True
    if(connection[1] and connection[2] and connection[3])
        threeD = true
    if((!connection[1] and connection[2] and connection[3]) or (connection[1] and !connection[2] and connection[3]) or (connection[1] and connection[2] and !connection[3])):
        twoD = True
    if ((connection[1] and !connection[2] and !connection[3]) or (!connection[1] and connection[2] and !connection[3]) or (!connection[1] and !connection[2] and !connection[3])):
        oneD = True
    return oneD,twoD,threeD

#根据连通性列表，判断某个结构的连通性,给定一个原子的半径，判断它是否是1D，2D，3D导通
def ConnValue(radius,connlist):
    connection = []
    for i in connlist:
        if radius <= connlist[i]["Rf"]:
            connection[i] = True
    if(connection[1] and connection[2] and connection[3])
        threeD = true
    if((!connection[1] and connection[2] and connection[3]) or (connection[1] and !connection[2] and connection[3]) or (connection[1] and connection[2] and !connection[3])):
        twoD = True
    if ((connection[1] and !connection[2] and !connection[3]) or (!connection[1] and connection[2] and !connection[3]) or (!connection[1] and !connection[2] and !connection[3])):
        oneD = True
    return oneD,twoD,threeD
    
#计算通道
def ChannelCom(filename, probe_rad, migrant=None, rad_flag=True, effective_rad=False, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = EffectiveRadCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    Channel.findChannelsInVornet(vornet,probe_rad,"test.zchan")

#计算ASA
def ASACom(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=False, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = EffectiveRadCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    asa_new("test.zsa",False,atmnet,probe_rad,probe_rad,num_sample)

#计算空隙网络
def VoidNetCom(filename, migrant=None, rad_flag=True, effective_rad=False, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = EffectiveRadCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    writeZVisFile("test.zvis", False, atmnet, vornet)
