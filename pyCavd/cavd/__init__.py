import os
import sys
import re
from cavd.netstorage import AtomNetwork
from cavd.netstorage import connection_values_list
from cavd.channel import Channel
from cavd.area_volume import asa_new
from cavd.netio import *
from cavd.ionic_radii import get_ionic_radii
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator

#获取特定结构中离子的有效半径
def EffectiveRadCom(filename):
    # stru = Structure.from_file(filename)
    # val_eval = ValenceIonicRadiusEvaluator(stru)
    # radii = val_eval.radii
    radii = get_ionic_radii(filename)
    radii_keys = list(radii.keys())
    
    #为了防止保存的半径信息无法匹配此处对半径信息做特殊处理，如Ag+的半径会保存为Ag、Ag+、Ag1+1
    for key in radii_keys:
        radii[re.sub('[^a-zA-Z]','',key)] = radii[key]
        if re.search('[A-Z][a-z]*\+|[A-Z][a-z]*\-', key) != None:
            s1 = re.sub('\+','1+',key)
            radii[re.sub('\-','1-',s1)] = radii[key]
    print(radii)
    return radii

def AllCom(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True, minRad=0.0, maxRad=0.0):
    radii = {}
    if rad_flag and effective_rad:
        radii = EffectiveRadCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    if migrant:
        os.remove(remove_filename)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    
    prefixname = filename.replace(".cif","")
    writeBIFile(prefixname+"_orgin.bi",atmnet,vornet)
    writeVaspFile(prefixname+"_orgin.vasp",atmnet,vornet,rad_store_in_vasp)
    writeVaspFile(prefixname+"_selected.vasp",atmnet,vornet,rad_store_in_vasp,minRad,maxRad)
    conn = connection_values_list(prefixname+".resex", vornet)
    oneD,twoD,threeD = ConnStatus(probe_rad, conn)
    Channel.findChannels(vornet,probe_rad,prefixname+".net")
    asa_new(prefixname+".zsa",False,atmnet,probe_rad,probe_rad,num_sample)
    writeZVisFile(prefixname+".zvis", rad_flag, atmnet, vornet)
    return conn,oneD,twoD,threeD

#计算指定结构的瓶颈和间隙
def BIComputation(filename, migrant=None, rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True,  minRad=0.0, maxRad=0.0):
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
        radii = EffectiveRadCom(filename)
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
        radii = EffectiveRadCom(filename)
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
        radii = EffectiveRadCom(filename)
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    if migrant:
        os.remove(remove_filename)
    prefixname = filename.replace(".cif","")
    #Channel.findChannelsInVornet(vornet,probe_rad,prefixname+".zchan")
    Channel.findChannels(vornet,probe_rad,prefixname+".net")

#计算ASA
def ASACom(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    if rad_flag and effective_rad:
        radii = EffectiveRadCom(filename)
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
        radii = EffectiveRadCom(filename)
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
