import os
import sys
import re
from zeo.netstorage import AtomNetwork
from zeo.channel import Channel
from zeo.area_volume import asa_new
from zeo.netio import *
from zeo.ionic_radii import get_ionic_radii
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator



#获取输入结构中的离子半径
def EffectiveRadCom(filename):
    # stru = Structure.from_file(filename)
    # val_eval = ValenceIonicRadiusEvaluator(stru)
    # radii = val_eval.radii
    radii = get_ionic_radii(filename)
    return radii

#计算某个结构的瓶颈、间隙和连通性
def Computation_new(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=False, rad_file=None, rad_store_in_vasp=True, minRad=0.0, maxRad=0.0):
    Ri="None"
    Rf="None"
    Rif="None"
    sucess=False
    
    #use pymatgen compute ionic radii or find in the cif file
    try:
        radii = {}
        if rad_flag and effective_rad:
#            stru = Structure.from_file(filename)
#            val_eval = ValenceIonicRadiusEvaluator(stru)
#            radii = val_eval.radii
            radii = get_ionic_radii(filename)
            #为防止cif文件中不含价态信息，额外存入不含价态信息的半径
            radii_keys = list(radii.keys())
            for key in radii_keys:
                radii[re.sub('[^a-zA-Z]','',key)] = radii[key]
            print(radii)
        if migrant:
            remove_filename = getRemoveMigrantFilename(filename,migrant)
        else:
            remove_filename = filename
        atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
        sucess,vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
        #delete temp file
        os.remove(remove_filename)
        
        #如果进行Voronoi分解成功，则将结果保存到相应文件中。
        if sucess:
            prefixname = filename.replace(".cif","")
            writeBIFile(prefixname+"_orgin.bi",atmnet,vornet)
            writeVaspFile(prefixname+"_orgin.vasp",atmnet,vornet,rad_store_in_vasp)
            writeVaspFile(prefixname+"_selected.vasp",atmnet,vornet,rad_store_in_vasp,minRad,maxRad)
            Ri,Rf,Rif = atmnet.through_VorNet(prefixname+".res",0)
            atmnet.calculate_free_sphere_parameters(prefixname+".resex")
    except IOError:
        print("cant compute present file. error file is ",filename)

#计算通道和ASA
def com(filename,probe_rad,num_sample):
    radii = get_ionic_radii(filename)
    atmnet = AtomNetwork.read_from_CIF(filename, radii, rad_flag=False)
    sucess,vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
    writeZVisFile("test.zvis", False, atmnet, vornet)
    Channel.findChannelsInVornet(vornet,probe_rad,"test.zchan")
    asa_new("test.zsa",False,atmnet,probe_rad,probe_rad,1000)


#计算某个结构的瓶颈、间隙和连通性
def Computation(filename, migrant=None, rad_flag=True, effective_rad=False, rad_file=None, rad_store_in_vasp=True, minRad=0.0, maxRad=0.0):
    Ri="None"
    Rf="None"
    Rif="None"
    sucess=False
	
    #use pymatgen compute ionic radii or find in the cif file
    try:
        radii = {}
        if rad_flag and effective_rad:
#            stru = Structure.from_file(filename)
#            val_eval = ValenceIonicRadiusEvaluator(stru)
#            radii = val_eval.radii
            radii = get_ionic_radii(filename)
            #为防止cif文件中不含价态信息，额外存入不含价态信息的半径
            radii_keys = list(radii.keys())
            for key in radii_keys:
                radii[re.sub('[^a-zA-Z]','',key)] = radii[key]
            print(radii)
        if migrant:
            remove_filename = getRemoveMigrantFilename(filename,migrant)
        else:
            remove_filename = filename
        atmnet = AtomNetwork.read_from_CIF(remove_filename, radii, rad_flag, rad_file)
        sucess,vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
        #delete temp file
        os.remove(remove_filename)
		#如果进行Voronoi分解成功，则将结果保存到相应文件中。
        if sucess:
            prefixname = filename.replace(".cif","")
            writeBIFile(prefixname+"_orgin.bi",atmnet,vornet)
            writeVaspFile(prefixname+"_orgin.vasp",atmnet,vornet,rad_store_in_vasp)
            writeVaspFile(prefixname+"_selected.vasp",atmnet,vornet,rad_store_in_vasp,minRad,maxRad)
            Ri,Rf,Rif = atmnet.through_VorNet(prefixname+".res",0)
            atmnet.calculate_free_sphere_parameters(prefixname+".resex")
        return sucess,Ri,Rf,Rif
    except IOError:
        print("cant compute present file. error file is ",filename)
        return sucess,Ri,Rf,Rif		
"""    
	except IOError:
        print("cant write information to result file. error file is ",filename)
        return sucess,Ri,Rf,Rif
    except ValueError:
        print("cant read or move migrant from cif file. error file is ",filename)
"""
#批处理计算
#added at 20180626
def Computation_batch(path, migrant=None, effective_rad=False, minRad=0.0, maxRad=0.0):
    filenames = BatchReadFilename(path,".cif")
    #print(filenames)
    output_path = path+"results/"
    for filename in filenames:
        filename = path+filename
        sucess,Ri,Rf,Rif = Computation(filename, migrant, rad_flag=True, effective_rad=effective_rad, rad_file=None, rad_store_in_vasp=True, minRad=minRad, maxRad=maxRad)
        print(filename+" compute complete1!")
    print("batch compute complete!")

#批量读取指定文件夹下指定类型文件
def BatchReadFilename(path,filetype):
    filenames=[]
    for i in os.listdir(path):
        if filetype in i:
            #filenames.append(i.replace(filetype,''))
            filenames.append(i)
    return filenames

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
        print("3D connected!")
    if((!connection[1] and connection[2] and connection[3]) or (connection[1] and !connection[2] and connection[3]) or (connection[1] and connection[2] and !connection[3])):
        print("2D connected!")
    if ((connection[1] and !connection[2] and !connection[3]) or (!connection[1] and connection[2] and !connection[3]) or (!connection[1] and !connection[2] and !connection[3])):
        print("1D connected!")

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
    
"""
#批处理计算
#added at 20180604
def Computation_batch(path, migrant=None, minRad=0.0, maxRad=0.0):
    #创建计算结果保存目录
    if not os.path.exists(path+"results"):
        os.mkdir(path+"results")
        print("create results directory successful !")
    else:
        print(path+"results already exit!")
    filenames = BatchReadFilename(path,".cif")
    #print(filenames)
    output_path = path+"results/"
	
    #result_file = open(output_path+"computation_result.txt","w")
    #connection_file = open(output_path+"connection_result.txt","w")
    for filename in filenames:
        filename = path+filename
        sucess_flag, Ri, Rf, Rif = Computation(filename, migrant, rad_flag=True, rad_file=None, rad_store_in_vasp=True, minRad=minRad, maxRad=maxRad)
        if sucess_flag:
            filestr = "{}    {}    {}    {}".format(filename, Ri, Rf, Rif)
            connection_file.write(filestr)
            connection_file.write("\n")
        sucess_flag = '{}'.format(sucess_flag)
        result_file.write(filename+"   "+ sucess_flag +'\n')
        result_file.write("\n")
        print(filename+" compute complete1!")
    print("batch compute complete!")

#计算某个结构的瓶颈和间隙
def BIComputation(filename, migrant=None, rad_flag=True, rad_file=None, rad_store_in_vasp=True, minRad=0.0):
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        renove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, rad_flag, rad_file)
    vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition()
    prefixname = filename.replace(".cif","")
    writeBIFile(prefixname+".bi",atmnet,vornet,minRad)
    writeVaspFile(prefixname+".vasp",atmnet,vornet,rad_store_in_vasp,minRad)
    
#计算某个结构的连通性
def Connection(filename, probe_size, migrant=None, rad_flag=True, rad_file=None):
    if migrant:
        remove_filename = getRemoveMigrantFilename(filename,migrant)
    else:
        remove_filename = filename
    atmnet = AtomNetwork.read_from_CIF(remove_filename, rad_flag, rad_file)
    prefixname = filename.replace(".cif","")
    return atmnet.through_VorNet(prefixname+".res",probe_size)
	
"""
