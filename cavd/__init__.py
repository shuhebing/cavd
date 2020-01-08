# -*- coding: UTF-8 -*-
import os
import re
from cavd.netio import *
from monty.io import zopen
from cavd.channel import Channel
from cavd.area_volume import asa_new
from cavd.cavd_consts import LOWER_THRESHOLD, UPPER_THRESHOLD
from cavd.netstorage import AtomNetwork, connection_values_list
from cavd.recovery import rediscovery, rediscovery_kdTree, rediscovery_byRad_kdTree
from cavd.get_Symmetry import get_symnum_sites, get_equivalent_vornet,get_labeled_vornet
from cavd.local_environment import CifParser_new, LocalEnvirCom, get_local_envir_fromstru
from cavd.channel_analysis import MigrationPaths

def outChannelToPOSCAR(filename, migrant, ntol=0.02, rad_flag=True, lower=0.0, upper=10.0, rad_dict=None):
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    print(stru)
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
    
    add_fcs_vornet = vornet.add_facecenters(faces)
    sym_vornet,voids =  get_labeled_vornet(add_fcs_vornet, sitesym)
    prefixname = filename.replace(".cif","")
    channels = Channel.findChannels2(sym_vornet, atmnet, lower, upper, prefixname+".net")
    
    migratPath = MigrationPaths(stru, migrant, channels)
    
    allPaths = migratPath.comAllPaths()
    return allPaths
    
def outVesta(filename, migrant, ntol=0.02, rad_flag=True, lower=None, upper=10.0, rad_dict=None):
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
    add_fcs_vornet = vornet.add_facecenters(faces)
    
    sym_vornet,voids =  get_labeled_vornet(add_fcs_vornet, sitesym, ntol)
    
    prefixname = filename.replace(".cif","")
   
    # calculate the connection values
    conn_val = connection_values_list(prefixname+".resex", sym_vornet)

    channels = Channel.findChannels2(sym_vornet, atmnet, lower, upper, prefixname+".net")
    dims = []
    for i in channels:
        dims.append(i["dim"])
    # output vesta file for visiualization
    Channel.writeToVESTA(channels, atmnet, prefixname)

    return dims, conn_val

# This function need to be updated.
def AllCom10(filename, minRad, maxRad, ntol=0.02, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]

    symm_number,symm_sybol = parser.get_symme()
    sitesym = parser.get_sym_opt()
    radii,migrant_radius,migrant_alpha,nei_dises,coordination_list = LocalEnvirCom(stru,migrant)
    atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, True, None)
    
    prefixname = filename.replace(".cif","")
    vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(True, ntol)
    # writeVaspFile(prefixname+"_origin_nofcs.vasp",atmnet,vornet)
    # spg_vornet,uq_voids = get_equivalent_vornet(vornet, 0.01)
    print("len fcs: ", len(fcs))
    symprec = 0.01
    sym_opt_num = len(sitesym)
    voids_num = len(vornet.nodes)

    # writeNETFile(prefixname+"_origin_nofcs.net",atmnet,vornet)
    add_fcs_vornet = vornet.add_facecenters(faces)
    # writeNETFile(prefixname+"_origin_addfcs.net",atmnet,add_fcs_vornet)
    # writeVaspFile(prefixname+"_origin_addfcs.vasp",atmnet,add_fcs_vornet)

    # spg_vornet,uq_voids = get_equivalent_vornet(add_fcs_vornet, 0.01)

    sym_vornet,voids =  get_labeled_vornet(add_fcs_vornet, sitesym, symprec)
    uni_voids_num = len(voids)

    voids_abs = []
    voids_rad = []
    for void in sym_vornet.nodes:
        voids_abs.append(void[2])
        voids_rad.append(void[3])
    # print("voids")
    # print(voids_abs)

    bottlenecks = []
    bottlenecs_rad = []
    for bt in sym_vornet.edges:
        frac_bt = bt[2]
        # frac_bt = [round(p%1.0, 6) for p in frac_bt]
        # frac_bt = [p%1.0 for p in frac_bt]
        # if frac_bt not in bottlenecks:
        bottlenecks.append(frac_bt)
        bottlenecs_rad.append(bt[3])
    # print("bottlenecks")
    # print(bottlenecks)
    
    # print("fcs",fcs)
    fcens = []
    fcens_rad = []
    for fc in fcs:
        fcens.append(fc[0])
        fcens_rad.append(fc[1])
    # print("facecenters")
    # print(facecenters)

    vorosites = [voids_abs, bottlenecks, fcens]
    vororad = [voids_rad, bottlenecs_rad, fcens_rad]
    # recover_rate, recover_state, migrate_mindis = rediscovery_kdTree(migrant,vorosites,stru)
    recover_rate, recover_state, migrate_mindis = rediscovery_byRad_kdTree(stru, migrant, vorosites, vororad)
    
    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet)

    conn_val = connection_values_list(prefixname+".resex", sym_vornet)
    dim_network,connect = ConnStatus(conn_val, minRad, maxRad)
    
    channels = Channel.findChannels2(sym_vornet,atmnet,minRad,maxRad,prefixname+"_select.net")
    
    # output vesta file for visiualization
    Channel.writeToVESTA(channels, atmnet, prefixname)

    dims_channel = []
    if len(channels)==0:
        dims_channel.append(0)
    else:
        for i in channels:
            dims_channel.append(i["dim"])
    
    writeVaspFile(prefixname+"_select.vasp",atmnet,sym_vornet,minRad,maxRad)
    # channels = Channel.findChannels(sym_vornet,atmnet,0,prefixname+"_0.net")
            
    return radii,symm_sybol,symm_number,symprec,voids_num,sym_opt_num,uni_voids_num,minRad,migrant_alpha,nei_dises,migrant_radius,conn_val,connect,dim_network,dims_channel,recover_rate,recover_state,migrate_mindis,coordination_list

# This function need to be updated.
def AllCom5(filename, standard, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    
    symm_number,symm_sybol = parser.get_symme()

    if rad_flag and effective_rad:
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

        # If the space group number consistent with the cif file cannot be obtained within the range of 0.01-0.10, 
        # use the largest space group number that appears in the process instead.
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
    dim_network,connect = ConnStatus(conn_val,minRad)
    writeVaspFile(prefixname+"_"+str(round(minRad,4))+".vasp",atmnet,sym_vornet,minRad,5.0)

    channels = Channel.findChannels(sym_vornet,atmnet,minRad,prefixname+"_"+str(round(minRad,4))+".net")

    # output vesta file for visiualization
    Channel.writeToVESTA(channels, atmnet, prefixname)

    dims_channel = []
    if len(channels)==0:
        dims_channel.append(0)
    else:
        for i in channels:
            dims_channel.append(i["dim"])
    return symm_sybol,symm_number,symm_num_vornet,symprec,conn_val,connect,dim_network,dims_channel,migrant_alpha,radii,minRad,nei_dises,recover_rate, recover_state, true_recover_dis,coordination_list

# This function need to be updated.
def AllCom(filename, minRad, maxRad, migrant=None, rad_flag=True, effective_rad=True, rad_file=None):
    radii = {}
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]

    symm_number,symm_sybol = parser.get_symme()
    sitesym = parser.get_sym_opt()
    radii,migrant_radius,migrant_alpha,nei_dises,coordination_list = LocalEnvirCom(stru,migrant)
    atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, True, None)
    
    prefixname = filename.replace(".cif","")
    vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(True)
   
    symprec = 0.01
    sym_opt_num = len(sitesym)
    voids_num = len(vornet.nodes)

    # writeNETFile(prefixname+"_origin_nofcs.net",atmnet,vornet)
    add_fcs_vornet = vornet.add_facecenters(faces)
    # writeNETFile(prefixname+"_origin_addfcs.net",atmnet,add_fcs_vornet)
    # writeVaspFile(prefixname+"_origin_addfcs.vasp",atmnet,add_fcs_vornet)

    # spg_vornet,uq_voids = get_equivalent_vornet(add_fcs_vornet, 0.01)

    sym_vornet,voids =  get_labeled_vornet(add_fcs_vornet, sitesym, symprec)
    uni_voids_num = len(voids)

    voids_abs = []
    voids_rad = []
    for void in sym_vornet.nodes:
        voids_abs.append(void[2])
        voids_rad.append(void[3])
    # print("voids")
    # print(voids_abs)

    bottlenecks = []
    bottlenecs_rad = []
    for bt in sym_vornet.edges:
        frac_bt = bt[2]
        # frac_bt = [round(p%1.0, 6) for p in frac_bt]
        # frac_bt = [p%1.0 for p in frac_bt]
        # if frac_bt not in bottlenecks:
        bottlenecks.append(frac_bt)
        bottlenecs_rad.append(bt[3])
    # print("bottlenecks")
    # print(bottlenecks)
    
    # print("fcs",fcs)
    fcens = []
    fcens_rad = []
    for fc in fcs:
        fcens.append(fc[0])
        fcens_rad.append(fc[1])
    # print("facecenters")
    # print(facecenters)

    vorosites = [voids_abs, bottlenecks, fcens]
    vororad = [voids_rad, bottlenecs_rad, fcens_rad]
    # recover_rate, recover_state, migrate_mindis = rediscovery_kdTree(migrant,vorosites,stru)
    recover_rate, recover_state, migrate_mindis = rediscovery_byRad_kdTree(stru, migrant, vorosites, vororad)
    
    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet)

    conn_val = connection_values_list(prefixname+".resex", sym_vornet)
    dim_network,connect = ConnStatus(conn_val, minRad, maxRad)
    
    channels = Channel.findChannels2(sym_vornet,atmnet,minRad,maxRad,prefixname+"_select.net")

    # output vesta file for visiualization
    Channel.writeToVESTA(channels, atmnet, prefixname)

    dims_channel = []
    if len(channels)==0:
        dims_channel.append(0)
    else:
        for i in channels:
            dims_channel.append(i["dim"])
    
    
    writeVaspFile(prefixname+"_select.vasp",atmnet,sym_vornet,minRad,maxRad)
    # channels = Channel.findChannels(sym_vornet,atmnet,0,prefixname+"_0.net")
            
    return radii,symm_sybol,symm_number,symprec,voids_num,sym_opt_num,uni_voids_num,minRad,migrant_alpha,nei_dises,migrant_radius,conn_val,connect,dim_network,dims_channel,recover_rate,recover_state,migrate_mindis,coordination_list

def bmd_com(filename, migrant, rad_flag=True, lower=None, upper=10.0, rad_dict=None, symprec=0.01):
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
    
    prefixname = filename.replace(".cif","")
    #for cst paper
    # vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(True)
    vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(False)
    sym_vornet,voids =  get_labeled_vornet(vornet, sitesym, symprec)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,sym_vornet)
    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet)
    
    # add_fcs_vornet = vornet.add_facecenters(faces)
    # sym_vornet,voids =  get_labeled_vornet(add_fcs_vornet, sitesym, symprec)
    
    # voids_abs = []
    # for void in sym_vornet.nodes:
        # voids_abs.append(void[2])

    # bottlenecks = []
    # for bt in sym_vornet.edges:
        # bottlenecks.append(bt[2])
    
    # fcens = []
    # for fc in fcs:
        # fcens.append(fc[0])

    # vorosites = [voids_abs, bottlenecks, fcens]
    # recover_rate, recover_state, migrate_mindis = rediscovery_kdTree(stru, migrant, vorosites)
    
    # for cst paper
    migrate_mindis = None

    minRad = 0.0
    if lower:
        minRad = lower
    else:
        standard = LOWER_THRESHOLD[migrant]
        minRad = standard*migrant_alpha*0.85
    conn_val = connection_values_list(prefixname+".resex", sym_vornet)
    dim_network,connect = ConnStatus(conn_val,minRad)
    writeVaspFile(prefixname+".vasp",atmnet, sym_vornet, minRad, upper)
    channels = Channel.findChannels(sym_vornet, atmnet, minRad, prefixname+".net")

    # output vesta file for visiualization
    Channel.writeToVESTA(channels, atmnet, prefixname)
    
    dims = []
    for i in channels:
        dims.append(i["dim"])

    return radii, minRad, conn_val, connect, dim_network, dims, migrate_mindis

# Calculate interstice and bottleneck for given structure.
def BIComputation(filename, migrant, rad_flag=True, lower=0.0, upper=0.0, rad_dict=None):
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]

    species = [str(sp).replace("Specie ","") for sp in stru.species]
    elements = [re.sub('[^a-zA-Z]','',sp) for sp in species]
    if migrant not in elements:
        raise ValueError("The input migrant ion not in the input structure! Please check it.")
    coord_list, effec_radii = get_local_envir_fromstru(stru)

    radii = {}
    if rad_flag:
        if rad_dict:
            radii = rad_dict
        else:
            radii = effec_radii
    
    atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, rad_flag)
    
    prefixname = filename.replace(".cif","")
    vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(True)
    add_fcs_vornet = vornet.add_facecenters(faces)
    writeVaspFile(prefixname+"_origin.vasp",atmnet,add_fcs_vornet)
    
    if lower and not upper:
        writeVaspFile(prefixname+"_selected.vasp",atmnet, add_fcs_vornet, lower, 10.0)
    if not lower and upper:
        writeVaspFile(prefixname+"_selected.vasp",atmnet, add_fcs_vornet, 0.0, upper)
    if lower and upper:
        writeVaspFile(prefixname+"_selected.vasp",atmnet, add_fcs_vornet, lower, upper)
    
# Calculate a list of connected values (Rf alongs the a, b, and c direction) for a structure.
def ConnValListCom(filename, migrant=None, rad_flag=True, rad_dict=None):
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    
    species = [str(sp).replace("Specie ","") for sp in stru.species]
    elements = [re.sub('[^a-zA-Z]','',sp) for sp in species]
    if migrant not in elements:
        raise ValueError("The input migrant ion not in the input structure! Please check it.")
    coord_list,effec_radii = get_local_envir_fromstru(stru)
    
    radii = {}
    if rad_flag:
        if rad_dict:
            radii = rad_dict
        else:
            radii = effec_radii
    
    atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, rad_flag)

    vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(True)
    add_fcs_vornet = vornet.add_facecenters(faces)

    prefixname = filename.replace(".cif","")
    conn = connection_values_list(prefixname+".resex",add_fcs_vornet)
    return conn

# Determine the connectivity of a structure.
# According to the given radius of the target ion to determine the dimension of the interstitial network.
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

# Determine the connectivity of a structure based on the list of connected values.
def ConnStatus(connlist, minRad, maxRad=10.0):
    connects = []
    for i in connlist:
        if minRad <= i and maxRad >= i:
            connects.append(True)
        else:
            connects.append(False)
    dim_net = connects.count(True)
    return dim_net,connects
    
# Compute the channel
def ChannelCom(filename, probe_rad = None, migrant=None, rad_flag=True, rad_dict=None, symprec=0.01):
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
       
    radii = {}
    if rad_flag:
        if rad_dict:
            radii = rad_dict
        else:
            coordination_list, effec_radii = get_local_envir_fromstru(stru)
            radii = effec_radii
    
    species = [str(sp).replace("Specie ","") for sp in stru.species]
    elements = [re.sub('[^a-zA-Z]','',sp) for sp in species]
    if migrant:
        if migrant not in elements:
            raise ValueError("The input migrant ion not in the input structure! Please check it.")
        atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, rad_flag)
    else:
        atmnet = AtomNetwork.read_from_CIF(filename, radii, rad_flag)
        
    vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(True)
    add_fcs_vornet = vornet.add_facecenters(faces)
    
    sitesym = parser.get_sym_opt()
    sym_vornet,voids =  get_labeled_vornet(add_fcs_vornet, sitesym, symprec)
    prefixname = filename.replace(".cif","")
    writeNETFile(prefixname+"_origin.net",atmnet,sym_vornet)

    channels = Channel.findChannels(sym_vornet, atmnet, probe_rad, prefixname+".net")
    
    # output vesta file for visiualization
    Channel.writeToVESTA(channels, atmnet, prefixname)

    dims = []
    for i in channels:
        dims.append(i["dim"])
    return dims