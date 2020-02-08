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
    
    allPaths = migratPath.comAllPaths(prefixname)
    
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

# The function used to calculate the ion transport descriptors.
def comDescriptors(filename, migrant, rad_flag=True, lower=0.0, upper=10.0, rad_dict=None, ntol=0.02):
    # Reading cif.
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)

    # Reading Symmetry_Int_Tables_number and Symmetry_space_group_name_H-M from ICSD cif.
    symm_number,symm_sybol = parser.get_symme()
   
    # Calculating ionic radii (fron Shannon effective ionic radii table), coordination environment (coordination number, minimal mobile-framework ions distance).
    stru = parser.get_structures(primitive=False)[0]
    effec_radii, nei_dises = LocalEnvirCom_new(stru, migrant)

    # Checking whether the structure contains given mobile ion.
    species = [str(sp).replace("Specie ","") for sp in stru.species]
    elements = [re.sub('[^a-zA-Z]','',sp) for sp in species]
    sitesym = parser.get_sym_opt()
    if migrant not in elements:
        raise ValueError("The input migrant ion not in the input structure! Please check it.")

    # Constructing a Atom network from the structure that mobile ion has been removed.
    radii = {}
    if rad_flag:
        if rad_dict:
            radii = rad_dict
        else:
            radii = effec_radii
    atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii, rad_flag)
    
    # Constructing a Voronoi network for given Atom network.
    vornet,edge_centers,fcs,faces = atmnet.perform_voronoi_decomposition(True, ntol)

    # Adding Voronoi face center into network to obatin the interstitial network.
    add_fcs_vornet = vornet.add_facecenters(faces)

    # Marking the calculated voids by Symmetry operation.
    sitesym = parser.get_sym_opt()
    sym_vornet,voids = get_labeled_vornet(add_fcs_vornet, sitesym, ntol)

    # Count the recovery rate by distance.
    voids_abs = []
    voids_rad = []
    for void in sym_vornet.nodes:
        voids_abs.append(void[2])
        voids_rad.append(void[3])

    bottlenecks = []
    bottlenecs_rad = []
    for bt in sym_vornet.edges:
        frac_bt = bt[2]
        bottlenecks.append(frac_bt)
        bottlenecs_rad.append(bt[3])

    fcens = []
    fcens_rad = []
    for fc in fcs:
        fcens.append(fc[0])
        fcens_rad.append(fc[1])

    vorosites = [voids_abs, bottlenecks, fcens]
    vororad = [voids_rad, bottlenecs_rad, fcens_rad]
    recover_rate, recover_state, migrate_mindis = rediscovery_byRad_kdTree(stru, migrant, vorosites, vororad)

    # Calculating connected threshold and dim of the interstitial network.
    prefixname = filename.replace(".cif","")
    conn_val = connection_values_list(prefixname+".resex", sym_vornet)
    dim_network,connect = ConnStatus(conn_val, lower, upper)
    
    # Calculating transport channels and dim of the transport channels.
    channels = Channel.findChannels2(sym_vornet,atmnet,lower,upper,prefixname+".net")
    dims_channel = []
    if len(channels)==0:
        pass
    else:
        for i in channels:
            dims_channel.append(i["dim"])
  
    return radii,symm_sybol,symm_number,nei_dises,conn_val,connect,dim_network,dims_channel,recover_rate
    
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
    channels = Channel.findChannels2(sym_vornet,atmnet,lower,upper,prefixname+".net")

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