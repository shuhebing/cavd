# -*- coding: UTF-8 -*-
import os
import re
from cavd.netio import *
from monty.io import zopen
from cavd.channel import Channel
from cavd.netstorage import AtomNetwork, connection_values_list
from cavd.local_environment import CifParser_new, LocalEnvirCom

def outVesta(filename, migrant, ntol=0.02, rad_flag=True, lower=0.0, upper=10.0, rad_dict=None):
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
	
    prefixname = filename.replace(".cif","")

    # compute the R_T
    conn_val = connection_values_list(prefixname+".resex", vornet)
    channels = Channel.findChannels2(vornet, atmnet, lower, upper, prefixname+".net")
    
    # output vesta file for visiualization
    Channel.writeToVESTA(channels, atmnet, prefixname)
    
    return conn_val

if __name__ == "__main__":
    # Calculate transport channel using given radii.
    rad_dict = {"Li1": 0.73, "C1": 0.06, "O1": 1.24, "O2": 1.23}
    conn_val = outVesta("icsd_16713.cif","Li",ntol=0.02, rad_flag=True, lower=0.0, upper=10.0, rad_dict)
    print(conn_val)


