import os
import re
from monty.io import zopen
from cavd.channel import Channel
from cavd.cavd_consts import LOWER_THRESHOLD, UPPER_THRESHOLD
from cavd.netstorage import AtomNetwork, connection_values_list
from cavd.recovery import rediscovery_byRad_kdTree
from cavd.get_Symmetry import get_labeled_vornet
from cavd.local_environment import CifParser_new, LocalEnvirCom_new
from cavd import ConnStatus

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
    
    # output *.vesta file for visiualization.
    Channel.writeToVESTA(channels, atmnet, prefixname)
  
    return radii,symm_sybol,symm_number,nei_dises,conn_val,connect,dim_network,dims_channel,recover_rate

def batchCom(dir, ions):
    for i in range(len(ions)):
        filenames=[]
        ionDir = dir + ions[i] + "/"
        if not os.path.exists(ionDir + "results"):
            os.mkdir(ionDir + "results")
        
        resultsDir = ionDir + "results/"
        comStatusFile = open(resultsDir + "comstatus_" + ions[i] + ".csv", "w")
        resultFile = open(resultsDir + "com_results_" + ions[i] + ".csv", "w")
        comStatusFile.write('Filename,status\n')
        resultFile.write('Filename\tSymmetry_space_group_name_H-M\tRT_a\tRT_b\tRT_c \tIND\tTCD\tAcc_a\tAcc_b\tAcc_c\tRecovery_status\tIonic_radii\t Symmetry_Int_Tables_number\tMinimal_mobile-framework_distance\n')
        
        for j in os.listdir(ionDir):
            if ".cif" in j:
                filenames.append(j)

        for filename in filenames:
            filename = ionDir + filename
            print(filename)
            try:
                radii,symm_sybol,symm_number,nei_dises,conn_val,connect,dim_network,dims_channel,recover_rate = comDescriptors(filename, ions[i], True, LOWER_THRESHOLD[ions[i]], UPPER_THRESHOLD[ions[i]])
                resultFile.write(filename)
                resultFile.write('\t')
                resultFile.write(symm_sybol)
                resultFile.write('\t')
                resultFile.write(str(conn_val[0]))
                resultFile.write('\t')
                resultFile.write(str(conn_val[1]))
                resultFile.write('\t')
                resultFile.write(str(conn_val[2]))
                resultFile.write('\t')
                resultFile.write(str(dim_network))
                resultFile.write('\t')
                resultFile.write(str(dims_channel))
                resultFile.write('\t')
                resultFile.write(str(connect[0]))
                resultFile.write('\t')
                resultFile.write(str(connect[1]))
                resultFile.write('\t')
                resultFile.write(str(connect[2]))
                resultFile.write('\t')
                resultFile.write(str(recover_rate))
                resultFile.write('\t')
                resultFile.write(str(radii))
                resultFile.write('\t')
                resultFile.write(str(symm_number))
                resultFile.write('\t')
                resultFile.write(str(nei_dises))
                resultFile.write("\n")

                comStatusFile.write(filename+','+'Success!'+'\n')
            except Exception as e:
                comStatusFile.write(filename+','+str(e)+'\n')
                continue
        print(ions[i]+" contained file compute completed!")
    print("All File compute completed!")

if __name__ == "__main__":
    ions = ["Li", "Na", "Mg", "Al"]
    dir = "cifs/"
    batchCom(dir, ions)
