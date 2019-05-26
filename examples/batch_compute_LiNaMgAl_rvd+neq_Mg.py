import os
import cavd

#数据挖掘获取的标准值
std_surf_list = [0.827138]
std_list = [2.062649]
ions = ["Mg"]

for i in range(len(ions)):
    filenames=[]
    path = "../../Li_Na_Mg_Al_rm_alloy/" + ions[i] + "/"
    if not os.path.exists(path+"results"):
        os.mkdir(path+"results")
        print("create results directory successful!")
    else:
        print(path+"results already exit!")
    output_path = path+"results/"
    com_status = open(output_path+"com_status_" + ions[i] + ".csv","w")
    com_status.write('filename,status\n')
    results = open(output_path+"results_" + ions[i] + ".csv","w")
    results.write('filename\tsymmSymbolInCif\tsymmNumInCif\tradii\tmin_dises\tsymprec\tsymmNumInVornet\tRT_a\tRT_b\tRT_c\talpha\tR_j\ta\tb\tc\tdimofNetwork\tdimofChannels\trecover_rate\trecover_state\ttrue_recover_dis\tCoordinates\n')
    for j in os.listdir(path):
        if ".cif" in j:
            filenames.append(j)

    for filename in filenames:
        filename = path+filename
        print(filename)
        try:
            symm_symbol_atmnt,symm_num_atmnt,symm_num_vornet,symprec,conn_val,connect,dim_network,dims_channel,migrant_alpha,radii,minRad,nei_dises,recover_rate,recover_state,true_recover_dis,coordination_list = cavd.AllCom5(filename,std_surf_list[i],ions[i],True,True,None)
            results.write(filename)
            results.write('\t')
            results.write(symm_symbol_atmnt)
            results.write('\t')
            results.write(str(symm_num_atmnt))
            results.write('\t')
            results.write(str(radii))
            results.write('\t')
            results.write(str(nei_dises))
            results.write('\t')
            results.write(str(symprec))
            results.write('\t')
            results.write(str(symm_num_vornet))
            results.write('\t')
            results.write(str(conn_val[0]))
            results.write('\t')
            results.write(str(conn_val[1]))
            results.write('\t')
            results.write(str(conn_val[2]))
            results.write('\t')
            results.write(str(migrant_alpha))
            results.write('\t')
            results.write(str(minRad))
            results.write('\t')
            results.write(str(connect[0]))
            results.write('\t')
            results.write(str(connect[1]))
            results.write('\t')
            results.write(str(connect[2]))
            results.write('\t')
            results.write(str(dim_network))
            results.write('\t')
            results.write(str(dims_channel))
            results.write('\t')
            results.write(str(recover_rate))
            results.write('\t')
            results.write(str(recover_state))
            results.write('\t')
            results.write(str(true_recover_dis))
            results.write('\t')
            whole_coords = []
            for coord in coordination_list:
                nei_list = []
                for nei in coord["coord_nei"]:
                    near_tmp = (nei[0].species_string,nei[1])
                    nei_list.append(near_tmp)
                tmp = (coord["label"],nei_list)
                whole_coords.append(tmp)
            results.write(str(whole_coords))
            results.write("\n")
            print(filename+" compute complete!")
            out = filename+','+'compute complete!'+'\n'
            com_status.write(out)
        except Exception as e:
            print(filename,str(e))
            com_status.write(filename+','+str(e)+'\n')
            continue
    print(ions[i]+" contained file compute completed!")
    print()
print("All File compute completed!")
