import os
import cavd

#数据挖掘获取的标准值
std_surf_list = [0.798601]
std_list = [2.068571]
ions = ["Li"]

for i in range(len(ions)):
    filenames=[]
    path = "./" + ions[i] + "/"
    if not os.path.exists(path+"results"):
        os.mkdir(path+"results")
        print("create results directory successful!")
    else:
        print(path+"results already exit!")
    output_path = path+"results/"
    com_status = open(output_path+"com_status_" + ions[i] + ".csv","w")
    com_status.write('filename,status\n')
    results = open(output_path+"results_" + ions[i] + ".csv","w")
    results.write('filename\tsymmSymbolInCif\tsymmNumInCif\tsym_opt_num\tvoids_num\tsymprec\tuni_voids_num\trecover_rate\trecover_state\tmigrate_mindis\n')
    for j in os.listdir(path):
        if ".cif" in j:
            filenames.append(j)

    for filename in filenames:
        filename = path+filename
        print(filename)
        try:
            symm_symbol_atmnt,symm_num_atmnt,symprec,voids_num,sym_opt_num,uni_voids_num,recover_rate,recover_state,migrate_mindis = cavd.AllCom8(filename,std_surf_list[i],ions[i],True,True,None)
            results.write(filename)
            results.write('\t')
            results.write(symm_symbol_atmnt)
            results.write('\t')
            results.write(str(symm_num_atmnt))
            results.write('\t')
            results.write(str(sym_opt_num))
            results.write('\t')
            results.write(str(voids_num))
            results.write('\t')
            results.write(str(symprec))
            results.write('\t')
            results.write(str(uni_voids_num))
            results.write('\t')
            results.write(str(recover_rate))
            results.write('\t')
            results.write(str(recover_state))
            results.write('\t')
            results.write(str(migrate_mindis))
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
