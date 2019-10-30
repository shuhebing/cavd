import os
import cavd
from monty.io import zopen
from cavd.local_environment import CifParser_new
#数据挖掘获取的标准值
ions = ["Li", "Na", "Mg", "Al"]

for i in range(len(ions)):
    filenames=[]
    path = "../../Li_Na_Mg_Al_cifs_order_revise/" + ions[i] + "/"
    if not os.path.exists(path+"results"):
        os.mkdir(path+"results")
        print("create results directory successful!")
    else:
        print(path+"results already exit!")
    output_path = path+"results/"
    com_status = open(output_path+"com_status_" + ions[i] + ".csv","w")
    com_status.write('filename,status\n')
    results = open(output_path+"results_" + ions[i] + ".csv","w")
    results.write('filename\tradii\tmigrant_radius\tmigrant_alpha\tnei_dises\tcoordination_list\n')
    for j in os.listdir(path):
        if ".cif" in j:
            filenames.append(j)

    for filename in filenames:
        filename = path+filename
        print(filename)
        radii = {}
    
        try:
            with zopen(filename, "rt") as f:
                input_string = f.read()
            parser = CifParser_new.from_string(input_string)
            stru = parser.get_structures(primitive=False)[0]

            radii,migrant_radius,migrant_alpha,nei_dises,coordination_list = cavd.LocalEnvirCom(stru,ions[i])
            results.write(filename)
            results.write('\t')
            results.write(str(radii))
            results.write('\t')
            results.write(str(migrant_radius))
            results.write('\t')
            results.write(str(migrant_alpha))
            results.write('\t')
            results.write(str(nei_dises))
            results.write('\t')
            results.write(str(coordination_list))
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
