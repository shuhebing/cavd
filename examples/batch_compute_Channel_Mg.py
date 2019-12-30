
"""
Calculate the transport channel for Mg-containing CIFs.

"""

import os
from cavd import outVesta
from cavd.cavd_consts import LOWER_THRESHOLD, UPPER_THRESHOLD

migrant = "Mg"

path = "./Li_Na_Mg_Al_cifs_order_revise_recover/"+migrant+"/"

if not os.path.exists(path+"results"):
    os.mkdir(path+"results")
    print("create results directory successful !")
else:
    print(path+"results already exit!")
    
output_path = path+"results/"
com_status = open(output_path+"channel_com_status_" + migrant + ".csv","w")
com_status.write('filename,status\n')
results = open(output_path+"com_results_" + migrant + ".csv","w")
results.write('filename\tDimofChannels\n')

filenames=[]
for i in os.listdir(path):
    if ".cif" in i:
        filenames.append(i)

for filename in filenames:
    filename = path+filename
    print(filename)
    try:
        dims_channel = outVesta(filename, migrant, ntol = 0.02, rad_flag=True, lower=LOWER_THRESHOLD[migrant], upper=UPPER_THRESHOLD[migrant], rad_dict=None)
        results.write(filename)
        results.write('\t')
        results.write(str(dims_channel))
        results.write("\n")
        print(filename+" compute complete!")
        out = filename+'\t'+'compute complete!'+'\n'
        com_status.write(out)
    except Exception as e:
        print(filename,str(e))
        com_status.write(filename+','+str(e)+'\n')
        continue
print(migrant+"-containing files compute completed!")
