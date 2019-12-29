
"""
Calculate the transport channel for Na-containing CIFs.

"""

import os
from cavd import outVesta
from cavd.netstorage import PerformVDError
from cavd.channel import FindChannelError
from cavd.local_environment import CoordComError
from cavd.cavd_consts import LOWER_THRESHOLD, UPPER_THRESHOLD

migrant = "Na"

path = "./Li_Na_Mg_Al_cifs_order_revise_recover/"+migrant+"/"

if not os.path.exists(path+"results"):
    os.mkdir(path+"results")
    print("create results directory successful !")
else:
    print(path+"results already exit!")
    
output_path = path+"results/"
result_file = open(output_path+"com_result_"+migrant+".txt","w")
result_file.write('filename\tProblem\n')
Rf_file = open(output_path+"Dim_"+migrant+".txt","w")
Rf_file.write('filename\tchannels_dim\n')

filenames=[]
for i in os.listdir(path):
    if ".cif" in i:
        filenames.append(i)

for filename in filenames:
    filename = path+filename
    print(filename)
    try:
        dims = outVesta(filename, migrant, ntol = 0.02, rad_flag=True, lower=LOWER_THRESHOLD[migrant], upper=UPPER_THRESHOLD[migrant], rad_dict=None)
        Rf_file.write(filename)
        for value in dims:
            Rf_file.write(str(value)+" ")
        Rf_file.write("\n")
        print(filename+" compute complete!")
        out = filename+'\t'+'compute complete!'+'\n'
        result_file.write(out)
    except AttributeError:
        print(filename," Have PARTITIAL occ or MIXED occ!")
        out = filename+'\t'+"Have PARTITIAL occ or MIXED occ!"+'\n'
        result_file.write(out)
        continue
    except IOError:
        print(filename, " Can't Open inputfile or Can't Write to outputfile.")
        out = filename+'\t'+"Can't Open inputfile or Can't Write to outputfile."+'\n'
        result_file.write(out)
        continue
    except PerformVDError:
        print(filename, " Can't Perform Voronoi Decompition.")
        out = filename+'\t'+"Can't Perform Voronoi Decompition."+'\n'
        result_file.write(out)
        continue
    except ValueError:
        print(filename, " Have MIXED occ!")
        out = filename+'\t'+" Have MIXED occ!"+'\n'
        result_file.write(out)
        continue
    except FindChannelError:
        print(filename, " Can't Find Channel in Voronoi Network.")
        out = filename+'\t'+"Can't Find Channel in Voronoi Network."+'\n'
        result_file.write(out)
        continue
    except KeyError:
        print(filename, " Compute radius failed when search radius information from  Shannon effective ionic radius table.")
        out = filename+'\t'+"Compute radius failed when search radius information from  Shannon effective ionic radius table."+'\n'
        result_file.write(out)
        continue
    except AssertionError:
        print(filename, " Compute radius failed when pymatgen try to read information.")
        out = filename+'\t'+"Compute radius failed when pymatgen try to read information."+'\n'
        result_file.write(out)
        continue
    except UnboundLocalError:
        print(filename, " Compute radius failed when pymatgen try to compute coordnum.")
        out = filename+'\t'+"Compute radius failed when pymatgen try to compute coordnum."+'\n'
        result_file.write(out)
        continue
    except ZeroDivisionError:
        print(filename, " Integer division or modulo by zero.")
        out = filename+'\t'+"Integer division or modulo by zero."+'\n'
        result_file.write(out)
        continue
    except CoordComError:
        print(filename, " Computer Coord Failed! Maybe have no CN!")
        out = filename+'\t'+"Computer Coord Failed! Maybe have no CN!"+'\n'
        result_file.write(out)
        continue

    except IndexError:
        print(filename, "IndexError: list index out of range! Maybe occuied in local_environment compute neighbor!")
        out = filename+'\t'+"Maybe occuied in local_environment compute neighbor!"+'\n'
        result_file.write(out)
        continue

print("All File compute completed!")
