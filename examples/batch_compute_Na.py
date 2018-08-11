import os
from zeo import AllCom
from zeo.netstorage import PerformVDError
from zeo.channel import FindChannelError

filenames=[]
path = "/home/yeanjiang/yaj/bi/Li_Na_Mg_Al_cifs/Na/"
if not os.path.exists(path+"results"):
    os.mkdir(path+"results")
    print("create results directory successful !")
else:
    print(path+"results already exit!")
output_path = path+"results/"
result_file = open(output_path+"computation_result.txt","w")
result_file.write("filename    a_Rf    b_Rf    c_Rf\n")
for i in os.listdir(path):
    if ".cif" in i:
        #filenames.append(i.replace(filetype,''))
        filenames.append(i)
output_path = path+"results/"
for filename in filenames:
    filename = path+filename
    try:
        conn,oneD,twoD,threeD = AllCom(filename, 0.5, 1000, migrant="Li", rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True, minRad=0.904, maxRad=1.356)
        result_file.write(filename)
        for i in conn:
            result_file.write("    "+str(i))
        result_file.write("\n")
        print(filename+" compute complete1!")
    except AttributeError:
        print("This cif file: ",filename,"have PARTITIAL occ or MIXED occ!")
        continue
    except IOError:
        print("Can't Open ", filename, " or Can't Write to outputfile.")
        continue
    except PerformVDError:
        print("Can't Perform Voronoi Decompition for ", filename)
        continue
    except ValueError:
        print("Can't Read structure information from cif file: ", filename)
        continue
    except FindChannelError:
        print("Find Channel in Voronoi Network for cif file: ", filename)
        continue
    except KeyError:
        print("Compute radius failed when search radius information from table. problem file: ", filename)
        continue
    except AssertionError:
        print("Compute radius failed when pymatgen try to read information from: ", filename)
        continue
    except UnboundLocalError:
        print("Compute radius failed when pymatgen try to compute coordnum from: ", filename)
        continue
print("All File compute completed!")
