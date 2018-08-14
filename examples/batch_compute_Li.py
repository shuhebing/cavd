import os
from zeo import AllCom
from zeo.netstorage import PerformVDError
from zeo.channel import FindChannelError

filenames=[]
path = "/home/yeanjiang/yaj/bi/Li_Na_Mg_Al_cifs/Li/"
if not os.path.exists(path+"results"):
    os.mkdir(path+"results")
    print("create results directory successful !")
else:
    print(path+"results already exit!")
output_path = path+"results/"
result_file = open(output_path+"com_result_Li.txt","w")
result_file.write('filename\tProblem\n')
Rf_file = open(output_path+"Rf_Li.txt","w")
Rf_file.write('filename\ta_Rf\tb_Rf\tc_Rf\toneD_Conn\ttwoD_Conn\tthreeD_Conn\n')
for i in os.listdir(path):
    if ".cif" in i:
        filenames.append(i)

for filename in filenames:
    filename = path+filename
    try:
        conn,oneD,twoD,threeD = AllCom(filename, 0.584, 1000, migrant="Li", rad_flag=True, effective_rad=True, rad_file=None, rad_store_in_vasp=True, minRad=0.584, maxRad=0.876)
        Rf_file.write(filename)
        for i in conn:
            Rf_file.write('\t'+str(i))
        Rf_file.write('\t'+str(oneD)+'\t'+str(twoD)+'\t'+str(threeD))
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
print("All File compute completed!")
