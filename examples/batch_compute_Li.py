"""
目前存在的问题：
1. 异常处理问题:
	a.异常cif文件处理（位置不全，含Li、C、O但仅给出Li的位置），输出错误信息”cif文件内容不全“。
	b.含部分占据、混占（无法计算配位数、半径），输出错误信息“含部分占据、混占，无法计算配位数”。
	c.非正常价态以及小数价态（无法从香农表中找到对应的价态，继而无法计算出半径），输出错误信息“非正常价态或非整数价态”。
2. 输出文件存放位置问题。
3. cif文件中价态表示问题：Ag+、Ag1+问题。
4. 接口问题，仔细思考设计对外提供的接口格式。
"""
from zeo import Computation_new
from zeo import com

filenames=[]
for i in os.listdir("/home/yeanjiang/yaj/bi/Li_Na_Mg_Al_cifs/Li"):
    if ".cif" in i:
        #filenames.append(i.replace(filetype,''))
        filenames.append(i)
output_path = path+"results/"
for filename in filenames:
    filename = path+filename
    try:
        conn,oneD,twoD,threeD = AllCom(filename, probe_rad, num_sample, migrant=None, rad_flag=True, effective_rad=False, rad_file=None, rad_store_in_vasp=True, minRad=0.0, maxRad=0.0)
        print(filename+" compute complete1!")
    except AttributeError:
        print("This cif file: ",filename,"have PARTITIAL occ or MIXED occ!")
    except IOError:
        print("Can't Open ", filename, " or Can't Write to outputfile.")
    except PerformVDError:
        print("Can't Perform Voronoi Decompition for ", filename)
    continue
    
    if not os.path.exists(path+"results"):
        os.mkdir(path+"results")
        print("create results directory successful !")
    else:
        print(path+"results already exit!")
    output_path = path+"results/"
    result_file = open(output_path+"computation_result.txt","w")
    result_file.write(filename)
    for i in conn:
       result_file.write("    "+i)
    result_file.write("\n")
    print(filename+" compute completed!")
print("All File compute completed!")
