import os
import sys
from zeo.netstorage import AtomNetwork
from zeo.netio import *

#计算某个结构的瓶颈和间隙
def BIComputation(filename, migrant=None, rad_flag=True, rad_file=None,rad_store_in_vasp=True, minRad=0.0):
	if migrant:
		remove_filename = getRemoveMigrantFilename(filename,migrant)
	else:
		renove_filename = filename
	atmnet = AtomNetwork.read_from_CIF(remove_filename, rad_flag, rad_file)
	vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition()
	prefixname = filename.replace(".cif","")
	writeBIfile(prefixname+".bi",vornet,minRad)
	writeVaspFile(prefixname+".vasp",atmnet,vornet,rad_store_in_vasp,minRad)
	
#计算某个机构的连通性
def Connection(filename, probe_size, migrant=None, rad_flag=True, rad_file=None):
	if migrant:
		remove_filename = getRemoveMigrantFilename(filename,migrant)
	else:
		remove_filename = filename
	atmnet = AtomNetwork.read_from_CIF(remove_filename, rad_flag, rad_file)
	prefixname = filename.replace(".cif","")
	return atmnet.through_VorNet(prefixname+".res",probe_size)

#batch read file function
def BatchReadFilename(path,filetype):
    filenames=[]
    for path_root,path_directory_name,files in os.walk(path):
        for i in files:
            if filetype in i:
                filenames.append(i.replace(filetype,''))
    return filenames

#批处理计算瓶颈和间隙
def BIComputation_batch(path, migrant=None, rad_flag=True, rad_file=None,rad_store_in_vasp=True, minRad=0.01):
	filenames = BatchReadFilename(path,".cif")
	print(filenames)
	if migrant:	
		for filename in filenames:
			remove_filename = getRemoveMigrantFilename(filename+".cif",migrant)
			atmnet = AtomNetwork.read_from_CIF(remove_filename, rad_flag, rad_file)
			vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition()
			writeBIfile(filename+".bi",vornet,minRad)
			writeVaspFile(filename+".vasp",atmnet,vornet,rad_store_in_vasp,minRad)
	else:	
		for filename in filenames:
			remove_filename = filename+".cif"
			atmnet = AtomNetwork.read_from_CIF(remove_filename, rad_flag, rad_file)
			vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition()
			writeBIfile(filename+".bi",vornet,minRad)
			writeVaspFile(filename+".vasp",atmnet,vornet,ad_store_in_vasp,minRad)

#批处理计算连通性
def Connection_batch(path, probe_size, migrant=None, rad_flag=True, rad_file=None):
	filenames = BatchReadFilename(path,".cif")
	result_file = open(path+"result.txt","w")
	if migrant:	
		for filename in filenames:
			remove_filename = getRemoveMigrantFilename(path+filename+".cif",migrant)
			atmnet = AtomNetwork.read_from_CIF(path+remove_filename, rad_flag, rad_file)
			result_file.write(filename+".cif    "+atmnet.through_VorNet(filename+".res",probe_size))
	else:	
		for filename in filenames:
			remove_filename = filename+".cif"
			atmnet = AtomNetwork.read_from_CIF(remove_filename, rad_flag, rad_file)
			result_file.write(filename+".cif    "+str(atmnet.through_VorNet(filename+".res",probe_size)))

"""
def main():
	BIComputation("Li2CO3-LDA.cif",migrant="Li")
	Connection("Li2CO3-LDA.cif",0.5)
	#BIComputation_batch("./cifs/",migrant="Li")
	#Connection_batch("./cifs/",0.5)
if __name__ == "__main__":
	sys.exit(main())
"""
