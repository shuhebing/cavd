import os
import sys
from cavd.high_accuracy import high_accuracy_atmnet
from cavd.netstorage import AtomNetwork
from cavd.netio import *
from matplotlib.font_manager import path

#batch read .cif filename
def batch_read_filename(path,filetype):
    filenames=[]
    for path_root,path_directory_name,files in os.walk(path):
        for i in files:
            if filetype in i:
                filenames.append(i.replace(filetype,''))
    return filenames

def main():
    path="./cifs/"
    filetype=".cif"
    filenames=batch_read_filename(path, filetype)
    for file in filenames:
        whole_filename=file+".cif"
        whole_radfilename=file+".rad"
        out_whole_file=file+".res"
        infile=path+whole_filename
        inradfile=path+whole_radfilename
        outfile=path+out_whole_file
        remove_filename = getRemoveMigrantFilename(infile,"Li")
        atmnet = AtomNetwork.read_from_CIF(remove_filename,rad_file=inradfile)
        atmnet.calculate_free_sphere_parameters(outfile)

if __name__ == "__main__":
    sys.exit(main())