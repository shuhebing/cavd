
from zeo.netstorage import AtomNetwork
from zeo.netio import *

remove_filename = getRemoveMigrantFilename("Li2CO3-LDA.cif","Li")
atmnet = AtomNetwork.read_from_CIF(remove_filename,False)
vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition()
writeBIfile("Li2CO3-LDA.bi",vornet,0.2)
writeVaspFile("Li2CO3-LDA.vasp",atmnet,vornet,True,0.2)
if atmnet.through_VorNet("Li2CO3-LDA.res",0.5):
	print("a sphere with radius id 0.5A can free through in this structure!")
else:
	print("a sphere with radius id 0.5A can not free through in this structure!")
