from cavd.netstorage import AtomNetwork
from cavd.graphstorage import DijkstraNetwork
from cavd.channel import Channel
from cavd.netio import writeAtomNetVaspFile

radii = {}
atmnet = AtomNetwork.read_from_CIF("icsd_16713.cif", radii, False, None)
vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
#a = DijkstraNetwork.from_VoronoiNetwork(vornet)
#print(a.lattice)
#print(a.nodes)
atmnet.writeAtomNetVaspFile("icsd_16713.vasp",True)
writeAtomNetVaspFile("icsd_16713.vasp2",atmnet,True)
channels1 = Channel.findChannelsInVornet(vornet,1.1,"icsd_16713.vmd")
channels = Channel.findChannels(vornet,1.1,"icsd_16713.net")

for i in channels:
    print(i.nodes)
    print(i.connections)
    print(i.nodes_deltapos)
    print(i.dimensionality)
    print(i.lattice)