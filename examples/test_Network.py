from cavd.netstorage import AtomNetwork
from cavd.netstorage import Atom
from cavd.netstorage import VoronoiNode
from cavd.netstorage import VoronoiEdge

a = Atom()
print(a.cart_coords)
print(a.frac_coords)
print(a.radius)
print(a.atom_type)
print(a.label)
print(a.specialID)
print(a.mass)
print(a.charge)

radii = {}
atmnet = AtomNetwork.read_from_CIF("icsd_16713.cif", radii, False, None)
print(atmnet.lattice_para)
print(atmnet.lattice_angle)
print(atmnet.lattice)
print(atmnet.atoms_num)
print(atmnet.atoms)

vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
print(vornet.nodes)
print(vornet.edges)

channels = Channel.findChannels(vornet,1.1,"icsd_16713.net")
for channel in channels:
    print(channel.nodes)
    print(channel.connections)
    print(channel.nodes_deltapos)