from zeo.netstorage import AtomNetwork
from zeo.netstorage import Atom
from zeo.netstorage import VoronoiNode
from zeo.netstorage import VoronoiEdge

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

b = VoronoiNode()
print(b.coords)
print(b.radius)

c = VoronoiEdge()
print(c.origin)
print(c.ending)
print(c.radius)
print(c.leng)
print(c.delta_uc)
print(c.bot_coords)

vornet,edge_centers,fcs = atmnet.perform_voronoi_decomposition(False)
print(vornet.nodes)
print(vornet.edges)