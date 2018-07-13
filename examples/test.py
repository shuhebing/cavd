import zeo
import time
#from get_radii import get_ionic_radii
#zeo.Computation_batch("./tests/cifs/",migrant="Li",minRad=0.2,maxRad=0.5)
#time_start = time.time()
zeo.com("./icsd_16713.cif")
#zeo.Computation("./tests/cifs/1+3 (2).cif",migrant="Li",rad_flag=True, rad_file=None, rad_store_in_vasp=True, minRad=0.76, maxRad=2.0)
#zeo.Computation("./tests/cifs/beta-Li3PS4.cif",migrant="Li",rad_flag=True, rad_file=None, rad_store_in_vasp=True, minRad=0.76, maxRad=2.0)
#time_end = time.time()
#print(time_end-time_start)
#zeo.Computation("./tests/cifs/Li2CO3-LDA.cif",migrant="Li",rad_flag=True, rad_file=None, rad_store_in_vasp=True, minRad=0.78, maxRad=5.0)
#zeo.Computation("./tests/cifs/icsd_467.cif", migrant="Na", rad_flag=True, rad_file=None, rad_store_in_vasp=True, minRad=0.5, maxRad=1.0)
#zeo.Computation("./tests/cifs/icsd_467copy1.cif", minRad=0.5, maxRad=1.0)
#zeo.Computation("./tests/cifs/icsd_467_removeNa_copy.cif", minRad=0.5, maxRad=1.0)
#zeo.Connection("./tests/cifs/Li2CO3-LDA.cif",0.5)
#zeo.Computation_batch("./tests/cifs/",migrant="Li")
#radii_dict = get_ionic_radii("./tests/cifs/Li2CO3-LDA.cif")
#print(radii_dict)
