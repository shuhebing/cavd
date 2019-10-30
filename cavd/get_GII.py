"""
Calculate GII

"""

from cavd.BVAnalysis import BVAnalysis
import numpy as np
from math import sqrt
import re
from cavd.local_environment import VoronoiNN_self
from cavd.local_environment import CifParser_new
from pymatgen.core.structure import Structure

from monty.io import zopen
import spglib
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def get_GII(filename, migrant):
    bvcla = BVAnalysis()
    bvcla.SetMoveIon(migrant)
    dM_sum = 0
    labels = []
    vnn = VoronoiNN_self(cutoff = 10.0)
    stru = Structure.from_file(filename)
    if not filename[-4:] == '.cif':    
        stru.add_oxidation_state_by_element(oxi_dic)

    #print(stru)
    for i in range(len(stru.sites)):
        site = stru.sites[i]
        # print(site.specie)
        coord = site.coords
        coord_no = vnn.get_cn(stru, i)

        neighbors = []
        oxi = int(site.specie.oxi_state)
        results = stru.get_neighbors(site, r=5.5)
        results = [(i[0],i[1]) for i in sorted(results, key=lambda s: s[1])]
        #print('center: ', (results))
        for i in range(len(results)):
            if (results[i][0].specie.oxi_state*oxi < 0):
                neighbors.append(results[i])

        bvs = 0
        for i in range(len(neighbors)):
            distance = neighbors[i][1]
            neig = neighbors[i][0].specie
            oxi_2 = int(neighbors[i][0].specie.oxi_state)
            #print(neig)
            #print('distance: ',distance)
            cent = re.sub('[^a-zA-Z]','',str(site.specie))
            neig = re.sub('[^a-zA-Z]','',str(neig))
            if oxi > 0:
                key="".join([str(cent),str(oxi),str(neig),str(oxi_2)])
            else:
                key="".join([str(neig),str(oxi_2),str(cent),str(oxi)])
            #print('key',key)
            if key in bvcla._BVparam:
                (r,b)=bvcla._BVparam[key][0]
                bv=np.exp((r-distance)/b)
                #print('bv:',bv)
                bvs = bv + bvs
        # print('bvs:',bvs)
        dM = bvs - abs(oxi)
        # print('dM:',dM)
        dM_sum = dM_sum + dM**2
        
    length = len(stru)
    #print(length)
    GII = sqrt(dM_sum/length)
    return GII

 


