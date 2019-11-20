"""
    Analyze the recovery rate of target ions.

    Update date：20191120
    Author：YAJ
    School of Computer Engineering and Science, ShangHai University 

"""
from pymatgen.core.sites import PeriodicSite
import numpy as np
from scipy.spatial.ckdtree import cKDTree

# Calculate the recovery rates of the lattice sites of mobile ions.
def rediscovery(migrate, vorosites, stru):
    labels = []
    recover_labels = []
    recover_state = {}
    true_recover_dis = {}
    
    # Site types: interstice, bottleneck, and center of face
    points_type = ["It","Bn","Fc"]


    for k in range(len(stru.sites)):
        site = stru.sites[k]
        label = site._atom_site_label
        if migrate not in label:
            continue
        if label not in labels:
            labels.append(label)
        if label in recover_labels:
            continue

        for pts_idx, pts in enumerate (vorosites):
            cp_tag = np.ones((len(pts), ), dtype=int)
            for pt_idx, pt in enumerate (pts):
                if cp_tag[pt_idx] != -1:
                    print("mobile:",site,"label",label)
                    print("void:",pt)
                    
                    #Use Ar as the temporary symbol for current gap.
                    tmp_site = PeriodicSite("Ar",pt,stru.lattice)
                    print(site.distance(tmp_site))

                    """
                        If a gap has been paired (within 0.5A) with a lattice position in the structure, 
                        the gap and all neighbor gaps (within 0.25A) are removed.
                    """
                    if site.distance(tmp_site) < 0.5:
                      
                        recover_labels.append(label)

                        true_recover_dis[str(label)] = (points_type[pts_idx]+str(pt_idx),site.distance(tmp_site))
                        cp_tag[pt_idx] = -1
                        
                        for pt_idx2, pt2 in enumerate (pts):
                            tmp_site2 = PeriodicSite("Ar",list(pt2),stru.lattice)
                            if tmp_site.distance(tmp_site2) < 0.25:
                                cp_tag[pt_idx2] = -1
                        break

    # Count the recovery rate of the current structure.
    recover_rate = len(recover_labels)/len(labels)
    for la in labels:
        if la in recover_labels:
            recover_state[str(la)] = True
        else:
            recover_state[str(la)] = False

    return recover_rate, recover_state, true_recover_dis

# Return the site type.
def get_point_tag(id, pts_len):
    vexs_len = pts_len[0]
    bts_len = pts_len[1]
    fcs_len = pts_len[2]

    if id < vexs_len:
        return "It" + str(id)
    elif id < vexs_len + bts_len:
        return "Bn" + str(id - vexs_len)
    elif id < vexs_len + bts_len + fcs_len:
        return "Fc" + str(id - vexs_len - bts_len)
    else:
        raise IndexError

def get_point_tag_onlyVertex(id, pts_len):
    vexs_len = pts_len[0]
    if id < vexs_len:
        return "It" + str(id)
    else:
        raise IndexError

# Calculate the recovery rates of the lattice sites of mobile ions by KD-Tree.
def rediscovery_kdTree(stru, migrate, vorosites, threshold = 0.5):
    recover_labels = []
    recover_state = {}
    migrate_mindis = {}
    
    migrate_pos_frac = np.around(np.array([site.frac_coords for site in stru.sites if migrate in site._atom_site_label], ndmin=2), 3)
    migrate_pos_frac %= 1.0
    migrate_pos_frac %= 1.0
    labels = [site._atom_site_label for site in stru.sites if migrate in site._atom_site_label]

    points = np.around(np.array(vorosites[0] + vorosites[1] + vorosites[2], ndmin=2), 3)
    points %= 1.0
    points %= 1.0
    vorositesKdTree = cKDTree(points)
    min_dis,min_ids = vorositesKdTree.query(migrate_pos_frac,k=1)

    for idx in range(len(min_ids)):
        if labels[idx] in recover_labels:
            continue
        tmp_site1 = PeriodicSite("Ar", migrate_pos_frac[idx], stru.lattice)
        tmp_site2 = PeriodicSite("Ar", points[min_ids[idx]], stru.lattice)
        pts_len = [len(vorosites[0]), len(vorosites[1]), len(vorosites[2])]
        pt_tag = get_point_tag(min_ids[idx], pts_len)
        migrate_mindis[str(labels[idx])] = (pt_tag, tmp_site1.distance(tmp_site2))
        if tmp_site1.distance(tmp_site2) <= threshold:
            recover_state[str(labels[idx])] = pt_tag
            recover_labels.append(labels[idx])
        else:
            recover_state[str(labels[idx])] = None

    recover_rate = len(recover_labels) / len(np.unique(labels))
    return recover_rate, recover_state, migrate_mindis

"""
    Find the Interstice, Bottleneck or Face Center that corresponding to the lattice position of mobile ions.
    Algorithm:
        1. Use KD-tree to find the nearest sites (Interstice, Bottleneck or Face Center) of 
        the specified lattice sites of mobile ion.
        2. Record the distance between them. If the distance is smaller than the size of 
        the nearest sites (Interstice, Bottleneck or Face Center), current mobile ions position is considered to be recoveried.
    Return:
        Recovery rate, recovery status of mobile ions, nearest sites (Interstice, Bottleneck or Face Center) and the  distance between them.
"""
def rediscovery_byRad_kdTree(stru, migrate, vorosites, vororad, threshold = 0.5):
    recover_labels = []
    recover_state = {}
    migrate_mindis = {}

    migrate_pos_frac = np.around(np.array([site.frac_coords for site in stru.sites if migrate in site._atom_site_label], ndmin=2), 3)
    migrate_pos_frac %= 1.0
    migrate_pos_frac %= 1.0
    expand_pos_frac = migrate_pos_frac
    # expand the migrant sites to 3*3*3
    for a in range(-1, 2):
        for b in range(-1, 2):
            for c in range(-1, 2):
                if a==b==c==0:
                    continue
                else:
                    expand_pos_frac = np.concatenate((expand_pos_frac,migrate_pos_frac+np.array([a,b,c])),axis=0)

    migrate_labels = [site._atom_site_label for site in stru.sites if migrate in site._atom_site_label]
    
    points = np.around(np.array(vorosites[0] + vorosites[1] + vorosites[2], ndmin=2), 3)
    points_rad = np.array(vororad[0] + vororad[1] + vororad[2])
    
    points %= 1.0
    points %= 1.0

    vorositesKdTree = cKDTree(points)
    min_dis,min_ids = vorositesKdTree.query(migrate_pos_frac,k=1)

    pts_len = [len(vorosites[0]), len(vorosites[1]), len(vorosites[2])]
    for idx in range(len(migrate_labels)):
        if migrate_labels[idx] in recover_labels:
            continue
        tmp_site1 = PeriodicSite("Ar", migrate_pos_frac[idx], stru.lattice)
        tmp_site2 = PeriodicSite("Ar", points[min_ids[idx]], stru.lattice)
        
        pt_tag = get_point_tag(min_ids[idx], pts_len)
        pt_rad = points_rad[min_ids[idx]]
        dis_st1_st2 = tmp_site1.distance(tmp_site2)
        migrate_mindis[str(migrate_labels[idx])] = (pt_tag, pt_rad, dis_st1_st2)
        
        # if dis_st1_st2 <= threshold or dis_st1_st2 <= pt_rad:
        if dis_st1_st2 <= threshold:
            recover_state[str(migrate_labels[idx])] = pt_tag
            recover_labels.append(migrate_labels[idx])
        else:
            recover_state[str(migrate_labels[idx])] = None

    recover_rate = len(recover_labels) / len(np.unique(migrate_labels))
    return recover_rate, recover_state, migrate_mindis

def rediscovery_byRad_kdTree_onlyVertex(stru, migrate, vorosites, vororad, threshold = 0.5):
    recover_labels = []
    recover_state = {}
    migrate_mindis = {}

    migrate_pos_frac = np.around(np.array([site.frac_coords for site in stru.sites if migrate in site._atom_site_label], ndmin=2), 3)
    migrate_pos_frac %= 1.0
    migrate_pos_frac %= 1.0
    expand_pos_frac = migrate_pos_frac
    
    # expand the migrant sites to 3*3*3
    for a in range(-1, 2):
        for b in range(-1, 2):
            for c in range(-1, 2):
                if a==b==c==0:
                    continue
                else:
                    expand_pos_frac = np.concatenate((expand_pos_frac,migrate_pos_frac+np.array([a,b,c])),axis=0)

    migrate_labels = [site._atom_site_label for site in stru.sites if migrate in site._atom_site_label]
    
    points = np.around(np.array(vorosites[0], ndmin=2), 3)
    points_rad = np.array(vororad[0])
    
    points %= 1.0
    points %= 1.0

    vorositesKdTree = cKDTree(points)
    min_dis,min_ids = vorositesKdTree.query(migrate_pos_frac,k=1)

    pts_len = [len(vorosites[0])]
    for idx in range(len(migrate_labels)):
        if migrate_labels[idx] in recover_labels:
            continue
        tmp_site1 = PeriodicSite("Ar", migrate_pos_frac[idx], stru.lattice)
        tmp_site2 = PeriodicSite("Ar", points[min_ids[idx]], stru.lattice)
        
        pt_tag = get_point_tag_onlyVertex(min_ids[idx], pts_len)
        pt_rad = points_rad[min_ids[idx]]
        dis_st1_st2 = tmp_site1.distance(tmp_site2)
        migrate_mindis[str(migrate_labels[idx])] = (pt_tag, pt_rad, dis_st1_st2)
        
        # if dis_st1_st2 <= threshold or dis_st1_st2 <= pt_rad:
        if dis_st1_st2 <= threshold:
            recover_state[str(migrate_labels[idx])] = pt_tag
            recover_labels.append(migrate_labels[idx])
        else:
            recover_state[str(migrate_labels[idx])] = None

    recover_rate = len(recover_labels) / len(np.unique(migrate_labels))
    return recover_rate, recover_state, migrate_mindis