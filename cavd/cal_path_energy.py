"""
    Calculating the energy profiles according to BVSE calculations.
    Author: Mi Penghui

"""

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
import networkx as nx
import math
import copy
import numpy.linalg as la
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from pymatgen.core.sites import PeriodicSite


def pathfind(start, end, V, n_images=21, dr=None, h=0.1, k=0.17, min_iter=100, max_iter=10000, max_tol=5e-6):
    # Set parameters
    if not dr:
        dr = np.array([1.0 / (V.shape[0] - 1), 1.0 / (V.shape[1] - 1), 1.0 / (V.shape[2] - 1)])
    else:
        dr = np.array(dr, dtype=float)
    keff = k * dr * n_images
    h0 = h
    
    g1 = np.linspace(0, 1, n_images)  
    s0 = start 
    s1 = end  
    s = np.array([g * (s1 - s0) for g in g1]) + s0  
    ds = s - np.roll(s, 1, axis=0) 
    ds[0] = (ds[0] - ds[0])  
    ls = np.cumsum(la.norm(ds, axis=1))  
    ls = ls / ls[-1]  
    fi = interp1d(ls, s, axis=0)  
    s = fi(g1)
    # print(s)
    # Evaluate initial distances (for elastic equilibrium)
    ds0_plus = s - np.roll(s, 1, axis=0)  
    ds0_minus = s - np.roll(s, -1, axis=0) 
    ds0_plus[0] = (ds0_plus[0] - ds0_plus[0])
    ds0_minus[-1] = (ds0_minus[-1] - ds0_minus[-1])
    dV = np.gradient(V)  

    # Evolve string
    for step in range(0, max_iter):
        if step > min_iter:
            h = h0 * np.exp(-2.0 * (step - min_iter) / max_iter)  # 逐步衰减步长以减少震荡
        else:
            h = h0
        # Calculate forces acting on string
        d = V.shape
        s0 = s
        edV = np.array([[dV[0][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]] / dr[0],

                         dV[1][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]] / dr[1],

                         dV[2][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]] / dr[2]] for pt in s])

        # Update according to force due to potential and string elasticity
        ds_plus = s - np.roll(s, 1, axis=0)
        ds_minus = s - np.roll(s, -1, axis=0)
        ds_plus[0] = (ds_plus[0] - ds_plus[0])
        ds_minus[-1] = (ds_minus[-1] - ds_minus[-1])

        Fpot = edV
        Fel = keff * (la.norm(ds_plus) - la.norm(ds0_plus)) * (ds_plus / la.norm(ds_plus))
        Fel += keff * (la.norm(ds_minus) - la.norm(ds0_minus)) * (ds_minus / la.norm(ds_minus))
        s = s - h * (Fpot + Fel)
     
        s[0] = s0[0]
        s[-1] = s0[-1]
        # Reparametrize string         
        ds = s - np.roll(s, 1, axis=0)
        ds[0] = (ds[0] - ds[0])
        ls = np.cumsum(la.norm(ds, axis=1))
        ls = ls / ls[-1]
        fi = interp1d(ls, s, axis=0)
        s = fi(g1)
        tol = la.norm((s - s0) * dr) / n_images / h
        if (tol > 1e9):
            s = [[0, 0, 0]]
            break
        if (step > min_iter and tol < max_tol):
            # print("Converged at step {}".format(step))
            break
        # if (step % 100 == 0):
        # print ("Step {} - ds = {}".format(step, tol))
    return s
def get_dis(struc,p1, p2):
    tmp_siteA = PeriodicSite("Ar",p1,struc.lattice)
    tmp_siteB = PeriodicSite("Ar",p2,struc.lattice)
    return tmp_siteA.distance(tmp_siteB)

def calpointenergy(energy, point):
    a=int(point[0])
    b=int(point[1])
    c=int(point[2])
    roundpoints=[energy[a][b][c]]
    if a+1<energy.shape[0]:
        roundpoints.append(energy[a+1][b][c])
        if b + 1 < energy.shape[1]:
            roundpoints.append(energy[a+1][b + 1][c])
        if c + 1 < energy.shape[2]:
            roundpoints.append(energy[a+1][b ][c+1])
            if b + 1 < energy.shape[1]:
                roundpoints.append(energy[a + 1][b + 1][c+1])
    if b + 1< energy.shape[1]:
        roundpoints.append(energy[a][b+1][c])
        if c + 1 < energy.shape[2]:
            roundpoints.append(energy[a][b+1][c + 1])
    if c + 1 < energy.shape[2]:
        roundpoints.append(energy[a][b][c+1])

    return min(roundpoints)
def pathenergy(energy,struc, p):
    energy = energy - np.amin(energy)
    energy_path = []
    for point_of_p in p:
        point_temp=[0.0,0.0,0.0]
        for i in range(len(point_temp)):
            if point_of_p[i]<0:
                 point_temp[i]=(point_of_p[i]+1) * (energy.shape[i]-1)
            elif point_of_p[i]<1:
                point_temp[i] = point_of_p[i] * (energy.shape[i]-1)
            else:
                point_temp[i] = (point_of_p[i]-1) * (energy.shape[i] - 1)
        energy_point_temp = calpointenergy(energy, point_temp)
        energy_path.append(energy_point_temp)
    dis_path=[0.0]
    tol=0
    for i in range(len(p)-1):
        dist = get_dis(struc,p[i], p[i+1])
        tol+=dist
        dis_path.append(tol)
    migs = np.zeros(shape=(len(dis_path), 2))
    for i in range(len(energy_path)):
        migs[i][1] = energy_path[i]
        migs[i][0] = dis_path[i]
    return migs

def cal_path_energy(file,start_point,end_point,images,degval,step,max_smi):
    struc = Structure.from_file(file+".cif")
    energy = np.load(file+".npy")
    energy = energy - np.amin(energy)
     
    start_f = start_point[1]
    end_f = end_point[1]
    start=np.array([start_f[0]*(energy.shape[0]-1),start_f[1]*(energy.shape[1]-1),start_f[2]*(energy.shape[2]-1)])
    end=np.array([end_f[0]*(energy.shape[0]-1),end_f[1]*(energy.shape[1]-1),end_f[2]*(energy.shape[2]-1)])
    p = pathfind(start, end, energy, n_images=images,
                              dr=[struc.lattice.a / (energy.shape[0]-1),
                                  struc.lattice.b / (energy.shape[1]-1),
                                  struc.lattice.c / (energy.shape[2]-1)], h=step, k=0.17, min_iter=100,
                              max_iter=max_smi, max_tol=5e-6)
    for p1 in p:
        p1[0] = round(p1[0] / (energy.shape[0] - 1), 4)
        p1[1] = round(p1[1] / (energy.shape[1] - 1), 4)
        p1[2] = round(p1[2] / (energy.shape[2] - 1), 4)
   
    p_e= pathenergy(energy,struc, p)
    
    energy_profile = open(file+"_path_"+start_point[0]+"-"+end_point[0]+"_"+str(images)+"_"+str(degval)+"_"+str(step)+"_"+str(max_smi)+".csv", 'w')
    
    for pi in range(len(p)):
        energy_profile.write(str(p[pi])+"\t"+str(p_e[pi])+"\n")
    print(p) 
    print(p_e) 


    for i in range(1,len(p)-1):
        struc.insert(0, 'H', p[i])
    struc.to("poscar",file+"_path_"+start_point[0]+"-"+end_point[0]+"_"+str(images)+"_"+str(degval)+"_"+str(step)+"_"+str(max_smi)+".vasp")
    
    xcoords=[]
    ycoords=[]
    for j in range(len(p_e)):
        xcoords.append(p_e[j][0])
        ycoords.append(p_e[j][1])

    poly = np.polyfit(xcoords, ycoords, deg=degval) 
    z = np.polyval(poly, xcoords)
    plt.figure(figsize=(6, 4))  
    plt.plot(xcoords, z, linewidth=3)
    plt.scatter(xcoords, ycoords, color='k', marker='o')
    plt.xlabel("Reaction Coordinate ") 
    plt.ylabel("Energy") 
    plt.savefig(file+"_path_"+start_point[0]+"-"+end_point[0]+"_"+str(images)+"_"+str(degval)+"_"+str(step)+"_"+str(max_smi)+".png")
    plt.show()
    
if __name__ == "__main__":
    start_point =("Li1",[0.50000, 0.75000, 0.37500])
    end_point = ("Li33",[0.41940, 0.91430, 0.30410])
    cal_path_energy("icsd_246817",start_point,end_point,11,15,0.1,50000)