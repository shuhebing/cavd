# testAllCom.py
import cavd
symm_symbol_atmnt,symm_num_atmnt,symm_num_vornet,symprec,conn_val,connect,dim_network,dims_channel,migrant_alpha,radii,minRad,nei_dises,recover_rate, recover_state,true_recover_dis,coordination_list= cavd.AllCom5("Li2CO3-LDA.cif",0.798601,"Li",True,True,None)
print("Symmetry symbol in Atom network: ", symm_symbol_atmnt)
print("Symmetry number in Atom network: ", symm_num_atmnt)
print("Symmetry number in Voronoi network: ", symm_num_vornet)
print("Distance tolerance in Cartesian coordinates to find crystal symmetry: ", symprec)

print(radii)

print("Rf_a: ",conn_val[0])
print("Rf_b: ",conn_val[1])
print("Rf_c: ",conn_val[2])

print("a: ", connect[0])
print("b: ", connect[1])
print("c: ", connect[2])

print("The Dimension of the transport network: ", dim_network)
print("The Dimension of the transport channels: ", dims_channel)

print(migrant_alpha)

print(minRad)
print(nei_dises)
print(coordination_list)

print("Recover rate: ", recover_rate)
print(recover_state)
print(true_recover_dis)




