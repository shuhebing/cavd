import cavd
radii, minRad, conn_val, connect, dim_network, dims_channel, migrate_mindis = cavd.outVesta("icsd_246817.cif", "Li", True, lower=0.5, upper=10.0, rad_dict=None, symprec=0.01)
print("radii: ",radii)

print("RT_a: ",conn_val[0])
print("RT_b: ",conn_val[1])
print("RT_c: ",conn_val[2])

print("a: ", connect[0])
print("b: ", connect[1])
print("c: ", connect[2])

print("The dimension of the transport network: ", dim_network)
print("The dimensions of the transport channels: ", dims_channel)

print("migrate_mindis: ", migrate_mindis)




