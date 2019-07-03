import cavd
radii, minRad, conn_val, connect, dim_network, dims_channel, migrate_mindis = cavd.bmd_com("icsd_009004.cif", "Li", True)

print("radii: ",radii)
print("channel threshold: ", minRad)

print("Rf_a: ",conn_val[0])
print("Rf_b: ",conn_val[1])
print("Rf_c: ",conn_val[2])

print("a: ", connect[0])
print("b: ", connect[1])
print("c: ", connect[2])

print("The dimension of the transport network: ", dim_network)
print("The dimension of the transport channels: ", dims_channel)

print("migrate_mindis: ", migrate_mindis)




