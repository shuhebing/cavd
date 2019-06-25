import cavd
conn_val = cavd.BIComputation("Li2CO3-LDA.cif", "Li", True)

print("Rf_a: ",conn_val[0])
print("Rf_b: ",conn_val[1])
print("Rf_c: ",conn_val[2])





