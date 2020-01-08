import cavd
dims, conn_val = cavd.outVesta("icsd_16713.cif", "Li", ntol = 0.02, rad_flag =  True, lower=0.5, upper=10.0, rad_dict=None)
dims, conn_val = cavd.outVesta("icsd_16713-merge.cif", "Li", ntol = 0.5, rad_flag = True, lower=0.5, upper=10.0, rad_dict=None)
dims, conn_val = cavd.outVesta("icsd_246817.cif", "Li", ntol = 0.02, rad_flag =  True, lower=0.5, upper=10.0, rad_dict=None)
dims, conn_val = cavd.outVesta("icsd_246817-merge.cif", "Li", ntol = 0.5, rad_flag = True, lower=0.5, upper=10.0, rad_dict=None)