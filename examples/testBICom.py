# testBICom.py
import cavd
rad_dict = {"Li1": 0.76, "C1": 0.04, "O1": 1.26, "O2": 1.24}
cavd.BIComputation("icsd_16713.cif", "Li", True, 0.5, 0.7, rad_dict)
