#Define constants for cavd
"""
 Lower threshold (LOWER_THRESHOLD) and upper threshold (UPPER_THRESHOLD) used to determine the ionic accessibility (for Li+, Na+, Mg2+, and Al3+). The criteria is determined by calculating the distances from the lattice sites of mobile ions to their framework atoms surface were calculated in 12,448 coordination environments. A kernel density estimate plot (Figure 4) is formed by computing a continuous probability distribution estimate (using Gaussian kernels with automatic bandwidth determination68,69) of the minimal distances. The calculated distances are distributed within a well-defined range, the thresholds are obtained from about 90% of the data. These thresholds are then used to determine the mobile ions accessibility of the voids in Li-, Na-, Mg- and Al-containing compounds, respectively. If the size r of the interstice or bottleneck satisfies the equation Tl ≤ r ≤ Tu, it is considered to be accessible.

"""

LOWER_THRESHOLD = {"Li":0.5267, "Na":0.9295, "Mg":0.5513, "Al":0.3447}
UPPER_THRESHOLD = {"Li":0.9857, "Na":1.3961, "Mg":1.0081, "Al":0.7307}