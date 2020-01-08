'''
Used to test channelAnalysis.py channel.

Created on 2019.5.8

@author: YeAnjiang
'''

import cavd

cavd.outChannelToPOSCAR("icsd_16713.cif","Li",ntol=0.02, rad_flag=True, lower=0.0, upper=10.0)

