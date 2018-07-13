from zeo.netinfo import readIonRadTable
from zeo.netinfo import lookupIonRadius
a={'C4+':0.3, 'O2-':1.23, 'Li+':1.06}
readIonRadTable(a)
b=lookupIonRadius('C4+')
print(b)
