# testConnStatusCom.py
import cavd
oneD,twoD,threeD = cavd.ConnStatusCom("Li2CO3-LDA.cif",0.5,"Li",True,True,None)
print(oneD)
print(twoD)
print(threeD)