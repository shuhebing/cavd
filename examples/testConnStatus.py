# testConnStatus.py
import cavd
conn = cavd.ConnValListCom("Li2CO3-LDA.cif","Li",True,True,None)
oneD,twoD,threeD = ConnStatus(0.5,conn)
print(oneD)
print(twoD)
print(threeD)