# testConnVal.py
import cavd
Ri,Rf,Rif = cavd.ConnValCom("Li2CO3-LDA.cif","Li",True,True,None)
print(Ri)
print(Rf)
print(Rif)
