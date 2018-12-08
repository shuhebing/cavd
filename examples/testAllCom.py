# testAllCom.py
import cavd
conn,oneD,twoD,threeD = cavd.AllCom("Li2CO3-LDA.cif",0.5,1000,"Li",True,True,None,True,0.5,0.7)
print(conn)
print(oneD)
print(twoD)
print(threeD)