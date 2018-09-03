import cavd
radii = cavd.EffectiveRadCom("./icsd_16713.cif")
print(radii)

cavd.BIComputation(filename="./icsd_16713.cif",migrant="Li",rad_flag=True,effective_rad=True,rad_file=None,rad_store_in_vasp=True,minRad=0.2,maxRad=0.8)
cavd.BIComputation(filename="./Li2CO3-LDA.cif",migrant="Li",rad_flag=False,effective_rad=True,rad_file=None,rad_store_in_vasp=True,minRad=0.2,maxRad=0.8)

Ri,Rf,Rif = cavd.ConnValCom("./icsd_16713.cif","Li",True,True,None)
print(Ri,Rf,Rif)
Ri1,Rf1,Rif1 = cavd.ConnValCom("./Li2CO3-LDA.cif","Li",False,True,None)
print(Ri1,Rf1,Rif1)

conn = cavd.ConnValListCom("./icsd_16713.cif","Li",True,True,None)
conn1 = cavd.ConnValListCom("./Li2CO3-LDA.cif","Li",False,True,None)
print(conn)
print(conn1)

oneD,twoD,threeD = cavd.ConnStatusCom("./icsd_16713.cif",0.4,"Li",True,True,None)
oneD1,twoD1,threeD1 = cavd.ConnStatusCom("./Li2CO3-LDA.cif",0.4,"Li",False,True,None)
print(oneD,twoD,threeD)
print(oneD1,twoD1,threeD1)

oneD2,twoD2,threeD2 = cavd.ConnStatus(0.4,conn)
oneD3,twoD3,threeD3 = cavd.ConnStatus(0.4,conn1)
print(oneD2,twoD2,threeD2)
print(oneD3,twoD3,threeD3)

cavd.ChannelCom("./icsd_16713.cif",0.2,"Li",True,True,None)
cavd.ChannelCom("./Li2CO3-LDA.cif",0.2,"Li",False,True,None)

cavd.ASACom("./icsd_16713.cif",0.5,1000,"Li",True,True,None)
cavd.ASACom("./Li2CO3-LDA.cif",0.5,1000,"Li",False,True,None)

cavd.VoidNetCom("./icsd_16713.cif","Li",True,True,None)
cavd.VoidNetCom("./icsd_16713.cif","Li",False,True,None)
#cavd.AllCom("./icsd_16713.cif",0.5,1000,"Li",True,True,None,True,0.584,0.876)