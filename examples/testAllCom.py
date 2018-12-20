# testAllCom.py
import cavd
conn,oneD,twoD,threeD,nei_dises,dims,voids = cavd.AllCom("Li2CO3-LDA.cif",0.5,1000,"Li",True,True,None,True,0.5,0.7)
print(conn)
Rf_file = open("Rf_Li.txt","w")
Rf_file.write('a_Rf\tb_Rf\tc_Rf\toneD_Conn\ttwoD_Conn\tthreeD_Conn\tmin_dises\tchannels_dim\tvoids\n')
for i in conn:
    Rf_file.write(str(i)+'\t')
Rf_file.write(str(oneD)+'\t'+str(twoD)+'\t'+str(threeD))
Rf_file.write('\t')
for key in nei_dises:
    Rf_file.write(str(key)+" "+str(nei_dises[key][0])+" "+str(nei_dises[key][1])+" ")
Rf_file.write('\t')
for value in dims:
    Rf_file.write(str(value)+" ")
Rf_file.write('\t')
for void in voids:
    Rf_file.write("("+str(void[0])+","+str(void[1])+","+str(void[2])+") ")
Rf_file.write("\n")