# load mSPA package
source("mSPA.R")

# load the ChromaToF data
pad1 = load.data(pattern="^Standard_*")
pad2 = load.data(pattern="^Deblank_*")

####################################
# peak alignment using "PAD" method

#
# peak alignment for Standard compounds
#
plist1 = list()
for(j in 1:9){
	tmp = mspa(padata=pad1[[2]],rid=j,tid=(j+1),method="PAD")
	cat("#1# ",j,"th Done!\n")
	plist1[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table1=align(pad1[[2]],plist1,method="peak",plot=F)
write.table(rt.table1$info$t1,file="PAD1_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table1$info$t2,file="PAD1_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table1=align(pad1[[2]],plist1,method="area")
write.table(area.table1$info,file="PAD1_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table1=align(pad1[[2]],plist1,method="name")
write.table(name.table1$info,file="PAD1_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

#
# peak alignment for spiked-in samples
#
plist2 = list()
for(j in 1:4){
	tmp = mspa(padata=pad2[[2]],rid=j,tid=(j+1),method="PAD")
	cat("#2# ",j,"th Done!\n")
	plist2[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table2=align(pad2[[2]],plist2,method="peak",plot=F)
write.table(rt.table2$info$t1,file="PAD2_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table2$info$t2,file="PAD2_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table2=align(pad2[[2]],plist2,method="area")
write.table(area.table2$info,file="PAD2_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table2=align(pad2[[2]],plist2,method="name")
write.table(name.table2$info,file="PAD2_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

####################################
# peak alignment using "PAS" method

#
# peak alignment for Standard compounds
#
plist1 = list()
for(j in 1:9){
	tmp = mspa(padata=pad1[[2]],rid=j,tid=(j+1),method="PAS")
	cat("#1# ",j,"th Done!\n")
	plist1[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table1=align(pad1[[2]],plist1,method="peak",plot=F)
write.table(rt.table1$info$t1,file="PAS1_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table1$info$t2,file="PAS1_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table1=align(pad1[[2]],plist1,method="area")
write.table(area.table1$info,file="PAS1_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table1=align(pad1[[2]],plist1,method="name")
write.table(name.table1$info,file="PAS1_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

#
# peak alignment for spiked-in samples
#
plist2 = list()
for(j in 1:4){
	tmp = mspa(padata=pad2[[2]],rid=j,tid=(j+1),method="PAS")
	cat("#2# ",j,"th Done!\n")
	plist2[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table2=align(pad2[[2]],plist2,method="peak",plot=F)
write.table(rt.table2$info$t1,file="PAS2_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table2$info$t2,file="PAS2_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table2=align(pad2[[2]],plist2,method="area")
write.table(area.table2$info,file="PAS2_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table2=align(pad2[[2]],plist2,method="name")
write.table(name.table2$info,file="PAS2_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

####################################
# peak alignment using "DW-PAS" method

#
# peak alignment for Standard compounds
#
plist1 = list()
for(j in 1:9){
	tmp = mspa(padata=pad1[[2]],rid=j,tid=(j+1),k=3,distance="maximum",method="DW-PAS")
	cat("#1# ",j,"th Done!\n")
	plist1[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table1=align(pad1[[2]],plist1,method="peak",plot=F)
write.table(rt.table1$info$t1,file="DW-PAS1_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table1$info$t2,file="DW-PAS1_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table1=align(pad1[[2]],plist1,method="area")
write.table(area.table1$info,file="DW-PAS1_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table1=align(pad1[[2]],plist1,method="name")
write.table(name.table1$info,file="DW-PAS1_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

#
# peak alignment for spiked-in samples
#
plist2 = list()
for(j in 1:4){
	tmp = mspa(padata=pad2[[2]],rid=j,tid=(j+1),k=15,method="DW-PAS")
	cat("#2# ",j,"th Done!\n")
	plist2[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table2=align(pad2[[2]],plist2,method="peak",plot=F)
write.table(rt.table2$info$t1,file="DW-PAS2_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table2$info$t2,file="DW-PAS2_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table2=align(pad2[[2]],plist2,method="area")
write.table(area.table2$info,file="DW-PAS2_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table2=align(pad2[[2]],plist2,method="name")
write.table(name.table2$info,file="DW-PAS2_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

####################################
# peak alignment using "SW-PAD" method

#
# peak alignment for Standard compounds
#
plist1 = list()
for(j in 1:9){
	tmp = mspa(padata=pad1[[2]],rid=j,tid=(j+1),rho=0.5,method="SW-PAD")
	cat("#1# ",j,"th Done!\n")
	plist1[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table1=align(pad1[[2]],plist1,method="peak",plot=F)
write.table(rt.table1$info$t1,file="SW-PAD1_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table1$info$t2,file="SW-PAD1_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table1=align(pad1[[2]],plist1,method="area")
write.table(area.table1$info,file="SW-PAD1_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table1=align(pad1[[2]],plist1,method="name")
write.table(name.table1$info,file="SW-PAD1_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

#
# peak alignment for spiked-in samples
#
plist2 = list()
for(j in 1:4){
	tmp = mspa(padata=pad2[[2]],rid=j,tid=(j+1),rho=0.93,distance="maximum",method="SW-PAD")
	cat("#2# ",j,"th Done!\n")
	plist2[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table2=align(pad2[[2]],plist2,method="peak",plot=F)
write.table(rt.table2$info$t1,file="SW-PAD2_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table2$info$t2,file="SW-PAD2_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table2=align(pad2[[2]],plist2,method="area")
write.table(area.table2$info,file="SW-PAD2_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table2=align(pad2[[2]],plist2,method="name")
write.table(name.table2$info,file="SW-PAD2_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

####################################
# peak alignment using "PAM" method

#
# peak alignment for Standard compounds
#
plist1 = list()
for(j in 1:9){
	tmp = mspa(padata=pad1[[2]],rid=j,tid=(j+1),w=0.5,method="PAM")
	cat("#1# ",j,"th Done!\n")
	plist1[[j]] = tmp$table
}
# the retention time table
#dev.new()
rt.table1=align(pad1[[2]],plist1,method="peak",plot=F)
write.table(rt.table1$info$t1,file="PAM1_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table1$info$t2,file="PAM1_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table1=align(pad1[[2]],plist1,method="area")
write.table(area.table1$info,file="PAM1_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table1=align(pad1[[2]],plist1,method="name")
write.table(name.table1$info,file="PAM1_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

#
# peak alignment for spiked-in samples
#
plist2 = list()
for(j in 1:4){
	tmp = mspa(padata=pad2[[2]],rid=j,tid=(j+1),w=0.05,distance="maximum",method="PAM")
	cat("#2# ",j,"th Done!\n")
	plist2[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table2=align(pad2[[2]],plist2,method="peak",plot=F)
write.table(rt.table2$info$t1,file="PAM2_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table2$info$t2,file="PAM2_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table2=align(pad2[[2]],plist2,method="area")
write.table(area.table2$info,file="PAM2_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table2=align(pad2[[2]],plist2,method="name")
write.table(name.table2$info,file="PAM2_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

####################################
# peak alignment using "OP-PAM" method

#
# peak alignment for Standard compounds
#
plist1 = list()
for(j in 1:9){
	tmp = mspa(padata=pad1[[2]],rid=j,tid=(j+1),method="OP-PAM")
	cat("#1# ",j,"th Done!\n")
	plist1[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table1=align(pad1[[2]],plist1,method="peak",plot=F)
write.table(rt.table1$info$t1,file="OP-PAM1_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table1$info$t2,file="OP-PAM1_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table1=align(pad1[[2]],plist1,method="area")
write.table(area.table1$info,file="OP-PAM1_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table1=align(pad1[[2]],plist1,method="name")
write.table(name.table1$info,file="OP-PAM1_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

#
# peak alignment for spiked-in samples
#
plist2 = list()
for(j in 1:4){
	tmp = mspa(padata=pad2[[2]],rid=j,tid=(j+1),method="OP-PAM")
	cat("#2# ",j,"th Done!\n")
	plist2[[j]] = tmp$table
}

# the retention time table
#dev.new()
rt.table2=align(pad2[[2]],plist2,method="peak",plot=F)
write.table(rt.table2$info$t1,file="OP-PAM2_rt1.txt",row.names=F,col.names=F,sep=" ",quote=F)
write.table(rt.table2$info$t2,file="OP-PAM2_rt2.txt",row.names=F,col.names=F,sep=" ",quote=F)

# the peak area table
area.table2=align(pad2[[2]],plist2,method="area")
write.table(area.table2$info,file="OP-PAM2_area.txt",row.names=F,col.names=F,sep=" ",quote=F)# the name table

# the compound name table
name.table2=align(pad2[[2]],plist2,method="name")
write.table(name.table2$info,file="OP-PAM2_name.txt",row.names=F,col.names=F,sep=";",quote=F)# the name table

