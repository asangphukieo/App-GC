#NH: Newly created to implement multiple alignment generation with gaps
malign <- function(padata,alist,method="peak",plot=FALSE){
    # alist: it should be arranged.
    # e.g.: (1,2), (2,3), (3,4), (4,5), ...
    # Its order sould be consistent with that of padata
    # method: peak, name, area
    
    len = length(alist)
    
    
    #AP: start by making the first column that contain all matched peaks in file 1
    t.mat = alist[[1]]
    name_ref = colnames(t.mat)[1]
    #print(name_ref)
    for(i in 2:len){
      ta = alist[[i]]
      t.pos = append(t.mat[,1],ta[,1])
      t.uniq= t.pos[!duplicated(t.pos)]
      t.uniq= t.uniq[order(t.uniq)]
      t.mat= data.frame(col1=t.uniq)
    }
    colnames(t.mat)= name_ref
    #print(t.mat)
    #NH: initially two columns (1,2)
    #wtable = t.mat
    #print(wtable[,1])

    wtable = alist[[1]]

    #AP: make first table ; fist column is all peaks in reference file that 
    #match at least 1 peak of other file
    t.temp = match(t.mat[,1],wtable[,1])
    wtable = cbind(t.mat,wtable[t.temp,2])
    colnames(wtable)[2] = colnames(alist[[1]])[2]
    #print(wtable)

    for(i in 2:len){
      ta = alist[[i]]
      
      #AP: match column  of wtable to column 1 
      t.pos = match(wtable[,1],ta[,1])
      #NH: selects rows from wtable that have a match against the lhs file (former rhs file)
      #NH: fuses with columns from ta from rhs file, only keeps peaks present in all files
      #       wtable = data.frame(wtable[which(!is.na(t.pos)),],ta[na.omit(t.pos),2])
      wtable = cbind(wtable,ta[t.pos,2])
      #NH: sets the name of the newly added column to the name of the second column of ta
      dimnames(wtable)[[2]][i+1] = colnames(ta)[2]
    }
    #wtable = wtable[rowSums(is.na(wtable[,1:len]))==0,]
    #NH: create alignment table
    rtable <- malignTable(wtable=wtable,padata=padata,method=method)
    list(pos=wtable,info=rtable)
}
#setwd('/home/corgi/Documents/Nong/msPeak/ABC_file')
#mspa.run(filePrefix='^mut_',filePrefixRef ='^mut_t1_c', plot = TRUE, dataPath = dir.ali.input, directory= dir.ali.main)


malignPeakTable <- function(wtable,padata,method="name") {
  wlen = dim(wtable)[2]
  rtable = matrix(NA,dim(wtable)[1],dim(wtable)[2])
  rtable = list()
  rtable$t1 = matrix(0,dim(wtable)[1],dim(wtable)[2])
  rtable$t2 = rtable$t1
  for(i in 1:wlen){
    tmp.wt = wtable[,i]
    rtable$t1[,i] = padata[[i]]$t1[tmp.wt]
    rtable$t2[,i] = padata[[i]]$t2[tmp.wt]
  }
  rtable
}
    

malignTable <- function(wtable,padata,method="name") {
  if(method=="peak") {
    rtable <- malignPeakTable(wtable,padata)  
    rtable
  }else{
    wlen = dim(wtable)[2]
    rtable = matrix(NA,dim(wtable)[1],dim(wtable)[2])
    for(i in 1:wlen){
      tmp.wt = wtable[,i]
      if(method=="name") {
        rtable[,i] = as.character(padata[[i]]$name[tmp.wt])
      }else if(method=="rt") {
        rtable[,i] = padata[[i]]$t1[tmp.wt]+padata[[i]]$t2[tmp.wt]
      }else if(method=="area") {
        rtable[,i] = padata[[i]]$area[tmp.wt]
      }else if(method=="index") {
        rtable[,i] = padata[[i]]$ridx[tmp.wt]
      }
    }
    rtable
  }
}
