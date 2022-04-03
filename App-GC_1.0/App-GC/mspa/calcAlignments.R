# -extracted common functionality from mSPA.example.R
# -modified progress reports to cat method name in fromt of step number
# by Nils.Hoffmann@CeBiTec.Uni-Bielefeld.DE (NH)
# on October 5th, 2012
calcAlignments <- function(outputPath, data, parameters) {
	#NH reset gc stats and free any used memory
	gc(reset = TRUE)
    #NH create output directory
    if(!file.exists(outputPath)) {
        output = dir.create(outputPath,recursive=TRUE)
    }
    # NH setting working directory
    #setwd(outputPath)
    
    cwd <- getwd()
    cat(cwd,file='log',append=F) #AS
    nmeth <- length(parameters)
    summary.table <- data.frame()
    stats.table <- data.frame()
    results <- list(summary.table=summary.table,stats.table=stats.table)
    for(i in 1:nmeth) {
        methodParams <- parameters[[i]]
        instanceOutputPath <- paste(outputPath,methodParams$name,i,sep="/")
	cat(instanceOutputPath,file='log',append=T) #AS
        #NH create output directory
        if(!file.exists(instanceOutputPath)) {
            dir.create(instanceOutputPath,recursive=TRUE)
        }
        setwd(instanceOutputPath)
	
        ptm <- proc.time()
        results <- runMethod(outputPath=instanceOutputPath,data=data,parameters=methodParams,results=results,plot=methodParams$plot)
        elapsed <- proc.time() - ptm
        userSysTime <- (elapsed[1]+elapsed[2])*1000.0*1000.0*1000.0
        #record peak memory usages
        memUsed <- gc(verbose=FALSE)
        #add Ncells and Vcells maximum usage
        runTime <- data.frame(field=c(paste("cputime_nanoseconds",userSysTime,sep=" = "),paste("maxUsedMemory_bytes",sum(memUsed[,6])*1024.0*1024.0,sep=" = ")))
        write.table(runTime,file="runtime.txt",row.names=F,col.names=F,sep=" ",quote=F)
    }
    setwd(cwd)
    cat(cwd,file='log',append=T) #AS
}

plotPwAssignments <- function(data) {
    cat("Plotting pairwise peak assignments!\n")
    #nrec <- length(data[[1]])-1
    #trt <- data[[1]]
    #names <- data[[3]]
    #MW: added to show matching of peaks
    #for(j in 1:nrec) {
    #   pdf(file=paste("idsMatchings_", data[[3]][j], ".pdf",sep=""), width=10, height=10)
    #   names<-trt[[j]][,c(4,1,2,3)]
    #   names$name <- as.character(names$name)
    #   plot(names[,2], names[,3], type="n",xlab="RT1 (min)",ylab="RT2 (s)")

        #plot names of compounds
    #   for (i in 1:dim(names)[1]) {
    #           text(names[i,2], names[i,3], names[i,1], cex=0.5)
    #   }

    #   uniqNames <- unique(names$name)
        #plot names of compounds with maximum areas
    #   for (i in 1:length(uniqNames)) {
    #           maxi <- names[which(names$area == max(names[which(names$name == uniqNames[i]),4]) & names$name == uniqNames[i]),]
    #           text(maxi$t1, maxi$t2, maxi$name, col=2, cex=0.5)
    #   }

     #   dev.off()
    #}
}

plotMultipleAssignments <- function(data,plist) {
    cat("Plotting multiple peak assignments!\n")
    pdf("multiplePeakAlignmentPeak.pdf",width=10,height=10)
    rt.table=align(data[[1]],plist,method="peak",plot=T)
    dev.off()
    #pdf("multiplePeakAlignmentArea.pdf",width=10,height=10)
    #rt.table=align(data[[1]],plist,method="area",plot=T)
    #dev.off()
}

# alingment of all peak lists
alignNH <- function(padata,alist,method="peak",plot=FALSE){
    # alist: it should be arranged.
    # e.g.: (1,2), (2,3), (3,4), (4,5), ...
    # Its order sould be consistent with that of padata
    # method: peak, name, area

    len = length(alist)
    wtable = alist[[1]]
    #NH: this method constructs only the multiple alignment against
    # all peaks in the first list
    for(i in 2:len){
        ta = alist[[i]]
        if(colnames(ta)[1]!=colnames(alist[[i-1]])[2]){
            stop("!!! Please make sure the order of alist !!!\n")
        }
        t.pos = match(wtable[,i],ta[,1])
        wtable = data.frame(wtable[which(!is.na(t.pos)),],ta[na.omit(t.pos),2])
        dimnames(wtable)[[2]][i+1] = colnames(ta)[2]
    }
    wlen = dim(wtable)[2]
    #NH: Added compound names list to determine multiple alignment
    #groups
    compoundNames = list()
    #TODO Alignmentgruppen aufbauen
    for(i in 1:wlen){
        compoundNames <- list(compoundNames,padata[[i]]$name)
    }
    compoundNames <- unique(compoundNames)
    if(method=="peak"){
        rtable = list()
        rtable$t1 = matrix(0,dim(wtable)[1],dim(wtable)[2])
        rtable$t2 = rtable$t1
        for(i in 1:wlen){
            tmp.wt = wtable[,i]
            rtable$t1[,i] = padata[[i]]$t1[tmp.wt]
            rtable$t2[,i] = padata[[i]]$t2[tmp.wt]
        }
    }else if(method=="name"){
        rtable = matrix("",dim(wtable)[1],dim(wtable)[2])
        for(i in 1:wlen){
            tmp.wt = wtable[,i]
            rtable[,i] = as.character(padata[[i]]$name[tmp.wt])
        }
    }else if(method=="area"){
        rtable = matrix(0,dim(wtable)[1],dim(wtable)[2])
        for(i in 1:wlen){
            tmp.wt = wtable[,i]
            rtable[,i] = padata[[i]]$area[tmp.wt]
        }
    #NH: added method="index" on Oct. 8th 2012
    }else if(method=="index") {
        rtable = matrix("",dim(wtable)[1],dim(wtable)[2])
        for(i in 1:wlen){
            tmp.wt = wtable[,i]
            #rt <- paste(padata[[i]]$t1[tmp.wt],padata[[i]]$t2[tmp.wt],sep=":")
            #rtable[,i] = paste(tmp.wt,sep="|")
            rtable[,i] = padata[[i]]$ridx[tmp.wt]
        }
    }

    if(plot){
        if(method=="peak"){
            plot(1:10,1:10,type="n",xlim=range(rtable$t1),ylim=range(rtable$t2),xlab="RT1",ylab="RT2")
            for(i in 1:wlen){
                points(rtable$t1[,i],rtable$t2[,i],col=i,pch=i)
            }
        }else if(method=="area"){
            plot(1:10,1:10,type="n",xlim=c(1,dim(rtable)[2]),ylim=range(rtable),xlab="Peak",ylab="Peak Area")
            for(i in 1:dim(rtable)[1]){
                points(1:wlen,rtable[i,],col=i,pch=i,lty=i,type="o")
            }
        }
    }


    list(pos=wtable,info=rtable)
}

runMethod <- function(outputPath,data,parameters,results,plot=FALSE) {
    # NH number of records in dataset
    # data[[1]] contains the raw peak lists
    # data[[2]] contains the condensed peak lists
    # data[[3]] contains the file names
    nrec <- length(data[[1]])-1
    methodName <- parameters$name
    stats.table <- results$stats.table
    summary.table <- results$summary.table

    ####################################
    # peak alignment using "data" method

    #
    # peak alignment for spiked-in samples
    #
    cat(methodName)
    plist = list()
    rlist = list()
    #NH added parameters list and summary data.frame
    for(j in 1:nrec){
        if(plot) {
            pdf(paste("pairwiseAlignment-",data[[3]][1],"-",data[[3]][j+1],".pdf",sep=""))#AP: data[[3]][1]
        }
        #data[[1]] contains the raw peak lists, data[[2]] contains reference data
        #tmp = mspa(padata=data[[1]],rid=j,tid=(j+1),method=parameters$name,w=parameters$w,k=parameters$k,rho=parameters$rho,distance=parameters$distance,similarity=parameters$similarity,anal=parameters$anal,plot=plot)

        #
        #AS: compare all files to file 1 (,which is used as ref)
        tmp = mspa(padata=data[[1]],
        		rid=1,
        		tid=(j+1),
			method=parameters$name,
			w=parameters$w,
			k=parameters$k,
			rho=parameters$rho,
			distance=parameters$distance,
			similarity=parameters$similarity,
			anal=parameters$anal,
			plot=plot)    

        if(plot){
            dev.off()
        }
        cat(".")
        plist[[j]] = tmp$table
        rlist[[j]] = tmp$ridxtable
        colnames(plist[[j]]) <- c(data[[3]][1],data[[3]][j+1])#AP:data[[3]][1]
        colnames(rlist[[j]]) <- c(data[[3]][1],data[[3]][j+1])#AP:data[[3]][1]
        #write.table(data[[1]][plist[[j]],c(j,j+1)],file=paste("pa-",j,"-names.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
        tmpSum <- tmp$summary
        tmpSum$tool <- methodName
        tmpSum$instance <- paste(data[[3]][1],"-",data[[3]][j+1],sep="")#AP:data[[3]][1]
        summary.table <- rbind(summary.table,tmpSum)
    }

    for(i in 1:length(rlist)) {
        cat("Creating pw alignments\n")
        #colnames(data[[2]][[i]])
        #df <- data.frame(x=data[[2]][[i]][plist[[i]],c("ridx")],y=data[[2]][[i+1]][plist[[i+1]],c("ridx")])
        #colnames(df) <- c(data[[3]][i],data[[3]][i+1])
        write.table(rlist[[i]],file=paste("pa-",i,".txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
        #write.table(df,file=paste("pa-",i,".txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
        
        #class(data[[2]][[i]]$ridx[c(unlist(plist[[i]][[data[[3]][i]]]))])
        #str(data[[2]][[i]]$ridx[c(unlist(plist[[i]][[data[[3]][i]]]))])
#         str(length(plist[[i]][[data[[3]][i]]]))
#         str(length(plist[[i]][[data[[3]][i+1]]]))
#         rows <- c()
#         for(j in 1:length(plist[[i]][[data[[3]][i]]])) {
#           if(!is.na(plist[[i]][[data[[3]][i]]][[j]]) && !is.na(plist[[i]][[data[[3]][i+1]]][[j]])) {
#             rows <- c(rows, j)
#           }
#         }
#         str(rows)
#         lhs <- data[[2]][[i]][rows,c("ridx")]
#         rhs <- data[[2]][[i+1]][rows,c("ridx")]
#         df <- data.frame(lhs,rhs)
#         colnames(df) <- c(data[[3]][i],data[[3]][i+1])
#         write.table(df,file=paste("pa-",i,".txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
    }

    # the compound name table
    name.table = malign(data[[1]],plist,method="name")
    col_name = gsub('.cdf_peak','',colnames(name.table$pos))

    name_align_out =name.table$info
    colnames(name_align_out) =col_name    
    write.table(name_align_out,file=paste("name.txt",sep=""),row.names=F,col.names=T,sep=",",quote=F)# the name table
        

    # the retention time table
    #dev.new()
    rt.table = malign(data[[1]],plist,method="peak",plot=F)
    RT_align_out = rt.table$info$t1
    RT2_align_out = rt.table$info$t2
    colnames(RT_align_out) = col_name
    colnames(RT2_align_out) = col_name
    write.table(RT_align_out,file=paste("rt1.csv",sep=""),row.names=F,col.names=T,sep=",",quote=F)
    write.table(RT2_align_out,file=paste("rt2.csv",sep=""),row.names=F,col.names=T,sep=",",quote=F)

    # the peak area table
    area.table=malign(data[[1]],plist,method="area")
    area_align_out =area.table$info
    colnames(area_align_out) =col_name
    write.table(area_align_out,file=paste("area.csv",sep=""),row.names=F,col.names=T,sep=",",quote=F)# the name table


    peak_num_aligned = name.table$info
    peak_num_aligned = cbind(round(rowMeans(rt.table$info$t1, na.rm=TRUE),digits=3),round(rowMeans(rt.table$info$t2, na.rm=TRUE),digits=3), peak_num_aligned)    
        
    colnames(peak_num_aligned) = c('mean_RT','mean_RT2',col_name)
    save(peak_num_aligned,file='save_ali_peak.Rdata')
    write.csv(peak_num_aligned,"peak_num_align.csv", row.names = FALSE, quote = FALSE)

    # NH: the peak index table
    index.table=malign(data[[1]],plist,method="index")
    # NH: set filenames as column names
    colnames(index.table$info) <- data[[3]]
    write.table(index.table$info,file="multiple-alignment.csv",row.names=F,col.names=T,sep="\t",quote=F)# the index table

    # NH: added stats table
    #stats = comp.name(name.table$info)
    #stats.table = rbind(stats.table,cbind(t(stats),NAME=methodName))
    #write.table(stats.table,file="STATS-TABLE.txt",sep="\t",quote=F,row.names=F)
    write.table(summary.table,file="SUMMARY-TABLE.txt",sep="\t",quote=F,row.names=F)

    # NH: added plotting of pairwise peak assignments
    if(plot) {
        plotPwAssignments(data=data)
        #plotMultipleAssignments(data=data,plist=plist)#AP:hide
    }

    # NH: added cat
    cat("done!\n")
    list(stats.table=stats.table,summary.table=summary.table,data=data,plist=plist)
}
