################################
# LOAD the ChromaToF DATA    #
# 1: MIXTURE (STANDARD) (#=10) #
# 2: SPIKE-IN (DEBLANK) (#=5)  #
################################

#Modifications done on June 22nd, 2012
#Nils.Hoffmann@CeBiTec.Uni-Bielefeld.DE
#-added path parameter to rt.import and load.data
#-added full.names=TRUE to allow relative paths
#The original sources were otherwise unaltered.

# import the data
rt.import <- function(path=".",pattern="^STANDARD*"){
    cat("Pattern: ",pattern,"\n")
    fname = dir(path=path,pattern=pattern,full.names=TRUE)
    cat("Files: ","\n")
    print(fname)
    rt = NULL
    lapply(1:length(fname),function(i) read.table(fname[i],sep='\t',header=TRUE))
}


# import the data
#AP: modified for fixing reference file
#hf=rt.import.setref(filePrefixRef=filePrefixRef,path='.',pattern=pattern)
rt.import.setref <- function(path=".",pattern="^STANDARD*",filePrefixRef){
    cat("Pattern: ",pattern,"\n")
    fname = dir(path=path,pattern=pattern,full.names=TRUE) #reorder this list , to set reference file
    cat("Files: ","\n")
    print(fname)

    rt = NULL
    rt=lapply(1:length(fname),function(i) read.table(fname[i],sep='\t',header=TRUE))

    if(filePrefixRef!=F){
        
        cat("Pattern reference file: ",filePrefixRef,"\n")
        fnameref = dir(path=path,pattern=filePrefixRef,full.names=TRUE) #reorder this list , to set reference file
        cat("Reference File for peak alignment: ","\n")
        ref_file = fnameref[1]#use only th first file as reference
        print(ref_file) 

        if (ref_file %in% fname){
             cat("Reference file was found in the input list...\n") 
             fname = fname[fname != ref_file]
             fname = append(ref_file,fname )

             rt = NULL
             rt=lapply(1:length(fname),function(i) read.table(fname[i],sep='\t',header=TRUE))

        }else{
             cat("Reference file was not found in the input list...\n") 
             fname = append(ref_file,fname )

             rtref = NULL
             rtref= read.table(ref_file,sep='\t',header=TRUE)             
             rt = append(list(rtref), rt)

        }
    }

    list("rt"=rt,"fname"=fname)

}



# extract the spec info
spec <- function(data=rt){
    len = length(data)
    nd = NULL
    nd = lapply(1:len,
                function(i){
                    sp=data[[i]]$spectra;
                    strsplit(as.character(sp),split=" ");
                }
            )
    sd = NULL
    sd = lapply(1:len,
                function(i){
                    td=nd[[i]];
                    len2=length(td);
                    lapply(1:len2,
                            function(j){
                                len3=length(td[[j]]);
                                td2=t(matrix(as.numeric(unlist(strsplit(td[[j]],split=":"))),2,len3));
                                td2[order(td2[,1]),];
                            }
                        )
                }
            )
    sd
}

#NH: plot of group stats
plotGroupStats <- function(outputPath=".",peaks,compoundGroupStats,compoundNames,sd1Thres=25.0,sd2Thres=0.05,maxX,maxY) {
    cat("Creating compound group statistics plots!\n")
    #NH: store working dir
    oldwd <- getwd()
    #NH: Plot showing decisionBoundary for stdev in t1 and t2
    library("ggplot2")
    library("RColorBrewer")
    #NH create output directory
    if(!file.exists(outputPath)) {
        output = dir.create(outputPath,recursive=TRUE)
    }
    #NH: set working dir
    setwd(outputPath)
    #NH: select items with non-zero stddev in either t1 or t2
    nitems <- length(unique(peaks$name))

    if(nitems!=length(unique(compoundNames))) {
        stop("Compound number mismatch!")
    }
    shapes <- rep(x=0:14,length=as.numeric(nitems))
    #nrow(compoundGroupStats[compoundGroupStats$t1Stdev>0 | compoundGroupStats$t2Stdev>0,])
    itemHeight <- 3
    itemWidth <- 7
    cols <- 3
    rows <- ceiling(nitems/cols)
    rowsPerPage <- min(rows,10)
    pages <- ceiling(rows/rowsPerPage)
    itemsPerPage <- rowsPerPage*cols
    #cat("Items:",nitems," cols: ",cols," rows: ",rows," rows per page: ",rowsPerPage," pages: ",pages," items per page: ",itemsPerPage,"\n")
    decisionBoundaryParameters <- data.frame(sd1Thres=sd1Thres,sd2Thres=sd2Thres)
    write.csv(decisionBoundaryParameters,file="decisionBoundaryParameters.txt")
    decisionBoundary <- data.frame()
    for(l in 1:100) {
        x <- (l-1)*(sd1Thres/100.0)
        y <- ellipse2d(x,a=sd1Thres,b=sd2Thres)
        decisionBoundary <- rbind(decisionBoundary, data.frame(x=x,y=y))
    }
    write.csv(decisionBoundary,file="decisionBoundaryAll.txt")
    write.csv(compoundGroupStats,file="compoundGroupStatsAll.txt")
    pdf("all-groups-rt1-rt2.pdf",width=20,height=12)
    #print(shapes)
    legend.position <- "bottom"
    if(length(compoundNames)>100) {
      legend.position <- "none"
    }
    glegend <- guide_legend(ncol=5,"Compound",title.position="top")
    print(
        ggplot(data=peaks,aes(x=t1,y=t2,colour=name,shape=name))+geom_point(alpha=0.5)+
        xlab("RT1")+ylab("RT2")+
        scale_shape_manual(values=shapes,limits=compoundNames,drop=FALSE,labels=compoundNames)+
        scale_color_manual(values = rep(c(brewer.pal(5, "Set1"),brewer.pal(9,"Set1")[7:9]),length=length(compoundNames)),drop=FALSE,limits=compoundNames,labels=compoundNames)+
        guides(shape = glegend,colour = glegend)+
        xlim(0,maxX)+ylim(0,maxY)+
        theme(legend.position=legend.position)
    )
    dev.off()
    offset <- 1
    for(i in 1:pages) {
        startItem <- offset
        itemsLeft <- nitems-startItem+1
        #cat("Items left: ",itemsLeft,"\n")
        rowsOnPage <- rowsPerPage
        colsOnPage <- cols
        if(itemsLeft <= itemsPerPage) {
            rowsOnPage <- ceiling(itemsLeft/cols)
            colsOnPage <- min(cols,itemsLeft)
        }
        #cat("Creating ",rowsOnPage," rows and ",colsOnPage," cols\n")
        pdf(paste("stddevs-rt1-rt2-all-",i,".pdf",sep=""),width=colsOnPage*itemWidth,height=itemHeight*rowsOnPage)
        selection <- seq(from=startItem,to=min(nitems,startItem+itemsPerPage-1),by=1)
        #cat("Items on page: ",length(selection),"\n")
        #cat("Creating page ",i," of ",pages," from ",startItem," to ",startItem+length(selection),"\n")
        print(
            ggplot(data=compoundGroupStats[selection,],
                aes(x=t1Stdev,y=t2Stdev,colour=suspicious,shape=suspicious,size=size))+
            scale_colour_manual(name="suspicious",values=c("FALSE"=rgb(5, 113, 176, max = 255),"TRUE"=rgb(202, 0, 32, max = 255)))+
            scale_size_area(name="Size")+
            geom_point(alpha=0.5)+
            annotate("ribbon",x=decisionBoundary$x,ymin=rep(0.0,length(decisionBoundary$y)),ymax=decisionBoundary$y,alpha=.2,fill="grey50")+
            xlab("sd(RT1)")+
            ylab("sd(RT2)")+
            scale_y_log10()+
            scale_x_log10()+
            facet_wrap(~name,ncol=colsOnPage,nrow=rowsOnPage,drop=TRUE)
        )
        dev.off()
        pdf(paste("rt1-rt2-all-",i,".pdf",sep=""),width=colsOnPage*itemWidth,height=itemHeight*rowsOnPage)
        #NH: Find matching positions in peaks for name name
        upos <- match(as.character(peaks$name),as.character(compoundGroupStats[selection,]$name))
        upos <- which(!is.na(upos))
        print(
            ggplot(data=peaks[upos,],aes(x=t1,y=t2,colour=filename,shape=suspicious))+
            scale_colour_hue()+
            geom_point(alpha=0.5,size=3)+
            xlab("RT1")+
            ylab("RT2")+
            xlim(0,maxX)+ylim(0,maxY)+
            facet_wrap(~name,ncol=colsOnPage,nrow=rowsOnPage,drop=TRUE)
        )
        dev.off()
        pdf(paste("t1-all-bpl-",i,".pdf",sep=""),width=colsOnPage*itemWidth,height=itemHeight*rowsOnPage)
        print(
            qplot(suspicious,t1,data=peaks[upos,],geom="boxplot",colour=filename)+coord_flip()+
            facet_wrap(~name,ncol=colsOnPage,nrow=rowsOnPage,drop=TRUE)+xlab("Outlier")+ylab("RT1")+ylim(0,maxX)
        )
        dev.off()
        pdf(paste("t2-all-bpl-",i,".pdf",sep=""),width=colsOnPage*itemWidth,height=itemHeight*rowsOnPage)
        print(
            qplot(suspicious,t2,data=peaks[upos,],geom="boxplot",colour=filename)+coord_flip()+
            facet_wrap(~name,ncol=colsOnPage,nrow=rowsOnPage,drop=TRUE)+xlab("Outlier")+ylab("RT2")+ylim(0,maxY)
        )
        dev.off()
        offset <- offset + length(selection)
    }
    #NH: reset working dir
    setwd(oldwd)
    #cat("Done plotting compound group stats!\n")
}

createCompoundIds <- function(urt,fnames) {
    ids <- data.frame()
    for(i in 1:length(urt)) {
        report <- urt[[i]]
        filename <- fnames[[i]]
        #for(j in 1:nrow(report)) {
        ids <- rbind(ids,data.frame(id=createCompoundId(filename,report$ridx,report$name),t1=report$t1,t2=report$t2,ridx=report$ridx, name=report$name, area=report$area, filename=filename))
        #}
    }
    #print(ids)
    ids
}

createCompoundId <- function(filename, ridx, name) {
  paste(filename,ridx,name,sep="-")
}

load.data.raw <- function(path=".",pattern) {
    rt = rt.import(path=path,pattern=pattern) # real data rt[[1]][,]
    nrt = spec(rt) # spectrum nrt[[1]][[1]][,]
    #NH: added from rt.import to keep filenames with dataPath, remove .csv ending
    fnames <- dir(path=path,pattern=pattern,full.names=FALSE)
    fnames <- sub(".csv","",fnames)

    #NH: added for plotting area
    maxX <- 0
    maxY <- 0

    # 1   2   3     4     5      6-
    # t1, t2, area, name, simil, spect...,

    #NH: a list of all compound names from all peak lists
    compoundNames <- list()
    rows <- c()
    trt = NULL
    for(i in 1:length(nrt)){
        tspec = NULL
        for(j in 1:length(nrt[[i]])){
            tspec = rbind(tspec,as.numeric(nrt[[i]][[j]][,2]))
        }
        trt[[i]] = data.frame(
                #NH: added to have original index
                ridx = seq(from=0,to=length(nrt[[i]])-1,by=1)
                ,t1 = (rt[[i]]$RT)
                ,t2 = (rt[[i]]$RT2)
                ,area = rt[[i]]$area
                #NH: modified to remove leading and trailing whitespace
                #,name = rt[[i]]$Name
                ,name = sub('^ +','',sub(' +$', '', rt[[i]]$peak_ID))
                #,simil = rt[[i]]$Simil
                ,tspec # spectrum # 6---
            )
    }
    #NH: remove duplicate peaks with equal area, equal name
    cat("Checking for duplicates\n")
    for(i in 1:length(trt)){
        uniques <- unique(trt[[i]]$area)
        total <- trt[[i]]$area
        matches <- match(uniques,total)
        if(length(uniques)!=length(total)) {
          delta <- setdiff(seq(1:nrow(trt[[i]])),matches)
          cat("Matches at ",delta,"\n")
          print(trt[[i]][delta,c("name","t1","t2","area")])
          cat("Removed ",length(total)-length(uniques)," peaks from file ",fnames[[i]],"\n")
        }
        trt[[i]] = trt[[i]][matches,]
        #NH: compoundNames
        compoundNames <- c(compoundNames,as.list(unique(as.character(trt[[i]]$name))))
        #rows <- c(rows,nrow(trt[[i]]))
        maxX <- max(maxX,trt[[i]]$t1)
        maxY <- max(maxY,trt[[i]]$t2)
    }
    #NH: Keep only unique names and sort
    compoundNames <- as.list(sort(unlist(unique(compoundNames))))
    list(rawPeakLists=trt,fileNames=fnames,compoundNames=compoundNames,maxX=maxX,maxY=maxY)
}

#AP:modified for setting single reference file
load.data.raw.setref <- function(path=".",pattern,filePrefixRef) {
    rt = rt.import.setref(path=path,pattern=pattern,filePrefixRef=filePrefixRef) # real data rt$rt[[1]][,]
    nrt = spec(rt$rt) # spectrum nrt[[1]][[1]][,]
    #NH: added from rt.import to keep filenames with dataPath, remove .csv ending

    #AP:
    fnames <- rt$fname
    fnames <- sub(".csv","",fnames)
    fnames <- sapply(1:length(fnames),function(i) tail(strsplit(fnames[i],'/')[[1]],n=1))

    #NH: added for plotting area
    maxX <- 0
    maxY <- 0

    # 1   2   3     4     5      6-
    # t1, t2, area, name, simil, spect...,

    #NH: a list of all compound names from all peak lists
    compoundNames <- list()
    rows <- c()
    trt = NULL
    for(i in 1:length(nrt)){
        tspec = NULL
        for(j in 1:length(nrt[[i]])){
            tspec = rbind(tspec,as.numeric(nrt[[i]][[j]][,2]))
        }
        trt[[i]] = data.frame(
                #NH: added to have original index
                ridx = seq(from=0,to=length(nrt[[i]])-1,by=1)
                ,t1 = (rt$rt[[i]]$RT)
                ,t2 = (rt$rt[[i]]$RT2)
                ,area = rt$rt[[i]]$area
                #NH: modified to remove leading and trailing whitespace
                #,name = rt[[i]]$Name
                ,name = sub('^ +','',sub(' +$', '', rt$rt[[i]]$peak_ID))
                #,simil = rt[[i]]$Simil
                ,tspec # spectrum # 6---
            )
    }
    #NH: remove duplicate peaks with equal area, equal name
    cat("Checking for duplicates by area (AP: ignore)\n")
    for(i in 1:length(trt)){
        uniques <- unique(trt[[i]]$name) #AP: Ignore: removed area duplication (change $area to $name)
        total <- trt[[i]]$name #AP: Ignore: removed area duplication (change $area to $name)
        matches <- match(uniques,total)
        if(length(uniques)!=length(total)) {
          delta <- setdiff(seq(1:nrow(trt[[i]])),matches)
          cat("Matches at ",delta,"\n")
          print(trt[[i]][delta,c("name","t1","t2","area")])
          cat("Removed ",length(total)-length(uniques)," peaks from file ",fnames[[i]],"\n")
        }
        trt[[i]] = trt[[i]][matches,]
        #NH: compoundNames
        compoundNames <- c(compoundNames,as.list(unique(as.character(trt[[i]]$name))))
        #rows <- c(rows,nrow(trt[[i]]))
        maxX <- max(maxX,trt[[i]]$t1)
        maxY <- max(maxY,trt[[i]]$t2)
    }
    #NH: Keep only unique names and sort
    compoundNames <- as.list(sort(unlist(unique(compoundNames))))
    list(rawPeakLists=trt,fileNames=fnames,compoundNames=compoundNames,maxX=maxX,maxY=maxY)
}

createPeaks <- function(trt,fnames) {
    #NH: a dataframe containing all peaks for statistics generation
    peaks <- data.frame()
    #cat("Creating peaks data.frame!\n")
    for(i in 1:length(trt)){
        #NH: For each data.frame/chromatof report
        utmp = trt[[i]]
        #cat("Processing report ",i,"\n")
        #NH: Add ridx, t1, t2, area and name to the peaks table
        peaks <- rbind(peaks,data.frame(ridx=utmp$ridx,t1=utmp$t1,t2=utmp$t2,area=utmp$area,name=utmp$name,filename=fnames[[i]],suspicious=FALSE))
    }
    cat("Loaded ",nrow(peaks)," peaks in total\n")
    peaks
}
# load the ChromaToF data
load.data <- function(path=".",pattern,plot=FALSE,sd1Thres=25.0,sd2Thres=0.05,outputPath=".",createReferenceAlignment=FALSE,variant="mSPA"){
    rawData <- load.data.raw(path=path,pattern=pattern)
    trt <- rawData$rawPeakLists
    fnames <- rawData$fileNames
    compoundNames <- rawData$compoundNames
    peaks <- createPeaks(trt,fnames)
    compoundGroupStats <- createCompoundGroupStats(peaks,compoundNames,sd1Thres,sd2Thres)
    peaks <- compoundGroupStats$peaks
    if(plot) {
        plotGroupStats(outputPath,peaks,compoundGroupStats$stats,compoundNames,sd1Thres,sd2Thres,rawData$maxX,rawData$maxY)
    }
    cat("Creating compound ids\n")
    # multiple peak merging
    results <- resolveAmbiguousPeaks(trt,fnames,compoundNames,compoundGroupStats$stats,variant=variant,sd1Thres,sd2Thres)
    urt <- results[[1]]
    if(createReferenceAlignment) {
        createReferenceAlignment(compoundNames,urt,outputPath,fnames,plot,rawData$maxX,rawData$maxY)
        cat("Creating custom reference ids\n")
        rawDataMgma <- load.data.raw(path=path,pattern=pattern)
        idsAll <- createCompoundIds(rawDataMgma$rawPeakLists,fnames)
        mgmaResults <- resolveAmbiguousPeaks(rawDataMgma$rawPeakLists,fnames,compoundNames,compoundGroupStats$stats,variant="mgma",sd1Thres,sd2Thres)
        createReferenceAlignment(compoundNames,mgmaResults[[1]],outputPath,fnames,plot,rawData$maxX,rawData$maxY,nameSuffix="mgma")
        urtCustom <- mgmaResults[[1]]
        ids1 <- createCompoundIds(urtCustom,fnames)
        cat("Creating ",variant," reference ids\n")
        ids2 <- createCompoundIds(urt,fnames)
        cat("Creating unions, intersections and set differences\n")
        u <- union(ids1$id,ids2$id)
        #idsAll[idsAll$id==u,] <-
        isec <- intersect(ids1$id,ids2$id)
        sdiff1 <- setdiff(ids1$id,ids2$id)
        sdiff2 <- setdiff(ids2$id,ids1$id)
        #NH: create pairwise intersections
        #all peaks
        A <- idsAll
        #our custom selected peaks
        B <- ids1
        #mSPA selected peaks
        C <- ids2

        #intersections
        #all with mgma
        nab <- intersect(A$id, B$id)
        #mgma and mspa
        nbc <- intersect(B$id, C$id)
        #all and mspa
        nac <- intersect(A$id, C$id)
        #intersection of all sets
        nabc <- intersect(nab, C$id)
        #a intersect b without c
        nabwoc <- setdiff(nab, C$id)
        #a intersect c without b
        nacwob <- setdiff(nac, B$id)

        #exclusive sets
        awobc <- setdiff(A$id, union(B$id,C$id))
        bwoac <- setdiff(B$id, C$id)
        cwoab <- setdiff(C$id, B$id)

        A$set <- "Unassigned"
        common <- which(!is.na(match(A$id,nabc)))
        if(length(common)>0) {
            A[common,]$set <- "Common"
        }
        mmdgExclusive <- which(!is.na(match(A$id,bwoac)))
        if(length(mmdgExclusive)>0) {
            A[mmdgExclusive,]$set <- "mgma Exclusive"
        }
        mspaExclusive <- which(!is.na(match(A$id,cwoab)))
        if(length(mspaExclusive)>0) {
            A[mspaExclusive,]$set <- "mSPA Exclusive"
        }
        A$set <- as.factor(A$set)
        #NH: store working dir
        oldwd <- getwd()
        #NH create output directory
        if(!file.exists(outputPath)) {
            output = dir.create(outputPath,recursive=TRUE)
        }
        #NH: set working dir
        setwd(outputPath)
        write.csv(A,"peak-set-partition.csv")
        cat("Found ",nrow(idsAll)," raw peaks\n")
        cat("Found ",length(common)," common peaks\n")
        cat("Unique to mgma: ",length(bwoac),"\n")
        cat("Unique to mSPA: ",length(cwoab),"\n")
        if(plot) {
            library("RColorBrewer")
            pdf("peak-set-partition.pdf",width=16,height=14)
            pal <- brewer.pal(5, "Set1")
            shapes <- rep(x=0:14,length=as.numeric(length(compoundNames)))
            legend.position <- "bottom"
            if(length(compoundNames)>100) {
              legend.position <- "none"
            }
            glegend <- guide_legend(ncol=5,title.position="top",title="Compound")
            print(
                ggplot(data=A,aes(x=t1,y=t2,colour=name,shape=name))+
                geom_point(alpha=0.5,size=3)+
                scale_shape_manual(values=shapes,limits=compoundNames,drop=FALSE,labels=compoundNames)+
                scale_color_manual(values = rep(c(brewer.pal(5, "Set1"),brewer.pal(9,"Set1")[7:9]),length=length(compoundNames)),drop=FALSE,limits=compoundNames,labels=compoundNames)+
                guides(shape = glegend,colour = glegend)+
                #scale_colour_manual(values=c("Unassigned"=pal[2],"Mmdg Exclusive"=pal[3],"mSPA Exclusive"=pal[4],"Common"=rgb(112,153,117,maxColorValue=255)))+
                xlab("RT1")+
                ylab("RT2")+
                xlim(0,rawData$maxX)+ylim(0,rawData$maxY)+
                facet_grid(~set)+
                theme(legend.position=legend.position)
            )
            dev.off()
        }
        setwd(oldwd)
    }
    #NH: modified return value to contain filenames and compoundNames
    list(trt,urt,fnames,compoundNames)
}

# load the ChromaToF data
#AP:modified for setting single reference file
load.data.setref <- function(path=".",pattern,filePrefixRef,plot=FALSE,sd1Thres=25.0,sd2Thres=0.05,outputPath=".",createReferenceAlignment=FALSE, variant="mSPA"){
    rawData <- load.data.raw.setref(path=path,pattern=pattern,filePrefixRef=filePrefixRef)
    trt <- rawData$rawPeakLists
    fnames <- rawData$fileNames
    compoundNames <- rawData$compoundNames
    peaks <- createPeaks(trt,fnames)
    compoundGroupStats <- createCompoundGroupStats(peaks,compoundNames,sd1Thres,sd2Thres)
    peaks <- compoundGroupStats$peaks
    #AS//
    #if(plot) {
    #    plotGroupStats(outputPath,peaks,compoundGroupStats$stats,compoundNames,sd1Thres,sd2Thres,rawData$maxX,rawData$maxY)
    #}
    #AS//
    cat("Creating compound ids\n")
    # multiple peak merging
    results <- resolveAmbiguousPeaks(trt,fnames,compoundNames,compoundGroupStats$stats,variant=variant,sd1Thres,sd2Thres)
    urt <- results[[1]]
    if(createReferenceAlignment) {
        createReferenceAlignment(compoundNames,urt,outputPath,fnames,plot,rawData$maxX,rawData$maxY)
        cat("Creating custom reference ids\n")
        rawDataMgma <- load.data.raw.setref(path=path,pattern=pattern,filePrefixRef=filePrefixRef)
        idsAll <- createCompoundIds(rawDataMgma$rawPeakLists,fnames)
        mgmaResults <- resolveAmbiguousPeaks(rawDataMgma$rawPeakLists,fnames,compoundNames,compoundGroupStats$stats,variant="mgma",sd1Thres,sd2Thres)
        createReferenceAlignment(compoundNames,mgmaResults[[1]],outputPath,fnames,plot,rawData$maxX,rawData$maxY,nameSuffix="mgma")
        urtCustom <- mgmaResults[[1]]
        ids1 <- createCompoundIds(urtCustom,fnames)
        cat("Creating ",variant," reference ids\n")
        ids2 <- createCompoundIds(urt,fnames)
        cat("Creating unions, intersections and set differences\n")
        u <- union(ids1$id,ids2$id)
        #idsAll[idsAll$id==u,] <-
        isec <- intersect(ids1$id,ids2$id)
        sdiff1 <- setdiff(ids1$id,ids2$id)
        sdiff2 <- setdiff(ids2$id,ids1$id)
        #NH: create pairwise intersections
        #all peaks
        A <- idsAll
        #our custom selected peaks
        B <- ids1
        #mSPA selected peaks
        C <- ids2

        #intersections
        #all with mgma
        nab <- intersect(A$id, B$id)
        #mgma and mspa
        nbc <- intersect(B$id, C$id)
        #all and mspa
        nac <- intersect(A$id, C$id)
        #intersection of all sets
        nabc <- intersect(nab, C$id)
        #a intersect b without c
        nabwoc <- setdiff(nab, C$id)
        #a intersect c without b
        nacwob <- setdiff(nac, B$id)

        #exclusive sets
        awobc <- setdiff(A$id, union(B$id,C$id))
        bwoac <- setdiff(B$id, C$id)
        cwoab <- setdiff(C$id, B$id)

        A$set <- "Unassigned"
        common <- which(!is.na(match(A$id,nabc)))
        if(length(common)>0) {
            A[common,]$set <- "Common"
        }
        mmdgExclusive <- which(!is.na(match(A$id,bwoac)))
        if(length(mmdgExclusive)>0) {
            A[mmdgExclusive,]$set <- "mgma Exclusive"
        }
        mspaExclusive <- which(!is.na(match(A$id,cwoab)))
        if(length(mspaExclusive)>0) {
            A[mspaExclusive,]$set <- "mSPA Exclusive"
        }
        A$set <- as.factor(A$set)
        #NH: store working dir
        oldwd <- getwd()
        #NH create output directory
        if(!file.exists(outputPath)) {
            output = dir.create(outputPath,recursive=TRUE)
        }
        #NH: set working dir
        setwd(outputPath)
        write.csv(A,"peak-set-partition.csv")
        cat("Found ",nrow(idsAll)," raw peaks\n")
        cat("Found ",length(common)," common peaks\n")
        cat("Unique to mgma: ",length(bwoac),"\n")
        cat("Unique to mSPA: ",length(cwoab),"\n")
        if(plot) {
            library("RColorBrewer")
            pdf("peak-set-partition.pdf",width=16,height=14)
            pal <- brewer.pal(5, "Set1")
            shapes <- rep(x=0:14,length=as.numeric(length(compoundNames)))
            legend.position <- "bottom"
            if(length(compoundNames)>100) {
              legend.position <- "none"
            }
            glegend <- guide_legend(ncol=5,title.position="top",title="Compound")
            print(
                ggplot(data=A,aes(x=t1,y=t2,colour=name,shape=name))+
                geom_point(alpha=0.5,size=3)+
                scale_shape_manual(values=shapes,limits=compoundNames,drop=FALSE,labels=compoundNames)+
                scale_color_manual(values = rep(c(brewer.pal(5, "Set1"),brewer.pal(9,"Set1")[7:9]),length=length(compoundNames)),drop=FALSE,limits=compoundNames,labels=compoundNames)+
                guides(shape = glegend,colour = glegend)+
                #scale_colour_manual(values=c("Unassigned"=pal[2],"Mmdg Exclusive"=pal[3],"mSPA Exclusive"=pal[4],"Common"=rgb(112,153,117,maxColorValue=255)))+
                xlab("RT1")+
                ylab("RT2")+
                xlim(0,rawData$maxX)+ylim(0,rawData$maxY)+
                facet_grid(~set)+
                theme(legend.position=legend.position)
            )
            dev.off()
        }
        setwd(oldwd)
    }
    #NH: modified return value to contain filenames and compoundNames
    list(trt,urt,fnames,compoundNames)
}

#NH: create group-wise statistics for compounds
createCompoundGroupStats <- function(peaks,compoundNames,sd1Thres,sd2Thres) {
    #NH: a dataframe holding stats for each peak group as determined by their name
    compoundGroupStats <- data.frame()
    totalSuspicious <- 0
    #cat("Checking ",length(compoundNames)," compound names\n")
    for(j in 1:length(compoundNames)) {
        name <- compoundNames[[j]]
        #NH: Find matching positions in peaks for name name
        upos <- match(as.character(peaks$name),as.character(name))
        upos <- which(!is.na(upos))
        size <- length(upos)
        t1Stdev <- 0
        t2Stdev <- 0
        t1Mean <- NA
        t2Mean <- NA
        t1Median <- t1Mean
        t2Median <- t2Mean
        size <- length(upos)
        areaStdev <- 0
        #NH: A suspicious group of peaks is defined by not lying within the
        #2D ellipse (fixed 2D Gaussian boundary) defined by sd1Thres and sd2Thres
        #and centered on 0,0
        suspicious <- FALSE
        if(size>1) {
            #cat("Found peak ",name," ",size," times\n")
            t1Stdev <- sd(peaks[upos,]$t1)
            t2Stdev <- sd(peaks[upos,]$t2)
            t1Mean <- mean(peaks[upos,]$t1)
            t2Mean <- mean(peaks[upos,]$t2)
            t1Median <- median(peaks[upos,]$t1)
            t2Median <- median(peaks[upos,]$t2)
            areaStdev <- sd(peaks[upos,]$area)
            #cat("t1Stdev: ",t1Stdev," t2Stdev: ",t2Stdev," sd1Thres: ",sd1Thres," sd2Thres: ",sd2Thres,"\n")
            suspicious <- !isInEllipse2d(t1Stdev,t2Stdev,0,0,sd1Thres,sd2Thres)
            if(suspicious) {
              totalSuspicious <- totalSuspicious+1
            }
            peaks[upos,]$suspicious <- suspicious
        }else if(size==1){
            #cat("Peak is unique: ",name,"\n")
            t1Mean <- peaks[upos,]$t1[1]
            t2Mean <- peaks[upos,]$t2[1]
            t1Median <- t1Mean
            t2Median <- t2Mean
            peaks[upos,]$suspicious <- TRUE
            totalSuspicious <- totalSuspicious+1
        }else{
          cat("Skipping unmatched peak ",name,"\n")
        }
        if(size!=0) {
          #cat("Adding compound group ",j,"\n")
          compoundGroupStats <- rbind(compoundGroupStats,data.frame(name=name,t1Stdev=t1Stdev,t2Stdev=t2Stdev,t1Mean=t1Mean,t2Mean=t2Mean,t1Median=t1Median,t2Median=t2Median,areaStdev=areaStdev,size=size,suspicious=suspicious))
        }
    }
    cat("Number of compound groups: ",length(compoundNames)," Suspicious: ",totalSuspicious,"\n")
    list(stats=compoundGroupStats,peaks=peaks)
}

#NH: remove groups that were marked as outliers
# resolve per-file peakname collisions using minimum distance to group median
# added variant: one of "mSPA", "magm" (max area, smallest distance to group median),
# Please note: !!! trt is altered by this method !!!
resolveAmbiguousPeaks <- function(trt,fnames,compoundNames,compoundGroupStats,variant="mSPA",sd1Thres,sd2Thres) {
    urt <- list()
    cat("Removing suspicious groups using '",variant,"' variant\n")
    ids <- data.frame()
    idsRemoved <- data.frame()
    for(i in 1:length(trt)){
        #NH: For each data.frame/chromatof report
        utmp = trt[[i]]
        cat("Processing report for ",fnames[[i]],"\n")
        #NH: For each compound name
        uniques <- 0
        maxArea <- 0
        suspicious <- 0
        minDist <- 0
        ftmp <- utmp
        rows <- nrow(utmp)
        removed <- 0
        #cat("Received ",nrow(utmp)," peaks\n")
        for(j in 1:length(compoundNames)){
            #NH: select the next row with name
            #cat("Locating ",compoundNames[[j]],"\n")
            upos = match(as.character(ftmp$name),as.character(compoundNames[[j]]))
            upos = which(!is.na(upos))
            compGroupStats <- compoundGroupStats[compoundGroupStats$name==compoundNames[[j]],]
            #cat("Found at positions ",upos,"\n")
            #cat("CompoundGroupStats for ",compGroupStats$name,"\n")
            #cat("Rows: ",nrow(ftmp)," original rows: ",rows," removed: ",removed," uniques: ",uniques,"\n")
            #NH: until there are no more rows left with that name
            #if(length(upos)>1){
            if(length(upos)>0){
                if(variant=="mSPA" || variant=="mspa") {
                    sutmp = ftmp[upos,]
                    tpos = which(sutmp$area==max(sutmp$area))[1]
                    dpos = upos[-tpos]
                    if(length(dpos)>0) {
                        idsRemoved <- rbind(idsRemoved, data.frame(id=createCompoundId(fnames[[i]],ftmp[dpos,]$id,ftmp[dpos,]$name)))
                        removed <- removed + length(dpos)
                        ftmp = ftmp[-dpos,]
                        if(length(dpos)==1) {
                          uniques <- uniques + 1
                        }else{
                          #cat("Removing peaks with non-maximal area: ",dpos,"\n")
                          maxArea <- maxArea+1
                        }

                    }
                    ftmp$name = factor(ftmp$name)
                    ids <- rbind(ids, data.frame(id=createCompoundId(fnames[[i]],ftmp[tpos,]$id,ftmp[tpos,]$name)))
                    #cat("Length: ",nrow(ftmp),"\n")
                    stopifnot(nrow(ftmp)==(rows-removed))
                } else if(variant=="mgma") {
                    stopifnot(nrow(ftmp)==(rows-removed))
                    if(compGroupStats$suspicious) {
                        cat("Removing potential outlier group: ",compoundNames[[j]]," at positions ",paste(upos,sep=","),"\n")
                        #NH: remove all other remaining rows with that name from ftmp
                        removed <- removed + length(upos)
                        idsRemoved <- rbind(idsRemoved, data.frame(id=createCompoundId(fnames[[i]],ftmp[upos,]$id,ftmp[upos,]$name)))
                        ftmp = ftmp[-upos,]
                        suspicious <- suspicious+(length(upos))
                    } else {
                        sutmp = ftmp[upos,]
                        tpos = which(sutmp$area==max(sutmp$area))[1]
                        dpos = upos[-tpos]
                        if(length(dpos)>0) {
                          idsRemoved <- rbind(idsRemoved, data.frame(id=createCompoundId(fnames[[i]],ftmp[dpos,]$id,ftmp[dpos,]$name)))
                          removed <- removed + length(dpos)
                          ftmp = ftmp[-dpos,]
                          if(length(dpos)==1) {
                            uniques <- uniques + 1
                          }else{
                            #cat("Removing peaks with non-maximal area: ",dpos,"\n")
                            maxArea <- maxArea+1
                          }
                        }
                        ftmp$name = factor(ftmp$name)
                        ids <- rbind(ids, data.frame(id=createCompoundId(fnames[[i]],ftmp[tpos,]$id,ftmp[tpos,]$name)))
                        stopifnot(nrow(ftmp)==(rows-removed))
                    }
                    spos = match(as.character(ftmp$name),as.character(compoundNames[[j]]))
                    spos = which(!is.na(spos))
                    ftmp$name = factor(ftmp$name)
                    stopifnot(length(spos)==0 || length(spos)==1)
                } else {
                  cat("Received unknown variant ",variant,"\n")
                  stop("!!! Unknown variant, valid ones are mSPA or mspa, and mgma !!!")
                }
            } else if(length(upos)==0) {
                #cat("Compound not found: ",compoundNames[[j]],"\n")
            } else {
                stop("!!! Unexpected result for compound name matching !!!")
            }
        }
        cat("Removed ",nrow(utmp)-nrow(ftmp)," ambiguous peaks of ",nrow(utmp)," peaks\n")
        cat("Found ",uniques," unique names of ",length(compoundNames)," compound names!\n")
        cat("Removed ",suspicious," peaks outside of decision boundary!\n")
        cat("Resolved reference peak via minimum euclidean distance to group median ",minDist," times!\n")
        cat("Resolved reference peak via maximum area ",maxArea," times!\n")
        if(nrow(ftmp)-length(unique(ftmp$name)) > 0) {
            print(ftmp$name)
            stop("!!! Mismatch in compound name vector lengths: Something went wrong during reference generation !!!")
        }
        #NH: add condensed 'reference' report one
        urt[[i]] <- ftmp
    }
    list(urt,ids,idsRemoved)
}

#NH: Added to write out reference multiple alignment
createReferenceAlignment <- function(compoundNames,urt,outputPath,fnames,plot=FALSE,maxX,maxY,nameSuffix="") {
    #compoundNames <- unique(compoundNames)
    #cat(length(compoundNames)," unique compound names\n")
    #print(compoundNames)
    cat("Creating reference alignment\n")
    groupMedianRt <- data.frame()
    referenceAlignment <- data.frame()
    referenceAlignmentRt <- data.frame()
    referenceAlignmentRt1 <- data.frame()
    referenceAlignmentRt2 <- data.frame()
    referenceAlignmentT1T2 <- data.frame()
    referenceAlignmentNames <- data.frame()
    for(j in 1:length(compoundNames)){
        indices <- c()
        t1t2 <- data.frame()
        rts <- c()
        rt1 <- c()
        rt2 <- c()
        cnames <- c()
        matchCnt <- 0
        for(i in 1:length(urt)){
            #cat("Checking compound ",compoundNames[[j]]," with report ",i,"\n")
            upos <- match(as.character(urt[[i]]$name),as.character(compoundNames[[j]]))
            upos <- which(!is.na(upos))
            #print(upos)
#               print(urt[[i]]$ridx[upos[1]])
            filename <- fnames[[i]]
            if(length(upos)==1){
                indices <- cbind(indices,as.character(urt[[i]][upos[1],]$ridx))
                rts <- cbind(rts,as.numeric(urt[[i]][upos[1],]$t1+urt[[i]][upos[1],]$t2))
                rt1 <- cbind(rt1,as.numeric(urt[[i]][upos[1],]$t1))
                rt2 <- cbind(rt2,as.numeric(urt[[i]][upos[1],]$t2))
                cnames <- cbind(cnames,compoundNames[[j]])
                t1t2 <- rbind(t1t2,data.frame(name=as.factor(compoundNames[[j]]),t1=urt[[i]][upos[1],]$t1,t2=urt[[i]][upos[1],]$t2,filename=filename))
                matchCnt <- matchCnt+1
            }else if(length(upos)>1){
                stop(paste("!!! Found more than one instance of ",compoundNames[[j]]," in ", filename,"! This indicates an error in reference generation !!!",sep=""))
                #cat("Found more than one occurrence for compound ",compoundNames[[j]]," in file ",fnames[[i]],"\n")
                #indices <- cbind(indices,as.character("+"))
            }else{
                #NH: This means the compound was not found in the current dataset
                indices <- cbind(indices,as.character("-"))
                rts <- cbind(rts,NA)#as.character("-"))
                rt1 <- cbind(rt1,NA)
                rt2 <- cbind(rt2,NA)
                cnames <- cbind(cnames,"-")
                #t1t2 <- rbind(t1t2,data.frame(name=compoundNames[[j]],t1=NA,t2=NA,filename=filename))
            }
        }
        #only allow rows with at least two matches (pair)
        if(matchCnt > 1) {
            referenceAlignment <- rbind(referenceAlignment,indices)
            referenceAlignmentRt <- rbind(referenceAlignmentRt,rts)
            referenceAlignmentNames <- rbind(referenceAlignmentNames,cnames)
            referenceAlignmentT1T2 <- rbind(referenceAlignmentT1T2,t1t2)
            referenceAlignmentRt1 <- rbind(referenceAlignmentRt1,rt1)
            referenceAlignmentRt2 <- rbind(referenceAlignmentRt2,rt2)
            groupMedianRt <- rbind(groupMedianRt,as.numeric(median(rts,na.rm=T)))
        }
    }
    colnames(referenceAlignment) <- fnames
    colnames(referenceAlignmentRt) <- fnames
    colnames(referenceAlignmentRt1) <- fnames
    colnames(referenceAlignmentRt2) <- fnames
    colnames(referenceAlignmentNames) <- fnames
    colnames(groupMedianRt) <- c("medianRt")
    cat(nrow(referenceAlignment)," unique alignment rows\n")
    #NH: store working dir
    oldwd <- getwd()
    #NH create output directory
    if(!file.exists(outputPath)) {
        output = dir.create(outputPath,recursive=TRUE)
    }
    #NH: restore name<->shape mapping
    referenceAlignmentT1T2.sorted <- referenceAlignmentT1T2[order(referenceAlignmentT1T2$name),]
    shapes <- rep(x=0:14,length=as.numeric(length(compoundNames)))
    #NH: set working dir
    setwd(outputPath)
    if(plot) {
        if(nameSuffix!="") {
            pdf(paste(paste("reference-groups-rt1-rt2",nameSuffix,sep="-"),".pdf",sep=""),width=20,height=12)
        }else{
            pdf("reference-groups-rt1-rt2.pdf",width=20,height=12)
        }
        #cat("ReferenceNames: ",paste(compoundNames,sep="|"),"\n")
        #cat("Names: ",paste(name,sep="|"),"\n")
        #print(shapes)
        legend.position <- "bottom"
        if(length(compoundNames)>100) {
          legend.position <- "none"
        }
        glegend <- guide_legend(ncol=5,"Compound",title.position="top")
        print(
            ggplot(data=referenceAlignmentT1T2,aes(x=t1,y=t2,shape=name,colour=name))+geom_point(alpha=0.5)+
            scale_shape_manual(values=shapes,limits=compoundNames,drop=FALSE,labels=compoundNames)+
            scale_color_manual(values = rep(c(brewer.pal(5, "Set1"),brewer.pal(9,"Set1")[7:9]),length=length(compoundNames)),drop=FALSE,limits=compoundNames,labels=compoundNames)+
            guides(shape = glegend,colour = glegend)+
            #scale_shape_manual(values=shapes,limits=compoundNames,drop=FALSE,labels=compoundNames)+
            #guides(shape = glegend,colour=glegend)+
            xlab("RT1")+ylab("RT2")+
            xlim(0,maxX)+ylim(0,maxY)+
            theme(legend.position=legend.position)
        )
        dev.off()
    }
    referenceAlignment.sorted <- referenceAlignment[order(groupMedianRt$medianRt),]
    referenceAlignmentRt.sorted <- referenceAlignmentRt[order(groupMedianRt$medianRt),]
    referenceAlignmentNames.sorted <- referenceAlignmentNames[order(groupMedianRt$medianRt),]
    referenceAlignmentRt1.sorted <- referenceAlignmentRt1[order(groupMedianRt$medianRt),]
    referenceAlignmentRt2.sorted <- referenceAlignmentRt2[order(groupMedianRt$medianRt),]
    if(nameSuffix!="") {
        write.table(referenceAlignment.sorted,file=paste(paste("reference-alignment",nameSuffix,sep="-"),".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
        write.table(referenceAlignmentRt.sorted,file=paste(paste("reference-alignment-rt",nameSuffix,sep="-"),".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
        write.table(referenceAlignmentRt1.sorted,file=paste(paste("reference-alignment-rt1",nameSuffix,sep="-"),".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
        write.table(referenceAlignmentRt2.sorted,file=paste(paste("reference-alignment-rt2",nameSuffix,sep="-"),".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
        write.table(referenceAlignmentNames.sorted,file=paste(paste("reference-alignment-names",nameSuffix,sep="-"),".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
    }else{
        write.table(referenceAlignment.sorted,file="reference-alignment.txt",sep="\t",quote=FALSE,row.names=FALSE)
        write.table(referenceAlignmentRt.sorted,file="reference-alignment-rt.txt",sep="\t",quote=FALSE,row.names=FALSE)
        write.table(referenceAlignmentRt1.sorted,file="reference-alignment-rt1.txt",sep="\t",quote=FALSE,row.names=FALSE)
        write.table(referenceAlignmentRt2.sorted,file="reference-alignment-rt2.txt",sep="\t",quote=FALSE,row.names=FALSE)
        write.table(referenceAlignmentNames.sorted,file="reference-alignment-names.txt",sep="\t",quote=FALSE,row.names=FALSE)
    }
    #NH: reset working dir
    setwd(oldwd)
}

# load the data
#pad = load.data(pattern="^Standard_*")
#pad = load.data(pattern="^Deblank_*")

#NH: Added for decision boundary definition
isInEllipse2d <- function(x,y,x0,y0,a,b) {
    z <- ((x-x0)^2/a^2) + ((y-y0)^2/b^2)
    isInside <- TRUE
    if(z>1) {
        isInside <- FALSE
    }
    #cat("Is outlier ",outlier," with value: ",z,"\n")
    isInside
}

#NH: Added for decision boundary definition
ellipse2d <- function(x,x0=0,a=1,b=1) {
    b*sqrt(1-((x-x0)/a)^2)
}

###################
#  PEAL ALIGNMENT #
###################

# similarity product
sim.pd <- function(d1,d2,ms.pos,sm="dot"){
    d1 = d1[,ms.pos]
    d2 = d2[,ms.pos]
    if(sm=="dot"){
        dd1 = as.matrix(d1)/sqrt(apply((as.matrix(d1))^2,1,sum))
        dd2 = as.matrix(d2)/sqrt(apply((as.matrix(d2))^2,1,sum))
        sim = dd1 %*% t(dd2)
    }else if(sm=="pearson"){
        sim = cor(t(d1),t(d2))
    }
}

# distance from R
dist.r <- function(d1,d2,dm="canberra"){
    # dm = "euclidean", "maximum", "manhattan", "canberra"

    d1 = d1
    d2 = d2
    len1 = dim(d1)[1]
    len2 = dim(d2)[1]

    da = rbind(d1,d2)
    da = as.matrix(da)
    tdist = as.matrix(dist(da,method=dm,diag=TRUE,upper=TRUE))[(1:len1),((len1+1):(len1+len2))]
    tdist
}

# accuracy calculation
cscore <- function(evalData,id1,id2,d1,d2){
    s1 <- evalData[[id1]]
    s2 <- evalData[[id2]]
    s1.name = as.character(s1$name)
    s2.name = as.character(s2$name)
    # match names of d1 to names of d2, should be unique, maximum length = |d1|
    total.pos = length(na.omit(match(s1.name,s2.name)))
    total.neg = length(s1.name)*length(s2.name) - total.pos

    d1 = d1[d1$nflag!=0,] # ref
    d2 = d2[d2$nflag!=0,] # tar
    d1 = d1[order(d1$nflag),]
    d2 = d2[order(d2$nflag),]
    d1$name = as.character(d1$name)
    d2$name = as.character(d2$name)
    total.match = dim(d1)[1]
    t.p = 0
    for(i in 1:total.match){
        if(d1$name[i] == d2$name[i]){
            t.p = t.p + 1
        }
    }
    f.p = total.match - t.p
    f.n = total.pos - t.p
    t.n = total.neg - f.p
    t.p.r = t.p/total.pos # t.p/t.p+f.p
    p.p.v = t.p/(t.p+f.p)
    f1 = 2*t.p.r*p.p.v/(t.p.r+p.p.v)

    rlt=c(t.p.r,p.p.v,t.p,f.p,f.n,t.n,f1)
    rlt
}

# Peak alignment
mspa <- function(
        padata = pad
        ,evalData = NULL
        ,rid = 1
        ,tid = 2
        ,method="PMWO"
        ,w = .5
        ,k = 5
        ,rho = .95
        ,distance="canberra"
        ,similarity="dot"
        ,vn = c("t1","t2")
        ,nseed = 3927
        ,sp.pos=NA
        ,anal=FALSE
        ,plot=FALSE
    ){
    # padata: peak alignment data; list
    # rid: index of the reference sample
    # tid: index of the target sample
    # distance (dm): "euclidean", "maximum", "manhattan", "canberra"
    # similarity (sm): "dot", "pearson"
    # method (pm): peak alignment method
    #   =1: Only Distance (PAD)
    #   =2: Only Similarity (PAS)
    #   =3: Distance window => Similarity (DW-PAS)
    #   =4: Similarity window => Distance (SW-PAD)
    #   =5: Mixture similarity (PAM)
    #   =6: Optimal version of Method 5 (OP-PAM)
    # anal: T if both data have the name information. i.e., identification information
    # sp.pos: positions of mass spectra in the dataset
    # vn: positions of RT1 and RT2
    # k (kdist): cutoff for the distance-based window (for pm=3)
    # rho (cv): cutoff for the similarity-based window (for pm=4)
    # w: weight for distance (for pm=5) and the initial weight for pm=6
    id1 = rid
    id2 = tid
    cv = rho
    kdist = k
    sm = similarity
    dm = distance
    if(method=="PAD"){
        pm = 1
    }else if(method=="PAS"){
        pm = 2
    }else if(method=="DW-PAS"){
        pm = 3
    }else if(method=="SW-PAD"){
        pm = 4
    }else if(method=="PAM"){
        pm = 5
    }else if(method=="OP-PAM"){
        pm = 6
    }
    w = as.numeric(w)

    set.seed(nseed)

    refd = padata[[id1]]
    tard = padata[[id2]]

    if(is.na(sp.pos)){
        sp.pos = sort(grep("X",colnames(refd)))
        sp.pos = sp.pos[1]:sp.pos[length(sp.pos)]
    }

    # index
    refd$ndx = c(1:dim(refd)[1])
    tard$ndx = c(1:dim(tard)[1])

    # flag
    refd$nflag = rep(0,dim(refd)[1])
    tard$nflag = rep(0,dim(tard)[1])

    # correlation
    refd$cor = rep(-2,dim(refd)[1])

    # distance
    refd$dist = refd$cor

    # underlying similarity (pm=1,4: distance; others: similarity)
    refd$sim = refd$cor

    # number of compounds
    lenr = dim(refd)[1]
    lent = dim(tard)[1]

    lntrans <- function(vdata){
        a = min(vdata)
        b = max(vdata)
        ndata = vdata/(b-a)-a/(b-a)
        ndata
    }

    # sort by RT1 and RT2
    refd = refd[order(refd[,vn[1]],refd[,vn[2]]),] # library
    tard = tard[order(tard[,vn[1]],tard[,vn[2]]),] # query

    # distance
    avec.dist = NULL
    avec.dist = dist.r(tard[,vn],refd[,vn],dm=dm)
    avec.dist = 1/(1+avec.dist)

    # spectra similarity
    amat.cor = NULL
    amat.cor = sim.pd(tard[,],refd[,],ms.pos=sp.pos,sm=sm)

    # search all the samples given a reference
    if(pm!=6){
        for(i in 1:lent){

            sim2 = -2

            vec.dist = c(avec.dist[i,])
            mat.cor = c(amat.cor[i,])

            if(pm==1){ # distance
                k.dist.pos = 1:length(vec.dist)
                mat.cor.dist = vec.dist
                max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
                max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
                max.cor.k.pos = k.dist.pos[max.cor.pos]
            }else if(pm==2){ # similarity
                cv.cor.pos = 1:length(mat.cor)
                mat.cor.dist = mat.cor
                max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
                max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
                max.cor.k.pos = cv.cor.pos[max.cor.pos]
            }else if(pm==3){ # distance-window
                if(kdist<=-1){ # i.e.: pm=1
                    k.dist.pos = 1:length(vec.dist)
                }else if(kdist<dim(tard)[1]){
                    k.dist.pos = which(vec.dist>=sort(vec.dist,decreasing=T)[kdist])
                }else{
                    k.dist.pos = 1:length(vec.dist)
                }

                smat.cor = mat.cor[k.dist.pos]
                cv.cor.pos = 1:length(smat.cor)
                mat.cor.dist = smat.cor[cv.cor.pos]
                max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
                max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
                max.cor.k.pos = k.dist.pos[cv.cor.pos[max.cor.pos]]
            }else if(pm==4){ # similarity-window
                if(cv<=-1){ # i.e.: pm=1
                    cv.cor.pos = 1:length(mat.cor)
                }else{
                    cv.cor.pos = which(mat.cor>=cv)
                }
                if(length(cv.cor.pos)==0)
                    next

                svec.dist = vec.dist[cv.cor.pos]
                k.dist.pos = 1:length(svec.dist)
                mat.cor.dist = svec.dist[k.dist.pos]
                max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
                max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
                max.cor.k.pos = cv.cor.pos[k.dist.pos[max.cor.pos]]
            }else if(pm==5){ # mixture
                k.dist.pos = 1:length(vec.dist)
                mat.cor.dist = w*vec.dist + (1-w)*mat.cor
                max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
                max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
                max.cor.k.pos = k.dist.pos[max.cor.pos]
            }

            sim2 = mat.cor.dist[max.cor.pos]

            if(refd$nflag[max.cor.k.pos]!=0){
                sim1 = refd$sim[max.cor.k.pos]
                if(sim2>sim1){
                    tard$nflag[tard$ndx==refd$nflag[max.cor.k.pos]] = 0
                    refd$nflag[max.cor.k.pos] = tard$ndx[i]
                    refd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
                    refd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
                    refd$sim[max.cor.k.pos] = sim2
                    tard$nflag[i] = tard$ndx[i]
                }
            }else{
                refd$nflag[max.cor.k.pos] = tard$ndx[i]
                refd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
                refd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
                refd$sim[max.cor.k.pos] = sim2
                tard$nflag[i] = tard$ndx[i]
            }
        }
    }else{

        fun <- function(para,edm="canberra"){
            avec.dist = NULL
            avec.dist = dist.r(tard[,vn],refd[,vn],dm=edm)
            avec.dist = 1/(1+avec.dist)

            irefd = refd
            itard = tard

            for(i in 1:lent){
                sim2 = -2

                vec.dist = c(avec.dist[i,])
                mat.cor = c(amat.cor[i,])

                k.dist.pos = 1:length(vec.dist)
                mat.cor.dist = para*vec.dist + (1-para)*mat.cor
                max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
                max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
                max.cor.k.pos = k.dist.pos[max.cor.pos]

                sim2 = mat.cor.dist[max.cor.pos]

                if(refd$nflag[max.cor.k.pos]!=0){
                    sim1 = refd$sim[max.cor.k.pos]
                    if(sim2>sim1){
                        itard$nflag[itard$ndx==irefd$nflag[max.cor.k.pos]] = 0
                        irefd$nflag[max.cor.k.pos] = itard$ndx[i]
                        irefd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
                        irefd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
                        irefd$sim[max.cor.k.pos] = sim2
                        itard$nflag[i] = itard$ndx[i]
                    }
                }else{
                    irefd$nflag[max.cor.k.pos] = itard$ndx[i]
                    irefd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
                    irefd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
                    irefd$sim[max.cor.k.pos] = sim2
                    itard$nflag[i] = itard$ndx[i]
                }
            }
            all.sim = irefd$sim[irefd$nflag!=0]
            all.cor = irefd$cor[irefd$nflag!=0]
            all.dist = irefd$dist[irefd$nflag!=0]

            s.sim = sum(all.sim)
            s.cor = sum(all.cor)
            s.dist = sum(all.dist)
            s.mix = s.sim + s.cor + s.dist
            ll = s.mix

            -ll
        }

        w1 = nlminb(w,fun,lower=0,upper=1,edm="euclidean")
        #cat("    *** euclidean = ",w1$par, "(",w1$obj,")\n")
        w2 = nlminb(w,fun,lower=0,upper=1,edm="maximum")
        #cat("    *** maximum = ",w2$par, "(",w2$obj,")\n")
        w3 = nlminb(w,fun,lower=0,upper=1,edm="manhattan")
        #cat("    *** manhattan = ",w3$par, "(",w3$obj,")\n")
        w4 = nlminb(w,fun,lower=0,upper=1,edm="canberra")
        #cat("    *** canberra = ",w4$par, "(",w4$obj,")\n")

        objs = c(w1$objective,w2$objective,w3$objective,w4$objective)
        ws = c(w1$par,w2$par,w3$par,w4$par)

        mim.pos = which(objs == min(objs))[1]
        w = ws[mim.pos]

        if(mim.pos==1){
            edm = "euclidean"
        }else if(mim.pos==2){
            edm = "maximum"
        }else if(mim.pos==3){
            edm = "manhattan"
        }else{
            edm = "canberra"
        }

        #cat("     *** ew=",w," at ",edm,"(",mim.pos,")\n")

        avec.dist = NULL
        avec.dist = dist.r(tard[,vn],refd[,vn],dm=edm)
        avec.dist = 1/(1+avec.dist)

        for(i in 1:lent){

            sim2 = -2

            vec.dist = c(avec.dist[i,])
            mat.cor = c(amat.cor[i,])

            k.dist.pos = 1:length(vec.dist)
            mat.cor.dist = w*vec.dist + (1-w)*mat.cor
            max.cor.dist = which(mat.cor.dist==max(mat.cor.dist))
            max.cor.pos = max.cor.dist[sample(1:length(max.cor.dist),1)]
            max.cor.k.pos = k.dist.pos[max.cor.pos]

            sim2 = mat.cor.dist[max.cor.pos]

            if(refd$nflag[max.cor.k.pos]!=0){
                sim1 = refd$sim[max.cor.k.pos]
                if(sim2>sim1){
                    tard$nflag[tard$ndx==refd$nflag[max.cor.k.pos]] = 0
                    refd$nflag[max.cor.k.pos] = tard$ndx[i]
                    refd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
                    refd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
                    refd$sim[max.cor.k.pos] = sim2
                    tard$nflag[i] = tard$ndx[i]
                }
            }else{
                refd$nflag[max.cor.k.pos] = tard$ndx[i]
                refd$cor[max.cor.k.pos] = mat.cor[max.cor.k.pos]
                refd$dist[max.cor.k.pos] = vec.dist[max.cor.k.pos]
                refd$sim[max.cor.k.pos] = sim2
                tard$nflag[i] = tard$ndx[i]
            }
        }
    }

    all.sim = refd$sim[refd$nflag!=0]
    all.cor = refd$cor[refd$nflag!=0]
    all.dist = refd$dist[refd$nflag!=0]

    s.sim = sum(all.sim)
    s.cor = sum(all.cor)
    s.dist = sum(all.dist)

    s.mix = s.sim + s.cor + s.dist
    ll = -s.mix


    refr = refd[refd$nflag!=0,]
    tarr = tard[tard$nflag!=0,]
    refr = refr[order(refr$nflag),]
    tarr = tarr[order(tarr$nflag),]

    atable = data.frame(refr$ndx,tarr$ndx)
    dimnames(atable)[[2]] = paste("S",c(id1,id2),sep="")
    ridxtable = data.frame(refr$ridx,tarr$ridx)
    dimnames(ridxtable)[[2]] = paste("S",c(id1,id2),sep="")

    if(plot){
        plot(refd[refd$nflag==0,vn[1]],refd[refd$nflag==0,vn[2]],cex=.5,pch=19,col=3
            ,xlim=c(min(refd[,vn[1]],tard[,vn[1]]),max(refd[,vn[1]],tard[,vn[1]]))
            ,ylim=c(min(refd[,vn[2]],tard[,vn[2]]),max(refd[,vn[2]],tard[,vn[2]]))
            ,xlab="RT1",ylab="RT2"
            )
        points(tard[tard$nflag==0,vn[1]],tard[tard$nflag==0,vn[2]],col=5,pch=19,cex=.5)

        points(refr[,vn[1]],refr[,vn[2]],cex=1,pch=3,col=1)
        points(tarr[,vn[1]],tarr[,vn[2]],col=2,pch=4,cex=1)

        for(i in 1:dim(refr)[1]){
            segments(refr[i,vn[1]],refr[i,vn[2]],tarr[i,vn[1]],tarr[i,vn[2]],col=4,lwd=.5)
        }
        legend("topleft",c("ref:nonmatch","tar:nonmatch","ref:match","tar:match","pair")
                ,lty=c(-1,-1,-1,-1,1),col=c(3,5,1,2,4),pch=c(19,19,3,4,-1)
                ,bty="n"
            )
    }

    if(anal){
        accur = cscore(evalData,id1,id2,refd,tard)
        summary = c(id1,id2,pm,dm,kdist,sm,cv,w,as.numeric(accur),ll)
        summary = data.frame(t(summary))
        dimnames(summary)[[2]] = c(
                    "id1","id2"
                    ,"method","distance","k","similarity","cv","w"
                    ,"tpr","ppv","tp","fp","fn","tn"
                    ,"f1","ll"
                )
    }else{
        summary = c(id1,id2,pm,dm,kdist,sm,cv,w,ll)
        summary = data.frame(t(summary))
        dimnames(summary)[[2]] = c(
                    "id1","id2"
                    ,"method","distance","k","similarity","cv","w"
                    ,"ll"
                )   }

    list(table=atable,summary=summary,ridxtable=ridxtable)
}

# alingment of all peak lists
align <- function(padata,alist,method="peak",plot=FALSE){
    # alist: it should be arranged.
    # e.g.: (1,2), (2,3), (3,4), (4,5), ...
    # Its order sould be consistent with that of padata
    # method: peak, name, area

    len = length(alist)
    wtable = alist[[1]]
    for(i in 2:len){
      ta = alist[[i]]
      ##AP: to allow alignment to file 1 only as ref
      #if(colnames(ta)[1]!=colnames(alist[[i-1]])[2]){
        #stop("!!! Please make sure the order of alist !!!\n")
      #}
      t.pos = match(wtable[,i],ta[,1])
      wtable = data.frame(wtable[which(!is.na(t.pos)),],ta[na.omit(t.pos),2])
      dimnames(wtable)[[2]][i+1] = colnames(ta)[2]
    }
    wlen = dim(wtable)[2]

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

# PPV calculation based on the name table
comp.name <- function(ntable){
    nt = factor(ntable)
    len = dim(ntable)
    total = len[1]
    correct = 0
    for(i in 1:len[1]){
        flag = 1
        for(j in 2:len[2]){
            if(ntable[i,j-1]!=ntable[i,j]){
                flag=0
                break
            }
        }
        if(flag==1){
            correct = correct + 1
        }
    }
    c(TP=correct,P=total,PPV=correct/total)
}
