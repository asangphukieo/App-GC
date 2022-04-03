
#GC method

#similar to removeDoubleMasses method in metaMS, but take 4 columns (mz, I, rt, Area)
#Since in nominal-mass GC data m/z values are rounded to whole numbers, in some cases an m/z value may occur multiple times - in that case the mean of the intensities is kept (and the retention times are averaged). removeDoubleMasses takes a list of spectra in the form of a three-column matrix, (mz, I, rt), summing the intensities and averaging the retention times for multiple identical masses. Not meant to be called directly by the user.  removeDoubleMasses also do sum the area of the peak with the same mz
removeDoubleMasses2 <- function(spclist) {
  lapply(spclist, function(pspc) {
    if (max(table(pspc[,1])) > 1) {
      if (ncol(pspc) == 4) {
        huhn <- cbind(aggregate(pspc[,2], list(pspc[,1]), FUN = sum),
                      aggregate(pspc[,3], list(pspc[,1]), FUN = mean)$x,
                      aggregate(pspc[,4], list(pspc[,1]), FUN = sum)$x )
      } else {
        huhn <- cbind(aggregate(pspc[,2], list(pspc[,1]), FUN = sum))
      }
      colnames(huhn) <- colnames(pspc)
      
      huhn
    } else {
      pspc
    }})
}

to.file_msp <- function(object, file = NULL,
                   settings = NULL, ndigit = 0, 
                   minfeat, 
                   minintens,
                   intensity = "maxo", 
                   secs2mins = TRUE,
                   to2D=FALSE,
                   mod_time=4, 
                   comment="") {
  if (!is.null(settings)) {
    intensity <- settings$intensityMeasure
    minfeat <- settings$minfeat
    minintens <- settings$minintens
  } else {
    intensity <- match.arg(intensity)
  }
  
  #
  if (class(object) == "xsAnnotate") { ## CAMERA annotation object
    allpks <- object@groupInfo
    minI <- minintens * max(allpks[, intensity])
    tooSmall <- which(allpks[, intensity] < minI)
    
    pspectra <- lapply(object@pspectra,
                       function(x) x[!x %in% tooSmall])
  } else { ## a peak table
    minI <- minintens * max(object[, intensity])
    allpks <- object[object[, intensity] >= minI,]
    pspectra <- split(1:nrow(allpks), allpks[,"rt"])
  }


  npeaks <- sapply(pspectra, length)
  pspectra <- pspectra[npeaks >= minfeat]
  #print(pspectra)
  #print('result')
  result <-
    removeDoubleMasses2(lapply(pspectra,
                              function(x)
                              cbind(mz = round(allpks[x,"mz"], digits = ndigit),
                                    allpks[x, c(intensity,"rt","into")] )))

  #print(pspectra[[9]])
  #print(result[[7]])

  ################
  if (!is.null(file)) {
    if (length(pspectra) > 0) {
      mylist=list() #list to export for peak alignment

      for (i in 1:length(pspectra)) {
        ## maximum of 1000 features per file
        #ofile <- paste(file, "_p", ceiling(i / 1000), ".msp", sep = "") ##make new file if exceed 1000 peaks
        ofile <- paste(file, ".msp", sep = "")
        ocsv <- paste(file, "_peak", ".csv", sep = "")

        #newfile <- (i %% 1000) == 1 ##make new file if exceed 1000 peaks
        
        idx <- pspectra[[i]]
        #pks <- allpks[idx,,drop = FALSE]
        pks <- result[[i]][order(result[[i]][,1]),,drop = FALSE]
        ### normalize to range 1-1000
        #pks[,2] <- 1000* as.numeric(pks[,2]) / as.numeric(max(pks[,2])) # normalize intensity to range 1-1000
        #pks[,3] <- 1000* as.numeric(pks[,3]) / as.numeric(max(pks[,3])) # normalize area to range 1-1000
        
	#cat(pks[,3],'\n',sep='')
	#print(pks[order(pks[,2]),,drop=FALSE] )

	#use RT of maximum peak as representative
	order_by_maxo= tail(pks[order(pks[,2]),,drop=FALSE],1)
	max_maxo_RT= order_by_maxo[,3]
	max_maxo = order_by_maxo[,2]
	max_maxo_into = order_by_maxo[,4]
	max_mz = floor(order_by_maxo[,1])
	
	if (to2D==TRUE) {
		#convert to 1st RT and 2nd RT
		RT_1= (max_maxo_RT%/%mod_time)*mod_time
		RT_2= max_maxo_RT - RT_1
		
	} else {
		RT_1= max_maxo_RT
		RT_2="0.00" #not use for 1D gc
	}

        if (secs2mins){
		
		RT_1= (RT_1/ 60)
		RT_unit='min'
		#ave_RT=paste((mean(pks[,3])/ 60),'min',sep='')

	}else{
		
		RT_unit='sec'
		#ave_RT = paste(mean(pks[,3]),'sec',sep='') 
	}
  
	area=max_maxo_into 

	#write msp file 
        cat("Name: peak_id_", i,'\n', sep = "",
            file = ofile, append = TRUE)

        cat("Synon: DATE: ",date(),'\n',sep = "",
            file = ofile, append = TRUE)

        cat("Synon: RT: ",RT_unit,RT_1,sep = "",
            file = ofile, append = TRUE)

	if (to2D==TRUE) {
        cat("Synon: RT2: sec",RT_2,'\n',sep = "",
            file = ofile, append = TRUE)
	}

        cat("Comments: ",comment,'\n',sep = "",
            file = ofile, append = TRUE)

        cat("\nNum Peaks in gr:", nrow(pks),'\n', file = ofile, append = TRUE)


        for (ii in 1:nrow(pks)) {

          if ((ii>1) & (ii%%6)==0 ){
              cat("\n", file = ofile, append = TRUE)             
              }

          cat(" ",round(pks[ii, 1], ndigit)," ",
              round(pks[ii, 2], ndigit), ";",
              file = ofile, append = TRUE, sep = "")
          }
        
        cat("\n\n", file = ofile, append = TRUE)

	#write dataframe file 
	peak_id=paste("peak_id_",i,sep="")

	spec_body=""
	for (ii in 1:nrow(pks)) {
          spec_body = paste(spec_body, paste( round(pks[ii, 1], ndigit),",",round(pks[ii, 2], ndigit), " ", sep = ""), sep = "")
          }
 
	#body_csv=""
	#for (ii in 1:nrow(pks)) {
        #  spec_body = paste(spec_body, paste( round(pks[ii, 1], ndigit),",",round(pks[ii, 2], ndigit), " ", sep = ""), sep = "")
        #  }                   
	#write csv file (use in MSPA peak alignment)
	#csv_write = paste( peak_id, RT_1, RT_2, max_mz, max_maxo, area, body_csv, sep='\t')
	#write.table( csv_write , raw_data ,sep = ""
	#		, row.names = FALSE, col.names = FALSE, append = TRUE,quote = FALSE)

	if (i == 1){
		ali_file = data.frame( 
					"peak_ID"=peak_id,
					"RT" = RT_1, 
					"RT2"= RT_2,
					"mass"= max_mz,
					"height"=max_maxo,
					"area"=area,
					"spectra"=spec_body
					)
	}else{
		ali_file = rbind(ali_file, data.frame(
					"peak_ID"=peak_id,
					"RT" = RT_1, 
					"RT2"= RT_2,
					"mass"= max_mz,
					"height"=max_maxo,
					"area"=area,
					"spectra"=spec_body
		))
	}
	
	write.table(ali_file,ocsv, row.names=FALSE, quote = FALSE,sep='\t')

      mylist[[file]] = ali_file
      }
    }
  }
  mylist
}



erah.to.file_msp <- function(object, file_ID="file_A",output_path="./") {  
  if (class(object) == "MetaboSet") { ## CAMERA annotation object
    cat(paste( "# Erah : Object contains ",length(names(object@Data@FactorList))," files\n",sep=""))
    
    mylist=list()
    
    for (file_i in names(object@Data@FactorList)){
        file=paste(output_path,file_i,sep="/")

        ofile <- paste(file, ".msp", sep = "")
        ocsv <- paste(file, "_peak", ".csv", sep = "")

	for (i in 1:length(object@Data@FactorList[[file_i]]$Spectra)){
               
	RT_unit='min'
	
	RT_1=object@Data@FactorList[[file_i]]$RT[i]
	ori_spec=object@Data@FactorList[[file_i]]$Spectra[i]
	
	SPEC=as.character(ori_spec)
	SPEC=gsub(" ",";",SPEC)
	SPEC=gsub(","," ",SPEC)
	SPEC=strsplit(SPEC,";")[[1]]
	
	Peak_Height=object@Data@FactorList[[file_i]]$"Peak Height"[i]
	Peak_Area=object@Data@FactorList[[file_i]]$"Area"[i]
	
	#write msp file 
        cat("Name: peak_id_", i,'\n', sep = "",
            file = ofile, append = TRUE)

        cat("Synon: DATE: ",date(),'\n',sep = "",
            file = ofile, append = TRUE)

        cat("Synon: RT: ",RT_unit,RT_1,"\n",sep= "",
            file = ofile, append = TRUE)

        cat("Num Peaks in gr:", length(SPEC),'\n', file = ofile, append = TRUE)

	max_MZ=0
	max_height=0
	spec_body=""
        for (ii in 1:length(SPEC)) {
		x_mz= as.numeric(strsplit(SPEC[ii]," ")[[1]][1])
		x_height= as.numeric(strsplit(SPEC[ii]," ")[[1]][2])
		
		if(x_height > max_height){
			max_height = x_height
			max_MZ = x_mz
		}

		if ((ii>1) & (ii%%6)==0 ){
		cat("\n", file = ofile, append = TRUE)             
		}

		cat(' ',x_mz,' ',x_height,";",file = ofile, append = TRUE, sep = "")
		
		spec_body = paste(spec_body,paste( x_mz, ",", x_height," ", sep = ""),sep="")
	
		}
        
        cat("\n\n", file = ofile, append = TRUE)
        
        peak_id=paste("peak_id_",i,sep="")
        RT_2="0.00"
	
        #for (ii in 1:length(SPEC)) {
        #  spec_body = paste(spec_body,paste( SPEC[ii], " ", sep = ""),sep="")
        #  }

	if (i == 1){
		ali_file = data.frame( 
					"peak_ID"=peak_id,
					"RT" = RT_1, 
					"RT2"= RT_2,
					"mass"= max_MZ,
					"height"=Peak_Height,
					"area"=Peak_Area,
					"spectra"=spec_body
					)
	}else{
		ali_file = rbind(ali_file, data.frame(
					"peak_ID"=peak_id,
					"RT" = RT_1, 
					"RT2"= RT_2,
					"mass"= max_MZ,
					"height"=Peak_Height,
					"area"=Peak_Area,
					"spectra"=spec_body
		))
	}
	
	write.table(ali_file,ocsv,row.names=FALSE,quote = FALSE,sep='\t')


      mylist[[file_i]] = ali_file
      }
    } 
    
    mylist
    
  } else { ## a peak table
    cat("# Erah : Error at Method erah.to.file_msp requires MetaboSet as object!")
  }
  
}



summary_result <- function( app_gc_path ='.', ref=NULL, name='by_consensus') {
	#load data from gcpipe output folder
	peak_detect_result = paste(app_gc_path,'/PEAK_DETECTION/save_peak.Rdata',sep='')	
	peak_ali_result = paste(app_gc_path,'/PEAK_ALIGNMENT/save_ali_peak.Rdata',sep='')
	peak_annot_result = paste(app_gc_path,'/PEAK_ANNOTATION/save_annotation_peak.Rdata',sep='')
	load(peak_detect_result) #variable:mypeak
	load(peak_ali_result) #variable:peak_num_aligned
	load(peak_annot_result) #varable: myannot
	#summary_table = data.frame()
	#remove reference sample from summary analysis
	if ( !is.null(ref) ){
		#ref is list of sample names that will be removed from the analysis
		peak_num_aligned = peak_num_aligned[, !colnames(peak_num_aligned) %in% ref]
		myannot = myannot[!names(myannot) %in% ref]
		mypeak = mypeak[!names(mypeak) %in% ref]		
	}
	#remove RT columns
	peak_num_aligned = peak_num_aligned[, !colnames(peak_num_aligned) %in% c("mean_RT","mean_RT2")]
	#collect sample name used for searching in peak annotation result(myannot)
	sample_name_list = colnames(peak_num_aligned)
	for ( peak_ali_ID in 1:nrow(peak_num_aligned)){
		#count peak in peak_ali 
		found_peak = !is.na(peak_num_aligned[peak_ali_ID,,drop=FALSE]) #TRUE or FALSE
		only_found_peak = peak_num_aligned[peak_ali_ID,,drop=FALSE][,found_peak==TRUE,drop=FALSE] #show only found peak_id
		n_peaks = length(only_found_peak)
		percent_peak_found = 100*(n_peaks/(length(sample_name_list)))
		#set to 2 digit
		percent_peak_found=format(round(percent_peak_found,2),nsmall=2)
		#collect peak information that is annotated in all samples
		#and calculate mean RT, mean Hight, mean Area from all samples
		#and calculate the maximum of NIST peak matching of similarity, MF, consensus Compound name from peak annotation of all peak from the peak_ali_ID
		mean_RT = c()
		mean_RT2= c()
		mean_Height = c()
		mean_Area =c()
		#max_Prob =c()
		max_MF_name=data.frame()		
		for (sample_i in colnames(only_found_peak)){
			annot =myannot[[sample_i]]
			match_peak_id = as.character(only_found_peak[,sample_i])
			match_annot_peak = annot[ annot$peak_ID==match_peak_id,]
			#check if peak anotation is found
			if (nrow(match_annot_peak) >0){
				mean_RT = c(mean_RT,as.numeric(as.character(match_annot_peak$RT[1])))
				mean_RT2= c(mean_RT2,as.numeric(as.character(match_annot_peak$RT2[1])))
				mean_Height = c(mean_Height,as.numeric(as.character(match_annot_peak$height[1])))
				mean_Area = c(mean_Area,as.numeric(as.character(match_annot_peak$area[1])))
				#max_Prob = c(max_Prob,match_annot_peak$Prob)
				new_MF_name =	data.frame(match_annot_peak$MF,match_annot_peak$annotation)
				max_MF_name = rbind(max_MF_name,new_MF_name)	
			}else{
				#in case of no annotation found, collect peak information from mypeak table instead
				annot = mypeak[[sample_i]]
				match_annot_peak = annot[ annot$peak_ID==match_peak_id,]
				mean_RT = c(mean_RT,as.numeric(as.character(match_annot_peak$RT[1])))
				mean_RT2= c(mean_RT2,as.numeric(as.character(match_annot_peak$RT2[1])))
				mean_Height = c(mean_Height,as.numeric(as.character(match_annot_peak$height[1])))
				mean_Area = c(mean_Area,as.numeric(as.character(match_annot_peak$area[1])))			
				new_MF_name =	data.frame(match_annot_peak.MF=NA,match_annot_peak.annotation=NA)
				max_MF_name = rbind(max_MF_name,new_MF_name)					
			}
		}		
		mean_RT=format(round(mean(mean_RT),2),nsmall=2)
		mean_RT2=format(round(mean(mean_RT2),2),nsmall=2)
		mean_Height=format(round(mean(mean_Height),2),nsmall=2)
		mean_Area=format(round(mean(mean_Area),2),nsmall=2)
		names(max_MF_name) = c('MF','annotation')		
		max_MF_name = max_MF_name[!is.na(max_MF_name[,1]),] #remove NA
		n_annot = nrow(max_MF_name)	
					
		if (n_annot > 0){	
			#calculate freqency of compound name
			max_MF_name=transform(max_MF_name,freq= ave(seq(nrow(max_MF_name)),annotation, FUN=length))		
			#select final compound name			
			if(name=="by_consensus"){
				max_freq = max_MF_name[order(max_MF_name$freq,decreasing=T),][1,"freq",drop=F]
				max_freq_max_MF_name = max_MF_name[max_MF_name$freq== as.numeric(max_freq),]
				comp_name = max_freq_max_MF_name[order(max_freq_max_MF_name$freq,decreasing=T),][1,]
				percent_consensus = (comp_name$freq / nrow(max_MF_name)) * 100.00
				percent_consensus = format(round(percent_consensus,2),nsmall=2)			
			}else if(name=="by_MF"){
				comp_name = max_MF_name[order(max_MF_name$MF,decreasing=T),][1,]
				percent_consensus = (comp_name$freq / nrow(max_MF_name)) * 100.00
				percent_consensus = format(round(percent_consensus,2),nsmall=2)			
			}
			sum_peaks = data.frame(mean_RT,mean_RT2, mean_Height, mean_Area, as.numeric(comp_name[1,"MF"]), comp_name[1,"annotation"], percent_peak_found, percent_consensus)
		}else{
			sum_peaks = data.frame(mean_RT,mean_RT2, mean_Height, mean_Area, NA, NA, percent_peak_found, NA)
			
		}	
		
		names(sum_peaks) = c('meanRT','meanRT2','meanHeight','meanArea','maxMF','Name','Percent_peak_ali','Percent_name_consensus')
		
		if (peak_ali_ID == 1){
			summary_table = sum_peaks
		}else{
			summary_table = rbind(summary_table,sum_peaks)
		}
	}
	#names(summary_table) = c('meanRT','meanRT2','meanHeight','meanArea','maxMF','Name','Percent_peak_ali','Percent_name_consensus')
	summary_table
}



## to convert detected peak from other platform to erah
## modeified from erah newExp function
Rdata2MetaboSet <- function(mypeak,from2D){	
	sampleID=names(mypeak)
	instrumental.dataframe = data.frame(sampleID)
	factors.list <- lapply(1:nrow(instrumental.dataframe), function(x){as.data.frame(NULL)})
		
	names(factors.list) <- as.vector(instrumental.dataframe$sampleID)
	ident.list <- as.data.frame(matrix(0,ncol=7, dimnames=list(row=0,col= c("AlignID", "tmean", "Name", "MatchFactor", "CAS", "Formula", "DB.Id"))))
	uni.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("Id", "FoldChangue", "pvalue"))))
	multi.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("Id", "CompoundsInvolved", "pvalue"))))

	al.par <- list()
	id.par <- list()
	soft.par <- list()
	
	stat.parameters <- new("MSResultsParameters", Alignment=al.par, Identification=id.par)
	statistics <- new("Statistics", Univariate = uni.stats, Multivariate = multi.stats)
	MS.Results <- new("Results", Parameters = stat.parameters, Identification = ident.list, Statistics = statistics )
	MS.Data <- new("Data", FeatureList = list(NULL), FactorList = factors.list, Parameters = list(NULL))
	MS.MetaData <- new("MetaData", Instrumental = instrumental.dataframe, Phenotype = as.data.frame(NULL), DataDirectory=".")

	sample.container <- new("MetaboSet", Info = "", Data = MS.Data, MetaData = MS.MetaData, Results = MS.Results)
	for(filename in names(mypeak)){
			rownames(mypeak[[filename]]) = NULL
		if (from2D==TRUE) {
			#RT_rev = (RT * 60.0) + RT2 #to convert to 1D RT in second unit
			RT_rev = ( mypeak[[filename]][["RT"]] * 60 ) + mypeak[[filename]][["RT2"]]
			#RT_rev = RT_rev/60
			mypeak[[filename]] = cbind(mypeak[[filename]], RT_rev)
			update_colname = c("ID","RT_rev","RT2","Mass","Peak Height","Area","Spectra","RT")
		}else{
			update_colname = c("ID","RT","RT2","Mass","Peak Height","Area","Spectra")
		}
	
		colnames(mypeak[[filename]]) = update_colname
		sample.container@Data@FactorList[[filename]] = mypeak[[filename]]
	}
	sample.container
}


## to convert detected peak from other platform to erah
## modeified from erah newExp function
Rdata2MetaboSet.prev <- function(mypeak){
	sampleID=names(mypeak)
	instrumental.dataframe = data.frame(sampleID)
	factors.list <- lapply(1:nrow(instrumental.dataframe), function(x){as.data.frame(NULL)})
		
	names(factors.list) <- as.vector(instrumental.dataframe$sampleID)
	ident.list <- as.data.frame(matrix(0,ncol=7, dimnames=list(row=0,col= c("AlignID", "tmean", "Name", "MatchFactor", "CAS", "Formula", "DB.Id"))))
	uni.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("Id", "FoldChangue", "pvalue"))))
	multi.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("Id", "CompoundsInvolved", "pvalue"))))

	al.par <- list()
	id.par <- list()
	soft.par <- list()
	
	stat.parameters <- new("MSResultsParameters", Alignment=al.par, Identification=id.par)
	statistics <- new("Statistics", Univariate = uni.stats, Multivariate = multi.stats)
	MS.Results <- new("Results", Parameters = stat.parameters, Identification = ident.list, Statistics = statistics )
	MS.Data <- new("Data", FeatureList = list(NULL), FactorList = factors.list, Parameters = list(NULL))
	MS.MetaData <- new("MetaData", Instrumental = instrumental.dataframe, Phenotype = as.data.frame(NULL), DataDirectory=".")

	sample.container <- new("MetaboSet", Info = "", Data = MS.Data, MetaData = MS.MetaData, Results = MS.Results)
	for(filename in names(mypeak)){
		colnames(mypeak[[filename]]) = c("ID","RT","RT2","Mass","Peak Height","Area","Spectra")
		sample.container@Data@FactorList[[filename]] = mypeak[[filename]]
	}
	sample.container
}




#wrapper of alignList to export only meanRT with RT, peak area or peak height of each aligned peak
#
erah.export.alignList <- function(ex.align, export=c("Area","RT","Height","ID"), from2D=FALSE, mod_time=4.0){
	aliTable = alignListRT(ex.align, export)
	aliTable = aliTable[,c(3,5:ncol(aliTable)),drop=F]
	header = colnames(aliTable)[2:ncol(aliTable)]
	colnames(aliTable) = c("mean_RT",header)
		
	if (from2D == TRUE) {
		#convert to 1st RT and 2nd RT
		mean_RT= (aliTable[,1,drop=F] %/% mod_time)*mod_time
		mean_RT2= round(aliTable[,1,drop=F] - mean_RT , digits=3)	
		mean_RT= round(mean_RT/60,digits=3) #in second unit
		
		aliTable = cbind(mean_RT2,aliTable[,c(2:ncol(aliTable)),drop=F])
		aliTable = cbind(mean_RT,aliTable)		
		header = colnames(aliTable)[3:ncol(aliTable)]
		colnames(aliTable) = c("mean_RT","mean_RT2",header)		
	}
	
	#Replace 0 with NA
	aliTable[aliTable==0] <- NA
	aliTable
}


#alignListRT was modified from alignList to export RT instead of Area
alignListRT <- function(object, export=c("Area","RT","Height","ID"),from2D=FALSE, mod_time=4.0) {
		
	if(!(any(unlist(lapply(object@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")

	if(export=="Height")
	{
		del.in <- which(colnames(object@Results@Alignment)=="Spectra")
		alList <- object@Results@Alignment[,-del.in]
		if(!is.null(object@Results@Parameters@Alignment$min.samples)){
			alList  <- alList[which(alList$FoundIn>=object@Results@Parameters@Alignment$min.samples),]
		}		
		return(alList)
	}else{		
		align.inds <- as.numeric(as.vector(object@Results@Alignment[,"AlignID"]))
		align.area <- lapply(object@Data@FactorList,function(x) {
			#search.for <- which(x$"AlignID" %in% align.inds)
			search.for <- (sapply(align.inds, function(i) which(as.numeric(as.vector(x$"AlignID"))==i)))

			if(class(search.for)=="list") 
			{
				fill.inds <- unlist(lapply(search.for, length))
				fill.vector <- rep(0,length(align.inds))
				if(export != "ID"){
					fill.vector[fill.inds==1] <- as.numeric(as.vector(x[[export]][unlist(search.for)]))
				}else{
					fill.vector[fill.inds==1] <- as.vector(x[[export]][unlist(search.for)])
				}
			}else{
				if(export != "ID"){
					fill.vector <- as.numeric(as.vector(x[[export]][search.for]))
				}else{
					fill.vector <- as.vector(x[[export]][search.for])
				}
			}
			fill.vector
		})
		
		if(!is.matrix(align.area))
		{
			del.list <- as.vector(which(unlist(lapply(align.area,length))==0))
			if(length(del.list)!=0) align.area <- align.area[-del.list]
			align.area <- unlist(align.area)
			dim(align.area) <- c(length(align.inds), length(align.area)/length(align.inds))
		}
		del.in <- which(colnames(object@Results@Alignment)=="Spectra")
		area.list <- cbind(object@Results@Alignment[,c("AlignID","Factor","tmean","FoundIn")],align.area)
		colnames(area.list) <- colnames(object@Results@Alignment[,-del.in])
		
		if(!is.null(object@Results@Parameters@Alignment$min.samples)){
			area.list  <- area.list[which(area.list$FoundIn>=object@Results@Parameters@Alignment$min.samples),]
		}		
		
		return(area.list)
	}
}




