##################################################
### App-GC pipeline for data processing of GCxGC data  ###
##################################################

########################################
# Author: Apiwat Sangphukieo           #  
# Version: v.1.2                       #
# Email: apiwat.bif@mail.kmutt.ac.th   #
########################################

# Load package
#library(xcms)
library(yaml)
library(ggplot2)

preprocess_2DGCdata <- function(parameter_file='parameter_config.yml') {
	#to extract RT start, RT end , sampling_rate, MZ min, MZmax from reference file and save to Rdata
	#and 
	cat('# App-GC : start 2D-GC data preprocess \n')
	ref.dir <- getwd()
	config = yaml.load_file(parameter_file)
	cdffiles=list.files(path = config$Input_path, pattern= config$Input_file ,full.names=FALSE, ignore.case=TRUE)
	cat('# App-GC : all input files in the path  \n')
	cat('# App-GC : file pattern =',config$Input_file,'\n')
	cat('# App-GC : ',cdffiles,'\n')
	############################

	if ( !is.null(config$GC2D$Peak_detection$reference) ){
		ref_file=config$GC2D$Peak_detection$reference
	}else{
		ref_file=cdffiles[1]
	}
	cat('# App-GC : extract information from reference file  : \n')
	cat('# App-GC : ',ref_file,'\n')

	nc <- xcmsRaw(paste(config$Input_path,'/',ref_file,sep='') )
	RTstart=min(nc@scantime)
	RTstop=max(nc@scantime)
	MZmin=min(nc@mzrange)
	MZmax=max(nc@mzrange)

	mzmin = MZmin
	mzmax = MZmax
	sampling_rate = config$GC2D$Peak_detection$sampling_rate	#from people in the lab, but it can be automaticcally calculated
	modulation_time = config$GC2D$Peak_detection$modulation_time #from people in the lab
	nc=''
	gc() 

	folder_name=gsub('.cdf$','',ref_file)
	temp_path <- paste(config$Output_folder,'/',folder_name,sep='') #make folder
	detected_param = load_2DGCdata(filename=paste(config$Input_path,'/',ref_file,sep=''), 
		mod_time=modulation_time, 
		sam_rate=sampling_rate, 
		exportTIC=T, 
		exportSIC=T, 
		temp_path=temp_path,
		rm_prev_tic=FALSE,
		profMat_method="binlin",
		profMat_step=1.0)

	##### plot reference file #####
	tic.fname = "tic_data.txt"
	# TIC information: original: RT2 by RT1 so transpose necessary
	atic = t(as.matrix(read.table(paste(temp_path,tic.fname,sep='/'))))
	tic.size = dim(atic)
	rownames(atic)=NULL #AS

	d2.df <- reshape2::melt(atic, c("x", "y"), value.name = "z")
	head(d2.df)

	save(RTstart,RTstop,MZmin,MZmax,sampling_rate,modulation_time,file='save_preprocess.Rdata')

	#png(paste(temp_path,'_2Dchrom.pdf',sep=''),res=300)
	p<- ggplot(data=d2.df,aes(x=x,y=y,fill=z) ) + labs(x = "First dimension retention time (s)", y = "Second dimension reteintion time (s)") + geom_tile()
	#p

	#ggsave(paste(temp_path,'_2Dchrom.pdf',sep=''))
	ggsave(paste(temp_path,'_2Dchrom.png',sep=''),dpi=300)

	#optional : for user interactive
	#options(browser = "/usr/bin/google-chrome")
	#ggplotly(p) %>% layout(xaxis=list(fixedrange=TRUE)) %>% layout(yaxis=list(fixedrange=TRUE))

}


#load 2D GCxGC data and export to 2d tic file (total ion chromatogram) and 2d sic file (single ion chromatogram) 
load_2DGCdata <- function(filename='ABC.cdf', profMat_method="binlin",profMat_step=1.0,mod_time=5, sampling_rate=FALSE, exportTIC=T, exportSIC=T, temp_path='./',rm_prev_tic=FALSE) {  
	cat('# App-GC : Remove previous TIC =',rm_prev_tic,'\n')
	if(rm_prev_tic==TRUE){
		cat('# App-GC : Remove previous TIC ',temp_path,'\n')
		unlink(temp_path,recursive=TRUE)
		#file.remove(temp_path)
	}


	if (!dir.exists(temp_path)){
	cat('# App-GC : load cdf file ...\n')

	dir.create(temp_path)
	nc <- xcmsRaw(filename) 
	RTstart=min(nc@scantime)
	RTstop=max(nc@scantime)
	MZmin=min(nc@mzrange)
	MZmax=max(nc@mzrange)
	

	cat('# App-GC : load data file ...\n')

	tic <- rawEIC(nc,rtrange=cbind(RTstart,RTstop),mzrange=cbind(MZmin,MZmax))
	#len_1d <- floor(sampling_rate * mod_time)
	if (sampling_rate==FALSE || is.na(sampling_rate) || sampling_rate == "NA" || is.null(sampling_rate) ){
		sampling_rate <- nc@scantime[1:2]
		sampling_rate <- 1 / abs(diff(sampling_rate)) 
		warning(paste('Detected sampling rate is >> ',sampling_rate,'\nIf this value is incorrect, you should put it manually!!'))
		
	}
	#profile matrix
	#"linbase"uses linear interpolation to impute values for empty elements 
	profmat <- profMat(nc, step = profMat_step, method = profMat_method)
	#profmat <- profMat(nc, step = 1.0, method = "binlin")#AP: edited
	
	#Warning: delete first row that shows 0 intensity
	profmat <- profmat[-1,]
	
	rm(nc) #release memory
	tic_length <- length(tic$scan)
	len_1d <- floor(as.numeric(sampling_rate) * as.numeric(mod_time))
	len_2d <- floor(tic_length / len_1d)
	#len_1d <- ceiling(as.numeric(sampling_rate) * as.numeric(mod_time))
	#len_2d <- ceiling(tic_length / len_1d)
	#cat(as.numeric(sampling_rate),mod_time,'\n')
	#len_1d <- floor(as.numeric(sampling_rate) * as.numeric(mod_time))
	#len_2d <- floor(tic_length / len_1d)
	#write.csv(tic, file="test1.csv")
	tic2D = matrix( as.integer(tic$intensity), nrow = len_1d, ncol = len_2d)
	rm(tic) #release memory

	cat('# App-GC : file path : ',filename,'\n')
	cat('# App-GC : Retention time start at (sec) : ',RTstart,'\n')
	cat('# App-GC : Retention time stop at (sec) : ',RTstop,'\n')
	cat('# App-GC : minimum mz : ',MZmin,'\n')
	cat('# App-GC : maximum mz : ',MZmax,'\n')
	cat('# App-GC : total scan : ',tic_length,'\n')
	cat('# App-GC : sampling rate : ',sampling_rate,'\n')

	if (exportSIC==T){

	t_profmat <- t(profmat)
	rm(profmat) #release memory
	c <- 1
	cat('# App-GC : construct SIC file ...\n')
	sic_file = paste(temp_path,"/sicdata_",c,".txt",sep="")
	file.create(sic_file)
	for (r in 1:nrow(t_profmat)){
		if ( (r%%len_1d) == 0 ){
			sic_file = paste(temp_path,"/sicdata_",c,".txt",sep="")
			file.create(sic_file)
			c <- (c+1)
		}
		write.table(t(as.integer(t_profmat[r,])),sic_file,sep = " "
			, row.names = FALSE, col.names = FALSE, append = TRUE,quote = FALSE)

	}
	}

	cat('# App-GC : construct TIC file ...\n')
	if (exportTIC==T){
		tic_file <- paste(temp_path,"/tic_data.txt",sep='')
		write.table(tic2D,tic_file,sep = " ", row.names = FALSE, col.names = FALSE)
	}
	#release memory
	rm(tic2D)
	rm(t_profmat)
	gc()

	save( RTstart,RTstop,MZmin,MZmax,tic_length,sampling_rate, file = paste(temp_path,"/save_parameters.RData",sep="") )
	return(c(RTstart=RTstart,RTstop=RTstop,MZmin=MZmin,MZmax=MZmax,tic_length=tic_length,sampling_rate=sampling_rate))

}else{
	cat('# App-GC : folder ' ,filename, 'exists, the TIC and SIC files in this folder will be used.','\n')
	cat('# App-GC : please remove this folder ', filename, ', if you change parameters (sampling number, modulation time).', '\n')
	load(paste(temp_path,"/save_parameters.RData",sep=""))

	return(c(RTstart=RTstart,RTstop=RTstop,MZmin=MZmin,MZmax=MZmax,tic_length=tic_length,sampling_rate=sampling_rate))
	#return(c(RTstart,RTstop,MZmin,MZmax,tic_length,sampling_rate))
	
}


}

# export detected peaks to msp file
# input : rdata from msPeak method
export.msp.RTcorrection <- function(
		rdata
		,list_peak_id
		,output_file='test'
		,pmin1=450	 # start of 1st RT in second unit
		,minute=F	#time unit of 1st RT in minute
		,export_peak=0
		,MZmin=1
		,n_top_peak=NULL ){	#0: export peaks according to the provided peak list, 1: export peaks before merging from rdata, 2:export peaks after merging from rdata
	#rdata e.g. "msPeak.output.cv0.95.odds10.m2.selected.rda"
	load(rdata, ex <- new.env())
	#ls.str(ex)
	ptab_data = ex$ptab.odd

	if (export_peak==0){
		list_peak_id <- list_peak_id
	}else if (export_peak==1){
		list_peak_id <- as.integer(row.names(ex$ptab1)) 
	}else{  list_peak_id <- as.integer(row.names(ex$ptab2)) 
	}

	if (minute==T) {
		time_unit='min'
	}else { time_unit='sec'}


	#raw data export
	raw_data = paste(output_file,'_peak.csv',sep='')

	csv_header = paste('peak_ID','RT', 'RT2', 'mass', 'height', 'area', 'spectra',sep='\t')
	write.table( csv_header , raw_data ,sep = "", row.names = FALSE, col.names = FALSE, append = FALSE,quote = FALSE) 
	
	#msp format export
	msp = paste(output_file,'.msp',sep='')
	file.create(msp)
	total_peak = length(list_peak_id)
	for (i in list_peak_id){
		#write msp file
		#get profil of each identified peak
		row_peak = ptab_data[ptab_data$id==i,]
		rt1_correct=row_peak$t1

		if (export_peak==1){
			peak_area <- row_peak$area
		}else{  peak_area <- row_peak$parea
		}
		
		if (pmin1>0) {
			rt1_correct=rt1_correct + pmin1
		}
		if (minute==T) {
			rt1_correct=(rt1_correct/60)
			time_unit='min'
		}
		rt1_final=paste(time_unit,rt1_correct,sep='')

		body = ''
		body_csv = ''
		spec_body = ''
		
		if (export_peak==1){
			mass = row_peak[7:length(row_peak)]
		}else{
			mass = row_peak[9:length(row_peak)]
		}
		
		start_MZ=MZmin
		header = sprintf("Name: peak_id_%s\nSynon: DATE: %s\nSynon: AREA: %s\nSynon: RT: %s\nSynon: RT2: sec%s\nComments: msPeakG method\nNum Peaks: %s",row_peak$id, date(), peak_area, rt1_final, row_peak$t2, length(mass))
				
		#max_peak_height=0
		#max_mz_peak_height=0
		#for (j in 1:length(mass) ){
		#	if ((j>1) & (j%%6)==0 ){
		#		body = paste(body,'\n',sep='')
		#	}
		#	#observe peak height
		#	if ( as.numeric(mass[j]) > max_peak_height ){
		#		max_peak_height = as.numeric(mass[j])
		#		max_mz_peak_height = start_MZ
		#	}
		#	body = paste(body,start_MZ,' ',mass[j],'; ',sep='')
		#	spec_body = paste(spec_body,start_MZ,',',mass[j],' ',sep='')
		#	body_csv = paste(body_csv,mass[j],sep=' ')
		#	start_MZ=start_MZ+1
		#}
		
		#######
		collect_mass = data.frame(mz=c(start_MZ),height=c(as.numeric(mass[1])))
		for (j in 2:length(mass) ){
			start_MZ=start_MZ+1
			collect_mass[nrow(collect_mass)+1,] = c(start_MZ,mass[j])			
		}	
		
		collect_mass = collect_mass[order(collect_mass$height,decreasing=TRUE),]
		if( !is.null(n_top_peak)){			
			collect_mass = collect_mass[1:as.numeric(n_top_peak),]
		}
		
		for (j in 1:nrow(collect_mass) ){
			if ((j>1) & (j%%6)==0 ){
				body = paste(body,'\n',sep='')
			}
			body = paste(body,collect_mass[j,'mz'],' ',collect_mass[j,'height'],'; ',sep='')
			spec_body = paste(spec_body,collect_mass[j,'mz'],',',collect_mass[j,'height'],' ',sep='')
			body_csv = paste(body_csv,collect_mass[j,'height'],sep=' ')
		}		
		#######
		
		txt = paste(header,'\n',body,'\n',sep='')
		write.table( txt , msp ,sep = " "
			, row.names = FALSE, col.names = FALSE, append = TRUE,quote = FALSE) 
	
		#write csv file (use in MSPA peak alignment)
		csv_write = paste( paste("peak_id_",row_peak$id,sep=''), rt1_correct, row_peak$t2, collect_mass[1,'mz'], collect_mass[1,'height'], peak_area, body_csv, sep='\t')
		write.table( csv_write , raw_data ,sep = ""
			, row.names = FALSE, col.names = FALSE, append = TRUE,quote = FALSE)
			
		#ali_file is used for erah alignment
		if ((exists('ali_file') && is.data.frame(get('ali_file')) ) == FALSE){
			ali_file = data.frame( 
						"peak_ID"=paste("peak_id_",row_peak$id,sep=''),
						"RT" = rt1_correct, 
						"RT2"= row_peak$t2,
						"mass"= collect_mass[1,'mz'],
						"height"=collect_mass[1,'height'],
						"area"=peak_area,
						"spectra"=spec_body
						)
		}else{
			ali_file = rbind(ali_file, data.frame(
						"peak_ID"=paste("peak_id_",row_peak$id,sep=''),
						"RT" = rt1_correct, 
						"RT2"= row_peak$t2,
						"mass"= collect_mass[1,'mz'],
						"height"=collect_mass[1,'height'],
						"area"=peak_area,
						"spectra"=spec_body
			))
		}			
		
			 
	}

	#ocsv <- paste(output_file, "_peak", ".csv", sep = "")
	#write.csv(ali_file,ocsv,row.names=FALSE,quote = FALSE)
	#release memory
	rm(ptab_data)
	rownames(ali_file)=NULL
	ali_file
}

# export detected peaks to msp file
# input : rdata from msPeak method
export.msp <- function(rdata,list_peak_id,output_file='test',MZmin=1 ){
	#rdata e.g. "msPeak.output.cv0.95.odds10.m2.selected.rda"
	load(rdata, ex <- new.env())
	#ls.str(ex)
	ptab_data = ex$ptab.odd

	#raw data export
	raw_data = paste(output_file,'.csv',sep='')
	#write.table( t(colnames(ptab_data)) , raw_data ,sep = ","
		#, row.names = FALSE, col.names = FALSE, append = FALSE,quote = FALSE) 
	csv_header = paste('Name','1st Dimension Time (s)','2nd Dimension Time (s)','Area','Spectra',sep=',')
	write.table( csv_header , raw_data ,sep = ""
			, row.names = FALSE, col.names = FALSE, append = FALSE,quote = FALSE) 
	
	#msp format export
	msp = paste(output_file,'.msp',sep='')
	file.create(msp)
	total_peak = length(list_peak_id)
	for (i in list_peak_id){
		#write msp file
		#get profil of each identified peak
		row_peak = ptab_data[ptab_data$id==i,]

		mass = row_peak[9:length(row_peak)]
		header = sprintf("Name: peak_id_%s\nSynon: DATE: %s\nSynon: AREA: %s\nSynon: RT: sec%s\nSynon: RT2: sec%s\nComments: msPeakG method\nNum Peaks: %s",row_peak$id, date(), row_peak$area, row_peak$t1, row_peak$t2, length(mass))

		body = ''
		body_csv = ''
		start_MZ=MZmin
		for (j in 1:length(mass) ){
			if ((j>1) & (j%%6)==0 ){
				body = paste(body,'\n',sep='')
			}
			body = paste(body,start_MZ,' ',mass[j],'; ',sep='')
			#body_csv = paste(body_csv,j,':',mass[j],' ',sep='')
			body_csv = paste(body_csv,mass[j],sep=' ')
			start_MZ=start_MZ+1
		}
		txt = paste(header,'\n',body,'\n',sep='')
		write.table( txt , msp ,sep = " "
			, row.names = FALSE, col.names = FALSE, append = TRUE,quote = FALSE) 
	
		#write csv file (use in peak alignment)
		csv_write = paste( paste("peak_id_",row_peak$id,sep=''), row_peak$t1, row_peak$t2, row_peak$area, body_csv, sep=',')
		write.table( csv_write , raw_data ,sep = ""
			, row.names = FALSE, col.names = FALSE, append = TRUE,quote = FALSE) 
	}
	#release memory
	rm(ptab_data)
}

#require msPeakG package
find.tic <- function(
                                tic.fname = "tic_data.txt" # TIC file name
                                ,col="grey"
                                ,p1=6 # modulation time of 1st RT
                                ,p2=0.005 # modulation time of 2nd RT
                                ,pmin1=0 # minimum of 1st RT
                                ,pmin2=0 # minimum of 2nd RT
                                ,plot=F
				,vlevel=15
){
        atic = t(as.matrix(read.table(tic.fname)))
        tic.size = dim(atic)
	#####

        #AS: resolved
	xr = (c(range(1:tic.size[1]))-1)*p1 
	yr = (c(range(1:tic.size[2]))-1)*p2
	xr[1] = xr[1]+pmin1
	yr[1] = yr[1]+pmin2

        if(plot){
                par(mfrow=c(2,1))
                sub.contourG(tic.fname=tic.fname,1:tic.size[1],1:tic.size[2],col=col,p1=p1,p2=p2,pmin1=pmin1,pmin2=pmin2,vlevel=vlevel)
                title("Entire Chromatogram")

        }

cat('# msPeakG :use full 2D TIC as a query with size >> ',tic.size[1],'x',tic.size[2],'\n')
list(idxt1=1:tic.size[1],idxt2=1:tic.size[2],range.t1=xr,range.t2=yr,tic.fname=tic.fname,p1=p1,p2=p2,pmin1=pmin1,pmin2=pmin2)

}

# import CDF data by pattern
cdf.import <- function(path=".",pattern="*.cdf"){
	cat("Pattern: ",pattern,"\n")
	fname = dir(path=path,pattern=pattern)
	cat('list of file input : ',fname,'\n')
	lapply(1:length(fname),function(i) gsub('.cdf','',fname[i]))
}




