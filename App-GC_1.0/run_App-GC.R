##################################################
### App-GC 1.0 : Pipline for data processing of GC-MS and GCxGC-MS data  ###
##################################################

########################################
# Author: Apiwat Sangphukieo           #  
# Email: apiwat.sang@cmu.ac.th         #
########################################

#### Install library #####
# 1. conda install -c r r-base=3.6.1
# 2. conda install -c bioconda bioconductor-metams
# 3. pip uninstall netCDF4 #in case libnetcdf error
# 4. conda install -c conda-forge "libnetcdf=4.6.2"
# 5. conda install -c r r-shiny
# 6. conda install -c conda-forge r-shinythemes
# 7. conda install -c conda-forge r-shinyfiles
# 8. conda install -c conda-forge r-shinyjs
# 9. conda install -c conda-forge r-plotly
# 10. conda install -c conda-forge r-htmlwidgets
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# 3. BiocManager::install("xcms")
# 4. BiocManager::install("CAMERA")
# 5. BiocManager::install("GCalignR")
# 6. BiocManager::install("filesstrings")
# 7. BiocManager::install("msm")
# 8. BiocManager::install("pracma")
# 9. BiocManager::install("doSNOW")
# 10. BiocManager::install("erah")

##########################

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  config_file='./running_App-GC_config.yml'
  #config_file='./parameter_config.yml'
  #stop("At least one argument must be supplied (input file).n", call.=FALSE)
}else{
  config_file=args[1]
}

#test_gc2dgui_config.yml
########################

#1D GC
library(GCalignR) 
library(metaMS)
library(filesstrings)
library(erah)

#2D and 1D GC
library(data.table)
library(yaml)
#conda install -c conda-forge r-filesstrings

call.dir <- getwd()
script.basedir <- paste( call.dir,"App-GC/",sep="/")

# Load packages
source(paste(script.basedir,"msPeakG.autoGCxGC.R",sep="/")) # 2) optional for msPeakG peak detection
source(paste(script.basedir, "mspa/mspa-run.R",sep="/")) # mSPA package for peak alignment
source(paste(script.basedir, "App-GC_GCxGC.R",sep="/")) # 
source(paste(script.basedir, "App-GC_GC.R",sep="/")) # 
source(paste(script.basedir, "erah_AP.R",sep="/")) # import erah version modified by AP

# App-GC : 
#######################################################################
#######################################################################
#config = yaml.load_file("parameter_config.yml")
config = yaml.load_file(config_file)

#######################################################################
# Load GC parameters from config file
#######################################################################
cat(paste('# App-GC : Start pipeline ... ',Sys.time(),'\n',sep=''))
cat('# App-GC : Read prarameters from file ... ',config_file,'\n')

if (config$mode=="1D"){
cat('# App-GC : 1D mode\n')


#######################################################################
#load files by fiding pattern
#######################################################################

cdffiles=list.files(path = config$Input_path, pattern= config$Input_file ,full.names=TRUE, ignore.case=TRUE)
print(cdffiles)
#cdffiles=c("FAME_CAL_4_27Jan20_Run28Jan20.cdf","FAME_CAL_4_27Jan20_Run28Jan20_2.cdf")

#######################################################################
#peak detecttion
#######################################################################
if (config$GC$Peak_detection$do_peak_detection == TRUE){

dir.create(config$Output_folder) # create output folder
cat(paste("# App-GC : peak detecttion - ",config$GC$Peak_detection$method," method \n",sep=""))

#Make peak detection ouput folder
peak_detection_folder1D=paste(config$Output_folder,"/PEAK_DETECTION",sep="")
dir.create(peak_detection_folder1D) 
	
if(config$GC$Peak_detection$method == "centWave" || config$GC$Peak_detection$method == "matchedFilter"){
	###########################################################
	#     peak detecttion - matchedFilter or centWave method  #
	###########################################################	
	
	####### Parameter setting #########
	#minintens=0.1 is in range 0 to 1 
	#centWave
	if(config$GC$Peak_detection$method == "centWave"){
	NONG.GC=metaMSsettings(protocolName="NONG",chrom="GC",
		PeakPicking=c(list(method= "centWave",
		ppm= config$GC$Peak_detection$centWaveParam$ppm,
		peakwidth=c(config$GC$Peak_detection$centWaveParam$peakwidth_min, config$GC$Peak_detection$centWaveParam$peakwidth_max),
		snthresh=config$GC$Peak_detection$centWaveParam$snthresh,
		BPPARAM = SnowParam(workers=config$GC$Peak_detection$worker) )))

	} else if(config$GC$Peak_detection$method == "matchedFilter"){
	#matchedFilter
	NONG.GC=metaMSsettings(protocolName="NONG",chrom="GC",
		PeakPicking=c(list(method = 'matchedFilter',
		fwhm = config$GC$Peak_detection$matchedFilterParam$fwhm, 
		step = config$GC$Peak_detection$matchedFilterParam$step, 
		steps = config$GC$Peak_detection$matchedFilterParam$steps, 
		mzdiff = config$GC$Peak_detection$matchedFilterParam$mzdiff, 
		max = config$GC$Peak_detection$matchedFilterParam$max,
		snthresh = config$GC$Peak_detection$matchedFilterParam$snthresh,
		BPPARAM = SnowParam(workers=config$GC$Peak_detection$worker) )))
	}
	#snthr = config$GC$Peak_detection$matchedFilterParam$snthr, 

	metaSetting(NONG.GC,"DBconstruction")<- list(minfeat=5,
		rttol=0.1,
		minintens=0.0,
		intensityMeasure="maxo",
		DBthreshold=0.8 )

	metaSetting(NONG.GC, "CAMERA")<- list(perfwhm=config$GC$Peak_detection$CAMERA_perfwhm)

	cat("parameters\n")
	metaSetting(NONG.GC, "PeakPicking")
	metaSetting(NONG.GC, "DBconstruction")
	metaSetting(NONG.GC, "CAMERA")	
	

	####### Run peak detection #########
	GCset <- peakDetection(cdffiles,settings=metaSetting(NONG.GC, "PeakPicking") ,convert2list = TRUE)
	allSamples <- lapply(GCset, runCAMERA, chrom = "GC",settings = metaSetting(NONG.GC, "CAMERA"))

	mypeak=list()
	
	#ignore
	#check start scan time (from file 1) for converting RT to 2D RT
	#if (config$GC$Peak_detection$to2D==TRUE){
	#	msFile <- openMSfile(cdffiles[1], backend = "netCDF")
	#	Fileinfo <- runInfo(msFile)
	#	startTime <- Fileinfo$dStartTime
	#}else{ 
	#	startTime = NA }
	
	for(f in 1:length(allSamples)){
		mypeak = append(mypeak, to.file_msp(allSamples[[f]], 
			file 	= names(allSamples)[f] ,
			settings 	= metaSetting(NONG.GC, "DBconstruction"),
			mod_time	= config$GC$Peak_detection$modulation_time,
			to2D		=config$GC$Peak_detection$to2D ,
			comment	=paste("peak detection method", config$GC$Peak_detection$method, sep='') ))
			
		file.copy(paste(names(allSamples)[f], ".msp", sep = ""), peak_detection_folder1D)
		file.copy(paste(names(allSamples)[f], "_peak", ".csv", sep = ""), peak_detection_folder1D)
		file.remove(paste(names(allSamples)[f], ".msp", sep = ""))	
		file.remove(paste(names(allSamples)[f], "_peak", ".csv", sep = ""))	
	}	
	save(mypeak,file='save_peak.Rdata')
	file.copy('save_peak.Rdata', peak_detection_folder1D)
	file.remove('save_peak.Rdata')	
	#sapply(allSamples.msp, length)
	
} else if (config$GC$Peak_detection$method == "ERah") {
	########################################
	#     peak detecttion - Erah method    #
	########################################
	#Make peak detection ouput folder
	peak_detection_folder1D=paste(config$Output_folder,"/PEAK_DETECTION",sep="")
	dir.create(peak_detection_folder1D) 
	peak_detection_folder1D.erah = paste(config$Output_folder,"/PEAK_DETECTION/TEMP_ERAH",sep="")
	peak_detection_folder1D.erah.cdf.exp1 = paste(config$Output_folder,"/PEAK_DETECTION/TEMP_ERAH/EXP1",sep="")

	dir.create(peak_detection_folder1D.erah) 
	dir.create(peak_detection_folder1D.erah.cdf.exp1) 
	cat(paste( "# Erah : copying cdf files to erah temp folder \n" ,sep="" ))
	#file.copy(cdffiles,peak_detection_folder1D.erah.cdf.exp1)
	file.symlink(cdffiles,peak_detection_folder1D.erah.cdf.exp1)

	path_to_exp=peak_detection_folder1D
	working_folder="/TEMP_ERAH"

	folder_cdf=paste(path_to_exp,working_folder,sep="")

	createdt(folder_cdf)

	ex <- newExp(instrumental =paste(folder_cdf,working_folder,"_inst.csv",sep=""),
		phenotype =paste(folder_cdf,working_folder,"_pheno.csv",sep="") )
		
	####### Erah parameter setting #########
	ex.dec.par <- setDecPar( min.peak.width = as.numeric(config$GC$Peak_detection$ERahParam$min.peak.width),  
		min.peak.height = as.numeric(config$GC$Peak_detection$ERahParam$min.peak.height) ,
		noise.threshold = as.numeric(config$GC$Peak_detection$ERahParam$noise.threshold )) 

	####### Erah deconvolution run #########
	ex<- deconvolveComp(ex, ex.dec.par)
	mypeak=list()	
	mypeak = append(mypeak,erah.to.file_msp( ex, output_path=peak_detection_folder1D ))
	save(mypeak,file=paste(peak_detection_folder1D,'/save_peak.Rdata',sep=''))
	#save(ex,file=paste(peak_detection_folder1D,'/save_peak_erah.Rdata',sep=''))	
	#file.move('save_peak.Rdata', peak_detection_folder1D)
	
}
}


#######################################################################
##### peak alignment #####
#######################################################################
if (config$GC$Peak_alignment$do_peak_alignment == TRUE){

load(paste(paste(config$Output_folder,"/PEAK_DETECTION",sep=""),'/save_peak.Rdata',sep=''))
peak_data = mypeak

pdf('Num_peaks.pdf')
passport <- check_input(data = peak_data,plot = T)
dev.off()

if(nrow(passport) > 1){

############# by GCalignR  ##############
if(config$GC$Peak_alignment$method == 'GCalignR'){
	peak_data_aligned <- align_chromatograms(data = peak_data, # input data
	    rt_col_name = "RT", # retention time variable name 
	    rt_cutoff_low = 0, # remove peaks below 15 Minutes
	    rt_cutoff_high = 100000000000000, # remove peaks exceeding 45 Minutes
	    reference = config$GC$Peak_alignment$reference, # NULL to choose automatically 
	    max_linear_shift = config$GC$Peak_alignment$GCalignRParam$max_linear_shift, # max. shift for linear corrections
	    max_diff_peak2mean = config$GC$Peak_alignment$GCalignRParam$max_diff_peak2mean, # max. distance of a peak to the mean across samples
	    min_diff_peak2peak = config$GC$Peak_alignment$GCalignRParam$min_diff_peak2peak, # min. expected distance between peaks 
	    blanks = NULL, # negative control
	    delete_single_peak = config$GC$Peak_alignment$GCalignRParam$delete_single_peak, # delete peaks that are present in just one sample 
	    write_output = NULL) # add variable names to write aligned data to text files

	#peak_data_aligned$aligned$RT # to access the aligned retention times
	#peak_data_aligned$aligned$hight # to access the aligned area data
	#peak_data_aligned$aligned$peak_num

	RT_aligned = peak_data_aligned$aligned$RT[order(peak_data_aligned$aligned$RT$mean_RT),]
	peak_num_aligned = peak_data_aligned$aligned$peak_ID[order(peak_data_aligned$aligned$RT$mean_RT),]
	peak_num_aligned[peak_num_aligned==0]<-NA


}else if(config$GC$Peak_alignment$method == 'ERah'){
	############# by ERah  ##############	
	ex.al.par<- setAlPar(min.spectra.cor = config$GC$Peak_alignment$ERahParam$min.spectra.cor,
				max.time.dist = config$GC$Peak_alignment$ERahParam$max.time.dist ,
				mz.range =40:600) #AP: need to edit
	ex.mypeak <- Rdata2MetaboSet(peak_data,from2D=FALSE)			
	ex.align <- alignComp(ex.mypeak, alParameters = ex.al.par)
	peak_num_aligned <- erah.export.alignList(ex.align, export="ID", from2D=FALSE)
	RT_aligned <- erah.export.alignList(ex.align, export="RT", from2D=FALSE)
}

write.csv(peak_num_aligned,"peak_num_align.csv", row.names = FALSE, quote = FALSE)

write.csv(RT_aligned,"peak_RT_align.csv", row.names = FALSE, quote = FALSE)
ali_folder_1D=paste(config$Output_folder,'/PEAK_ALIGNMENT',sep='')
dir.create(ali_folder_1D)

#save aligned peaks
#>> mean_RT Run1 Run2 Run3
save(peak_num_aligned,file='save_ali_peak.Rdata')
file.copy('save_ali_peak.Rdata', ali_folder_1D)
file.copy('peak_RT_align.csv', ali_folder_1D)
file.copy('peak_num_align.csv', ali_folder_1D)
file.copy('Num_peaks.pdf', ali_folder_1D)

file.remove('save_ali_peak.Rdata')
file.remove('peak_RT_align.csv')
file.remove('peak_num_align.csv')
file.remove('Num_peaks.pdf')

}else {cat("# Peak alignment : need more than 1 input file for performing alignment!\n")} #end check peak alignment input (passport)
} #end peak alignment

} #end GC method



#######################################################################
# 2D GC 
#######################################################################
if (config$mode=="2D"){
cat('# App-GC : 2D mode\n')
#######################################################################
#load files by fiding pattern
#######################################################################

cdffiles=list.files(path = config$Input_path, pattern= config$Input_file ,full.names=FALSE, ignore.case=TRUE)
print(cdffiles)
#######################################################################
dir.create(config$Output_folder) # create output folder

#######################################################################
# Load 2D GC parameters from config file
#######################################################################

##### selected region #####
## step 1 : define parameters for msPeak (peak detection algoritm) , ##

##### selected region #####
## step 1 : detect parameter from sample file and define parameters for msPeak (peak detection algoritm) ##
if ( !is.null(config$GC2D$Peak_detection$reference) ){
	ref_file=config$GC2D$Peak_detection$reference
}else{
	ref_file=cdffiles[1]
}
cat('# App-GC : extract information from reference file  : ',ref_file,'\n')
#cat(ref_file,'\n')

ref_file_path=paste(config$Input_path,'/',ref_file,sep='')
ref_folder_name=gsub('.cdf$','',ref_file)
ref_folder_name <- paste(config$Output_folder,'/',ref_folder_name,sep='') #

if(is.na(config$GC2D$Peak_detection$sampling_rate) || is.null(config$GC2D$Peak_detection$sampling_rate) || config$GC2D$Peak_detection$sampling_rate == "NA" ){	#from people in the lab, but it can be automaticcally calculated
	sampling_rate = FALSE
}else {
	sampling_rate = config$GC2D$Peak_detection$sampling_rate
}
modulation_time = config$GC2D$Peak_detection$modulation_time #from people in the lab

cat('# App-GC : generate TIC of reference file \n')
detected_param = load_2DGCdata(filename=ref_file_path, 
	mod_time=modulation_time, 
	sampling_rate=sampling_rate, 
	exportTIC=T, 
	exportSIC=T, 
	temp_path=ref_folder_name,
	profMat_method=config$GC2D$Peak_detection$profMat_method,
	profMat_step=config$GC2D$Peak_detection$profMat_step,
	rm_prev_tic=config$GC2D$Peak_detection$rm_prev_tic)
sampling_rate = detected_param['sampling_rate']	#from people in the lab, but it can be automaticcally calculated in load_2DGCdata method
RTstart=detected_param['RTstart']
RTstop=detected_param['RTstop']


full_scan=config$GC2D$Peak_detection$full_scan
if (full_scan == TRUE){
	config$GC2D$Peak_detection$useXY=FALSE #AS
	cat('# App-GC : Runing full scan mode is not generally recommended. \n')
	cat('# App-GC : You might consider to set time_scan to increase computation process. \n')
	RT1Dstart= (RTstart/60)	# in minute unit
	RT1Dstop= (RTstop/60)
	RT2Dstart= 0.1	#in range of 1 to modulation time , in second
	RT2Dstop= modulation_time	#in range of 1 to modulation time , in second
}else{
	if (config$GC2D$Peak_detection$useXY != TRUE){
	RT1Dstart= config$GC2D$Peak_detection$time_scan$RT1Dstart	# in minute unit
	RT1Dstop= config$GC2D$Peak_detection$time_scan$RT1Dstop
	RT2Dstart= config$GC2D$Peak_detection$time_scan$RT2Dstart	#in range of 1 to modulation time , in second
	RT2Dstop= config$GC2D$Peak_detection$time_scan$RT2Dstop	#in range of 1 to modulation time , in second
	}
}

#######################################################################
# peak detection
#######################################################################
PEAK_2D_DETECT_FOLDER=paste(config$Output_folder,'/PEAK_DETECTION',sep='')

if(config$GC2D$Peak_detection$do_peak_detection==T){
	### for loop to run each CDF file ### 	
	dir.create(PEAK_2D_DETECT_FOLDER)
	
mypeak=list()
for (file in cdffiles){
	cat('# App-GC : running file input ...' ,file,'\n')
	##### load single GCxGC-MS file in CDF format and write out the tic file and sic files in temp_path
	folder_name=gsub('.cdf$','',file)
	temp_path <- paste(config$Output_folder,'/',folder_name,sep='') #make folder

	## step 2 : loading data and extract total ion chromatogram and single ion chromatogram##
	
	detected_param = load_2DGCdata(filename=paste(config$Input_path,'/',file,sep=''), 
		mod_time=modulation_time, 
		sampling_rate=sampling_rate, 
		exportTIC=T, 
		exportSIC=T, 
		temp_path=temp_path,
		profMat_method=config$GC2D$Peak_detection$profMat_method,
		profMat_step=config$GC2D$Peak_detection$profMat_step,
		rm_prev_tic=config$GC2D$Peak_detection$rm_prev_tic)

	RTstart=detected_param['RTstart']
	RTstop=detected_param['RTstop']
	MZmin=detected_param['MZmin']

	
	######## step 3 : convert time to xy plot ###################
	cat('# App-GC : detected sampling rate ' ,sampling_rate,'\n')
	if (config$GC2D$Peak_detection$useXY != TRUE){
	#to convert time (minute) to xy region to select for peak detection
		x_start= (((RT1Dstart*60)-RTstart)/modulation_time )+1
		x_stop= ((RT1Dstop*60)-RTstart)/modulation_time
		y_start= (RT2Dstart*sampling_rate)
		y_stop= (RT2Dstop*sampling_rate)
	}else{
		#use xy region from user input to select for peak detection
		x_start= config$GC2D$Peak_detection$x1	# in minute unit
		x_stop= config$GC2D$Peak_detection$x2
		y_start= config$GC2D$Peak_detection$y1	#in range of 1 to modulation time , in second
		y_stop= config$GC2D$Peak_detection$y2	#in range of 1 to modulation time , in second
	}


	## step 4 : define region for msPeak (peak detection algoritm) ##
	cat('# App-GC : User-defined X start',x_start,'\n')
	cat('# App-GC : User-defined X stop',x_stop,'\n')
	cat('# App-GC : User-defined Y start',y_start,'\n')
	cat('# App-GC : User-defined Y stop',y_stop,'\n')
	sdg = find.subG2(subrt1=x_start:x_stop,
		subrt2= y_start:y_stop,
		tic.fname= paste(temp_path,"/tic_data.txt",sep=''),
		plot=T,
		p1= modulation_time,
		p2= (1/sampling_rate),
		pmin1=RTstart,
		pmin2=0,
		vlevel=50)

	#### step 5 : peak detection , the script is msPeakG.autoGCxGC.R #############
	setwd(temp_path)
	cat('# App-GC : running peak detection with msPeakG with ' ,config$GC2D$Peak_detection$method,' method \n')
	if(config$GC2D$Peak_detection$method == "NEB"){
		ngb_method=FALSE
	}else{
		ngb_method=TRUE
	}
	cat('# App-GC : running peak detection with multiple cpu processing ...' ,file,'\n')
	peaks <- msPeakG.multicore(sdg, 
		zcv=config$GC2D$Peak_detection$zcv,
		method=config$GC2D$Peak_detection$method_parallel,
		cpu=config$GC2D$Peak_detection$worker,
		verbose=T,
		cv=config$GC2D$Peak_detection$cv,
		out.fname="msPeakG",
		sic.name=paste(temp_path,"/sicdata",sep=''),
		tic.fname=paste(temp_path,"/tic_data.txt",sep=''), 
		ngb=ngb_method,
		merge=config$GC2D$Peak_detection$peak_merge,
		den=config$GC2D$Peak_detection$method) #use NEB method or fft
	msPeakG_output <- "msPeakG.rda"
	list_peak_id <- as.integer(row.names(peaks$aptable)) #list peak after peak merging

	####### step 6 :  export peak to msp file ################
	cat('# App-GC : export peak to file ...' ,file,'\n')
	if (config$GC2D$Peak_detection$peak_merge==FALSE){
		#0: export peaks according to the provided peak list, 1: export peaks before merging from rdata, 2:export peaks after merging from rdata
		export_peak=1 
	}else{
		export_peak=0
	}
	peak_out = export.msp.RTcorrection(msPeakG_output, 
		list_peak_id, 
		output_file=file, 
		export_peak=export_peak,
		pmin1=0,
		minute=T,
		MZmin=MZmin,
		n_top_peak=20) #pmin1 is start first dimention retention time for scanning, #list_peak_id = list of selected peaks to export to msp file

	mypeak[[folder_name]] = peak_out #save peak to Rdata
	
	###### step 7 : plot peak output
	png(filename="peaks_plot.png", width = 1500, height = 1200)
	tplotG(peaks,merge=config$GC2D$Peak_detection$peak_merge)
	dev.off()
	setwd('../')
	gc(reset = TRUE)
	
	file.copy(msPeakG_output, PEAK_2D_DETECT_FOLDER)
	file.copy(paste(temp_path,'/',file,'_peak.csv',sep=''), PEAK_2D_DETECT_FOLDER)
	file.copy(paste(temp_path,'/',file,'.msp',sep=''), PEAK_2D_DETECT_FOLDER)
} #end for loop


save(mypeak, file=paste(PEAK_2D_DETECT_FOLDER,'/save_peak.Rdata',sep=''))

} #end peak detection

####################################################################
##### GCxGC peak alignment between samples #####
####################################################################
peak_alignment<-config$GC2D$Peak_alignment$do_peak_alignment
if(peak_alignment==TRUE){

dir.ali.main <- paste(call.dir,'/PEAK_ALIGNMENT',sep='')
dir.create(dir.ali.main)
ali_folder_2D=paste(config$Output_folder,'/PEAK_ALIGNMENT',sep='')
dir.create(ali_folder_2D)

load(paste(PEAK_2D_DETECT_FOLDER,'/save_peak.Rdata',sep=''))
peak_data = mypeak
	
pdf('Num_peaks.pdf')
passport <- check_input(data = peak_data,plot = T)
dev.off()
	
if(config$GC2D$Peak_alignment$method == 'mSPA'){
	dir.ali.input <- paste(dir.ali.main,'/RAW_PEAK',sep='')
	dir.create(dir.ali.input)

	# copy all detected peak files (.csv) into alignment folder
	PEAK_2D_DETECT_FOLDER = paste(config$Output_folder,'/PEAK_DETECTION',sep='')
	if (dir.exists(PEAK_2D_DETECT_FOLDER)){
		for (i in cdffiles ){
			cdf_name = gsub('.cdf','',basename(i))
			file.copy(paste(PEAK_2D_DETECT_FOLDER,'/',i,"_peak.csv",sep=''),dir.ali.input)
		}		
	}else{
		cat('# App-GC : cannot find peak detection output, please rerun peak detection step before peak alignment! \n')
	}

	#mSPA method for homogenous data
	#filePrefixRef: pattern of single reference file
	file_ref=config$GC2D$Peak_alignment$reference
	if ( !is.null(file_ref)){
		file_ref=paste('^',file_ref,sep='')
	}else {
		file_ref=FALSE
	}

	mspa.run(filePrefix='./*'
		,filePrefixRef =file_ref
		,plot = TRUE
		,dataPath = dir.ali.input
		,directory= dir.ali.main
		,variant = config$GC2D$Peak_alignment$mSPAParam$method
		,w = config$GC2D$Peak_alignment$mSPAParam$w
		,k = config$GC2D$Peak_alignment$mSPAParam$k
		,rho = config$GC2D$Peak_alignment$mSPAParam$rho
		,distance = config$GC2D$Peak_alignment$mSPAParam$distance #One of 'maximum', 'manhattan', 'canberra', or 'euclidean'
		,arraySimilarity = config$GC2D$Peak_alignment$mSPAParam$arraySimilarity # similarity (sm): "dot", "pearson"
		,performance = FALSE
		,createReferenceAlignment = config$GC2D$Peak_alignment$mSPAParam$createReferenceAlignment
		,referenceAlignmentVariant = config$GC2D$Peak_alignment$mSPAParam$referenceAlignmentVariant #"Create given reference alignment variant. Variants are 'mSPA', 'mgma' (recommended)
		,sd1Thres = config$GC2D$Peak_alignment$mSPAParam$sd1Thres
		,sd2Thres = config$GC2D$Peak_alignment$mSPAParam$sd2Thres)


	setwd(call.dir)

	# copy output to PEAK_ALIGNMENT folder
	output_ali_folder=paste(dir.ali.main,'/',config$GC2D$Peak_alignment$mSPAParam$method,'/1/',sep='')
	dir.create(paste(config$Output_folder,'/PEAK_ALIGNMENT',sep=''))
	copy_files <- list.files(output_ali_folder, r = T)
	file.copy(paste(output_ali_folder,copy_files,sep=''), ali_folder_2D)
	file.remove(paste(output_ali_folder,copy_files,sep=''))
	unlink(dir.ali.main,recursive=TRUE)
	
} else if(config$GC2D$Peak_alignment$method == 'ERah'){
	############# by ERah  ##############
	ex.al.par<- setAlPar(min.spectra.cor = config$GC2D$Peak_alignment$ERahParam$min.spectra.cor,
				max.time.dist = config$GC2D$Peak_alignment$ERahParam$max.time.dist ,
				mz.range =40:600)
	ex.mypeak <- Rdata2MetaboSet(peak_data,from2D=TRUE)			
	ex.align <- alignComp(ex.mypeak ,alParameters =ex.al.par)
	peak_num_aligned <- erah.export.alignList(ex.align, export="ID", from2D=TRUE)
	RT_aligned <- erah.export.alignList(ex.align, export="RT", from2D=TRUE)

	write.csv(peak_num_aligned,"peak_num_align.csv", row.names = FALSE, quote = FALSE)

	write.csv(RT_aligned,"peak_RT_align.csv", row.names = FALSE, quote = FALSE)

	#save aligned peaks
	#>> mean_RT Run1 Run2 Run3
	save(peak_num_aligned,file='save_ali_peak.Rdata')
	file.copy('save_ali_peak.Rdata', ali_folder_2D)
	file.copy('peak_RT_align.csv', ali_folder_2D)
	file.copy('peak_num_align.csv', ali_folder_2D)
	file.copy('Num_peaks.pdf', ali_folder_2D)
	file.remove('save_ali_peak.Rdata')
	file.remove('peak_RT_align.csv')
	file.remove('peak_num_align.csv')
	file.remove('Num_peaks.pdf')	
}
} #end of peak alignment
} #end 2D GC mode

#######################################################################
##### Peak annotation  #####
#######################################################################
if (config$Peak_annotation$do_peak_annotation == TRUE) {
#annotate: all ".msp" files in OUTPUT/PEAK_DETECTION folder
#Output: ".annot" files in OUTPUT/PEAK_ANNOTATION folder

annot_folder=paste(config$Output_folder,'/PEAK_ANNOTATION',sep='')
dir.create(annot_folder)

list_msp=Sys.glob(paste(config$Output_folder,"/PEAK_DETECTION/*.msp",sep=""))
myannot=list()

#MSP_OUT="FAME_CAL_4_27Jan20_Run28Jan20_p1.msp"
#wine ../2019_02_22_MSPepSearch_x64/MSPepSearch64.exe Q /MinInt 50 /MinMF 700 /HITS 1 /MAIN ../NIST_library/mainlib /REPL ../NIST_library/replib/ /INP ./$MSP_OUT /OUTTAB ./$MSP_OUT.annot
#wine /home/corgi/Documents/Nong/2019_02_22_MSPepSearch_x64/MSPepSearch64.exe Q /MinInt 250 /MinMF 700 /HITS 1 /MAIN /home/corgi/Documents/Nong/NIST_library/mainlib /REPL /home/corgi/Documents/Nong/NIST_library/replib /INP ../PEAK_DETECTION/MixAminoAcid_split100_03082020.cdf.msp /OUTTAB ./MSP_OUT.annot

path_MSPepSearch= config$Peak_annotation$path_MSPepSearch
MinInt= config$Peak_annotation$MinInt
MinMF= config$Peak_annotation$MinMF
HITS= config$Peak_annotation$HITS
path_MAIN= config$Peak_annotation$path_MAIN
path_REPL= config$Peak_annotation$path_REPL

for( input_msp in list_msp){

	cat(" # App-GC : Run annotation ...",input_msp,"\n")
	output_annot=paste(input_msp,".annot",sep="")
	search_command= paste("wine ",path_MSPepSearch," Q /MinInt ",MinInt," /MinMF ",MinMF," /HITS ",HITS," /MAIN ",path_MAIN," /REPL ",path_REPL," /INP ",input_msp," /OUTTAB ",output_annot,sep="")
	cat('# App-GC : excute ... ',search_command,'\n')

	system(search_command, TRUE)

	#combine annotation to csv file (i)
	input_csv = gsub('.msp','_peak.csv',input_msp)

	output_annot_csv = gsub('.msp','_annot.csv',input_msp)

	system(paste('sed -i s/"\'"//g ',output_annot,sep=''), TRUE)

	read_input_csv = read.table(input_csv,header=T,sep='\t')
	#read_output_annot=read.table(output_annot,comment='>', header=T,sep='\t')
	command_read_annot=paste("grep '>' -v ",output_annot," |cut -f1-9",sep="")
	read_output_annot = fread(cmd=command_read_annot)
	
	#combine annotations to csv file
	#cat('Name\t1st Dimension Time (min)\t2nd Dimension Time (sec)\tArea\tMass\tHeight\tProb (%)\tMF\tAnnotation\n',file=output_annot_csv,append=F)

	annot_data=data.frame()
	count_peak=0
	for(x in read_input_csv$peak_ID){
		match_annot = read_output_annot[which(read_output_annot$Unknow==x),][,c(9,8,7)]
		query=read_input_csv[which(read_input_csv$peak_ID==x),][,1:6]
		#print(match_annot)
		
		if(nrow(match_annot)>=1){
			first_RT=round(as.numeric(query[2]), 3)

			sec_RT=round(as.numeric(query[3]), 3)
			area=query[4]	
			mass=query[5]
			height=query[6]
			for(i in 1:nrow(match_annot)){
				
				annotation=match_annot[i,]$Name
				#print(annotation)
				#annotation=gsub("\'","",match_annot[i,]$Name)
				
				prob=match_annot[i,2]
				MF=match_annot[i,3]
				#cat(paste(x,first_RT,sec_RT,area,mass,height,prob,MF,annotation,'\n', sep='\t'), file=output_annot_csv, append=T)
				count_peak=count_peak+1
				if (count_peak == 1){
					
					annot_data = data.frame( 
								"peak_ID"=x,
								"RT" = first_RT, 
								"RT2"= sec_RT,
								"mass"= mass,
								"height"=height,
								"area"=area,
								"prob"=prob,
								"MF"=MF,
								"annotation"=annotation)
				}else{
					annot_data = rbind(annot_data, data.frame(
								"peak_ID"=x,
								"RT" = first_RT, 
								"RT2"= sec_RT,
								"mass"= mass,
								"height"=height,
								"area"=area,
								"prob"=prob,
								"MF"=MF,
								"annotation"=annotation))
				}
								
				
			}
			
		}
	}
	
	write.csv(annot_data,output_annot_csv,row.names=FALSE,quote = FALSE)
	
	basename_file=gsub('.msp','',basename(input_msp))
	basename_file=gsub('.cdf','',basename_file)
	myannot[[basename_file]] = annot_data

	file.copy(output_annot, annot_folder) 
	file.copy(output_annot_csv, annot_folder)
	file.remove(output_annot)
	file.remove(output_annot_csv)
}

		
#save peak annotation
save(myannot,file='save_annotation_peak.Rdata')
file.copy('save_annotation_peak.Rdata', annot_folder)
file.remove('save_annotation_peak.Rdata')		
} #end peak annotation


# Peak summary
if (config$Peak_summary$do_peak_summary == TRUE) {

summary_folder=paste(config$Output_folder,'/SUMMARY',sep='')
dir.create(summary_folder)
output_summary_file = paste(summary_folder,'/peak_summary.txt',sep='')

summary_table = summary_result(app_gc_path=config$Output_folder, ref=config$Peak_summary$reference, name=config$Peak_summary$select_name)
write.table(summary_table,file=output_summary_file,quote = FALSE,sep='\t', row.name=FALSE)

}

file.copy(config_file, config$Output_folder)
file.remove(config_file)

cat(paste('# App-GC : Finish pipeline ... ',Sys.time(),'\n',sep=''))
############################################
############################################
############################################

#Usage
##Rscript run_App-GC.R parameter_config.yml

