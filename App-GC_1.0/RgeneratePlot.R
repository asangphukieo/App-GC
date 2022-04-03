

# Load package
library(xcms)
library(yaml)
library(ggplot2)

generatePlot <- function(Input_path='',Input_file='',reference=FALSE,sampling_rate=200,modulation_time=5) {
	#to extract RT start, RT end , sampling_rate, MZ min, MZmax from reference file and save to Rdata
	#and 
	cat('# GCxPipe : start 2D-GC data preprocess \n')
	ref.dir <- getwd()
	#config = yaml.load_file(parameter_file)
	cdffiles=list.files(path = Input_path, pattern= Input_file ,full.names=FALSE, ignore.case=TRUE)
	cat('# GCxPipe : all input files in the path \n')
	cat('# GCxPipe : file pattern =',Input_file,'\n')
	cat('# GCxPipe : ',cdffiles,'\n')
	############################

	if ( reference != FALSE){
		ref_file=reference
	}else{
		ref_file=cdffiles[1]
	}
	cat('# GCxPipe : extract information from reference file  : \n')
	cat('# GCxPipe : ',ref_file,'\n')

	nc <- xcmsRaw(paste(Input_path,'/',ref_file,sep='') )
	RTstart=min(nc@scantime)
	RTstop=max(nc@scantime)
	MZmin=min(nc@mzrange)
	MZmax=max(nc@mzrange)

	mzmin = MZmin
	mzmax = MZmax
	sampling_rate = sampling_rate	#from people in the lab, but it can be automaticcally calculated
	modulation_time = modulation_time #from people in the lab
	nc=''
	gc() 

	folder_name=gsub('.cdf$','',ref_file)
	temp_path <- paste(ref.dir,'/',folder_name,sep='') #make folder
	detected_param = load_2DGCdata(filename=paste(Input_path,'/',ref_file,sep=''), 
		mod_time=modulation_time, 
		sam_rate=sampling_rate, 
		exportTIC=T, 
		exportSIC=T, 
		temp_path=temp_path)

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
	ggplotly(p, source="graph") %>% layout(xaxis=list(fixedrange=TRUE)) %>% layout(yaxis=list(fixedrange=TRUE))

}
