
#modified by Apiwat Sangphukieo from mspa-evaluation 
#make callable function
# retrieve parameter from the input file

call.dir.mspa <- paste(getwd(),'/App-GC/mspa',sep='')
script.mspa <- call.dir.mspa

mSPA.script <- paste(sep="/", script.mspa, "mSPA.R")
malign.script <- paste(sep="/",script.mspa,"malign.R")
calcAlignments.script <- paste(sep="/",script.mspa, "calcAlignments.R")

# load mSPA package
source(mSPA.script)
source(malign.script)
source(calcAlignments.script)

#mspa.run(dataPath = "../../INPUT/mSPA_Dataset_I")
#AP: add filePrefixRef to set a single reference file from the list of input, used for alignment
mspa.run <- function(
		filePrefix = "^Standard_*"
		,filePrefixRef ='filePrefixRef'
		,dataPath = "../data/mSPA_Dataset_I"
		,variant = "PAD"
		,w = 0.5
		,k = 5
		,rho = 0.95
		,distance = "canberra" #One of 'maximum', 'manhattan', 'canberra', or 'euclidean'
		,arraySimilarity = "dot"
		,plot = FALSE
		,performance = FALSE
		,createReferenceAlignment = FALSE
		,referenceAlignmentVariant = "mgma" #"Create given reference alignment variant. Variants are 'mSPA', 'mgma' (recommended)
		,sd1Thres = 25.0
		,sd2Thres = 0.05
		,directory = "mspa-output" #output dir
		) {



	#"mspa-output"
	outputPath = directory
	#paste("^",filePrefix,"*",sep="")
	pattern = filePrefix
	#NH create output directory
	if(!file.exists(outputPath)) {
	    output = dir.create(outputPath,recursive=TRUE)
	}

	cat("sd1Thres=",sd1Thres," sd2Thres=",sd2Thres,"\n")

	if(arraySimilarity == "linCorr") {
	    arraySimilarity = "pearson"
	}

	# NH added relative path parameter
	# load the ChromaToF data
	# changed name to data
	#AP: add filePrefixRef to set refferent file from the list of input
	data = load.data.setref(path=dataPath
		,pattern=pattern
		,filePrefixRef=filePrefixRef
		,sd1Thres=sd1Thres
		,sd2Thres=sd2Thres
		,plot=plot
		,outputPath=outputPath
		,createReferenceAlignment=createReferenceAlignment
		,variant=referenceAlignmentVariant)

	if(createReferenceAlignment){
	    if(!is.null(warnings()) && length(warnings())>0){
		warnings()
	    }else{
	      cat("Wrote reference alignment, exiting!\n")
	    }
	    quit(save="no",status=0)
	}

	toolParams <- list(name=variant,w=w,k=k,rho=rho,distance=distance,similarity=arraySimilarity,anal=performance,plot=plot)
	parameters <- list(toolParams=toolParams)

	# NH factored out common functionality
	print("run alignment......")
	
	calcAlignments(outputPath=outputPath,data=data,parameters=parameters)
	if(!is.null(warnings()) && length(warnings())>0){
	    warnings()
	}

}
