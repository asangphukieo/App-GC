#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
    make_option("--filePrefix", default="Standard_",
    help="Prefix of filenames which are looked up in dataPath. [default \"%default\"]"),
    make_option("--dataPath", default="../data/mSPA_Dataset_I",
        help="Path relative to current directory containing input files. [default
        \"%default\"]"),
    make_option("--variant", default="PAD", help="Name of the method variant to use. One of 'PAD','PAS','DW-PAS','SW-PAD','PAM', or 'OP-PAM'. [default
        \"%default\"]"),
    make_option("--w", default=0.5, help="Value of the w parameter. [default
        \"%default\"]"),
    make_option("--k", default=5, help="Value of the k parameter. [default
        \"%default\"]"),
    make_option("--rho", default=0.95, help="Value of the rho parameter. [default
        \"%default\"]"),
    make_option("--distance", default="canberra", help="Value of the distance parameter. One of 'maximum', 'manhattan', 'canberra', or 'euclidean'. [default
        \"%default\"]"),
    make_option("--arraySimilarity", default="dot", help="Value of the similarity parameter. One of 'dot', or 'linCorr'. [default
        \"%default\"]"),
    make_option("--plot", default=FALSE, help="Create diagnostic plots. [default
        \"%default\"]"),
    make_option("--performance", default=TRUE, help="Evaluate performance of pairwise alignments. [default
        \"%default\"]"),
    make_option("--createReferenceAlignment", default=FALSE, help="Create reference alignment table and exit to console. [default
        \"%default\"]"),
    make_option("--referenceAlignmentVariant", default="mSPA", help="Create given reference alignment variant. Variants are 'mSPA', 'mgma' (recommended) [default
        \"%default\"]"),
    make_option("--sd1Thres", default=25.0, help="Create reference alignment with given standard deviation threshold for rt1. [default
        \"%default\"]"),
    make_option("--sd2Thres", default=0.05, help="Create reference alignment with given standard deviation threshold for rt2. [default
        \"%default\"]"),
    make_option("--directory", default="mspa-output",
        help="Base directory for output. [default \"%default\"]")
)
options(error=traceback)
opt <- parse_args(OptionParser(option_list=option_list))
call.dir <- getwd()
script.name <- "mspa-evaluation.R"
script.basedir <- call.dir
initial.options <- commandArgs(trailingOnly = FALSE)
#if we were called as a script from the command line
#--file gives the script to execute
if(length(grep("--file",commandArgs(trailingOnly = FALSE)))>0) {
    cat("We were called as a script and not sourced!\n")
    file.arg.name <- "--file="
    script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    if(length(grep("^/{1}",script.name))==0) {
      cat("Relative path\n")
      script.basedir <- dirname(paste(sep="/",call.dir,script.name))
    }else{
      cat("Absolute path\n")
      script.basedir <- dirname(script.name)
    }
}
mSPA.script <- paste(sep="/", script.basedir, "mSPA.R")
malign.script <- paste(sep="/",script.basedir,"malign.R")
calcAlignments.script <- paste(sep="/",script.basedir, "calcAlignments.R")

#"../data/mSPA_Dataset_I"
dataPath = opt$dataPath
#"mspa-output"
outputPath = opt$directory
#paste("^",opt$filePrefix,"*",sep="")
pattern = opt$filePrefix
#NH create output directory
if(!file.exists(outputPath)) {
    output = dir.create(outputPath,recursive=TRUE)
}

# NH cat data and output path
#cat("Data path: ",dataPath,"\n")
#cat("Output path: ",outputPath,"\n")

# load mSPA package
source(mSPA.script)
source(malign.script)

sd1Thres=opt$sd1Thres
sd2Thres=opt$sd2Thres

cat("sd1Thres=",sd1Thres," sd2Thres=",sd2Thres,"\n")

arraySimilarity <- opt$arraySimilarity
if(arraySimilarity == "linCorr") {
    arraySimilarity = "pearson"
}

# NH added relative path parameter
# load the ChromaToF data
# changed name to data
data = load.data(path=dataPath,pattern=pattern,sd1Thres=sd1Thres,sd2Thres=sd2Thres,plot=opt$plot,outputPath=outputPath,createReferenceAlignment=opt$createReferenceAlignment,variant=opt$referenceAlignmentVariant)

if(opt$createReferenceAlignment){
    if(!is.null(warnings()) && length(warnings())>0){
        warnings()
    }else{
      cat("Wrote reference alignment, exiting!\n")
    }
    quit(save="no",status=0)
}

toolParams <- list(name=opt$variant,w=opt$w,k=opt$k,rho=opt$rho,distance=opt$distance,similarity=arraySimilarity,anal=opt$performance,plot=opt$plot)

parameters <- list(toolParams=toolParams)

# NH factored out common functionality
source(calcAlignments.script)
calcAlignments(outputPath=outputPath,data=data,parameters=parameters)
if(!is.null(warnings()) && length(warnings())>0){
    warnings()
}
