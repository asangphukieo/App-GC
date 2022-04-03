

# Load package
library(xcms)
library(ggplot2)
library(plotly)
source('generatePlot.R')
source('GCpipe/auto_GCxGC.R')

#Rscript Run_generatePlot.R Input_path=/home/corgi/Documents/Nong/TEST_DATA_paper/NONG_AUM_GCxGCproject/GCxGCMS/2D_AminoAcidSTD modulation_time=5

args <- commandArgs(trailingOnly = TRUE)
Input_path <- "/home/corgi/Documents/Nong/TEST_DATA_paper/NONG_AUM_GCxGCproject/GCxGCMS/2D_AminoAcidSTD"
Input_file <- "*.cdf"
reference <- "FALSE"
sampling_rate <- FALSE
modulation_time <- 5
out <- ""
profMat_method='binlin'
profMat_step=1.0
rm_prev_tic=FALSE
ref.dir='./'


for(i in 1:length(args)){
  temp<-strsplit(args[i],"=")
  temp<-temp[[1]]
  temp1<-temp[1]
  temp2<-temp[2]
  assign(as.character(temp1),temp2)
}

generatePlot(Input_path=Input_path,
	Input_file=Input_file,
	reference=reference,
	sampling_rate=sampling_rate,
	modulation_time=modulation_time,
	html=TRUE,
	profMat_method=profMat_method,
	profMat_step=profMat_step,
	rm_prev_tic=rm_prev_tic,
	ref.dir=ref.dir)
