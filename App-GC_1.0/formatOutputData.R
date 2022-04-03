formatOutputData <- function(ali_Rdata,annot_Rdata){
  load(ali_Rdata)
  load(annot_Rdata)
  cat("# App-GC : generate alignment heatmap\n")
  cat("# App-GC : ",ali_Rdata,"\n")
  cat("# App-GC : ",annot_Rdata,"\n")
  #load('/home/gcpipe/Documents/GCPIPE/gcgui/OUTPUT/test10/PEAK_ALIGNMENT/save_ali_peak.Rdata')
  #load('/home/gcpipe/Documents/GCPIPE/gcgui/OUTPUT/test10/PEAK_ANNOTATION/save_annotation_peak.Rdata')
  
  format_RT=format(round(as.numeric(peak_num_aligned[,"mean_RT"]),3),nsmall=2)
  
  #if RT2 available
  if("mean_RT2" %in% colnames(peak_num_aligned)){
    mean_RT2=format(round(as.numeric(peak_num_aligned[,"mean_RT2"]),3),nsmall=2)
    
    format_RT=paste(format(round(as.numeric(peak_num_aligned[,"mean_RT"]),3),nsmall=2),mean_RT2,sep=',')
    peak_num_aligned=peak_num_aligned[,colnames(peak_num_aligned)!="mean_RT2"]
    
  }
  
  rownames(peak_num_aligned) = format_RT
  mean_RT = format_RT
  peak_num_aligned = peak_num_aligned[,colnames(peak_num_aligned)!="mean_RT"]
  peak_num_aligned = cbind(peak_num_aligned, mean_RT)
  
  library(dplyr)
  myannot_merged = bind_rows(myannot, .id = "run") #merge annotation list to a dataframe
  myannot_grouped <- myannot_merged %>%
    group_by(run,peak_ID) %>% 
    summarise(annotation=paste(unique(annotation), collapse=" | "), .groups = 'keep') #group annotation by run and peak id
  peak_num_aligned_long = reshape2::melt(peak_num_aligned, id.vars=c("mean_RT")) #convert to long-form dataframe
  colnames(peak_num_aligned_long) = c("mean_RT","run","peakid")
  join_data = left_join(peak_num_aligned_long, myannot_grouped, by = c("run" = "run", "peakid" = "peak_ID")) #join annotation and peak-aligned
  join_data$aligned = as.factor(!is.na(join_data$peakid)) #convert factors
  return(join_data)
}


