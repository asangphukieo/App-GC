source("util.R", local=TRUE)
source("generateDat.R", local=TRUE)
source("generatePlot.R", local=TRUE)
source("generateHeatmap.R", local=TRUE)
source("formatOutputData.R", local=TRUE)


server <- function(input, output, session) {
  ## -- Render UI
  output$ui <- renderUI({
    if (is.null(input$navbar))
      return()
    output$summary <- renderText({paste("Ready to run a job ...")}) #reset
    selected_area$coord <- NULL #reset
    switch(input$navbar,
           "homepage" = source("pages/homepage.R", local=TRUE)[1],
           "gcpage" = source("pages/gcpage.R", local=TRUE)[1],
           "gc2dpage" = source("pages/gc2dpage.R", local=TRUE)[1],
           "helppage" = source("pages/helppage.R", local=TRUE)[1]
    )
  })
  ## -- Reactive variables
  Input_navbar_var <- reactive({#mode
    navbar <- input$navbar
    if (navbar == "gcpage"){
      "1D"
    }else if(navbar == "gc2dpage"){
      "2D"
    }else{
      return(NULL) 
    }
  })
  gc_config_param <- reactive({#gc configuration parameters
    paste0("mode: ", Input_navbar_var(),
           "\nInput_path : ", input$Input_path_txt,
           "\nInput_file : \"", input$Input_file,"\"",
           "\nOutput_folder : ", input$output_folder_txt,
           "\nGC:",
           "\n Peak_detection : ",
           "\n  do_peak_detection : ", input$do_peak_detection,           
           "\n  method : ", input$method_peak_gc,
           "\n  worker : ", input$worker, 
           "\n  CAMERA_perfwhm : ", input$CAMERA_perfwhm, 
           "\n  to2D : ", input$to2D, 
           "\n  modulation_time : ", input$to2Dmodulation_time,                       
           "\n  centWaveParam : ",       
           "\n   snthresh : ", input$snthresh,          
           "\n   ppm : ", input$ppm,
           "\n   peakwidth_min : ", input$peakwidth_min,
           "\n   peakwidth_max : ", input$peakwidth_max,
           "\n  matchedFilterParam : ",    
           "\n   snthresh : ", input$snthresh,
           "\n   fwhm : ", input$fwhm,
           "\n   step : ", input$step,
           "\n   steps : ", input$steps,
           "\n   mzdiff : ", input$mzdiff,
           "\n   max : ", input$max,
           "\n  erahParam : ",
           "\n   noise.threshold : ", input$noise.threshold,
           "\n   min.peak.height : ", input$min.peak.height,
           "\n   min.peak.width : ", input$min.peak.width,
           "\n Peak_alignment :",
           "\n  do_peak_alignment : ", input$do_peak_alignment,
           "\n  method : ", input$method_align_gc,           
           "\n  reference : ", input$reference_txt,
           "\n  GCalignRParam : ", input$GCalignRParam,       
           "\n   max_linear_shift : ", input$max_linear_shift,
           "\n   max_diff_peak2mean : ", input$max_diff_peak2mean,
           "\n   min_diff_peak2peak : ", input$min_diff_peak2peak,               
           "\n   delete_single_peak : ", input$delete_single_peak,
           "\n  ERahParam : ", input$ERahParam,
           "\n   min.spectra.cor : ", input$min.spectra.cor,
           "\n   max.time.dist : ", input$max.time.dist,                      
           "\nPeak_annotation : ",
           "\n do_peak_annotation : ", input$do_peak_annotation,
           "\n path_MSPepSearch : ", input$path_MSPepSearch_txt,
           "\n MinInt : ", input$MinInt,
           "\n MinMF : ", input$MinMF,
           "\n HITS : ", input$HITS,
           "\n path_MAIN : ", input$path_MAIN_txt,
           "\n path_REPL : ", input$path_REPL_txt,
           "\nPeak_summary : ",
           "\n do_peak_summary : ", input$do_peak_summary,
           "\n reference : ", input$reference_sum,
           "\n select_name : ", input$select_name)
  })
  gc2d_config_param <- reactive({#gc2d configuration parameters
    paste0("mode: ", Input_navbar_var(),
           "\nInput_path : ", input$Input_path_txt,
           "\nInput_file : \"", input$Input_file,"\"",
           "\nOutput_folder : ", input$output_folder_txt,
           "\nGC2D:",
           "\n Peak_detection : ",
           "\n  do_peak_detection : ", input$do_peak_detection,
           "\n  method : ", input$method_peak_gc2d,
           "\n  method_parallel : ", input$method_parallel,
           "\n  zcv : ", input$zcv,
           "\n  cv : ", input$cv,     
           "\n  peak_merge : ", input$peak_merge,      
           "\n  reference : ", input$reference_gc2d,	#AS: added
           "\n  worker : ", input$worker,
           "\n  sampling_rate : ", input$sampling_rate,
           "\n  modulation_time : ", input$modulation_time,
           "\n  profMat_method : ", input$profMat_method,       
           "\n  profMat_step : ", input$profMat_step,
           "\n  rm_prev_tic : ", input$rm_prev_tic,                         
           "\n  full_scan : ", input$full_scan,
           "\n  time_scan :",
           "\n   RT1Dstart : ", input$RT1Dstart,
           "\n   RT1Dstop : ", input$RT1Dstop,
           "\n   RT2Dstart : ", input$RT2Dstart,
           "\n   RT2Dstop : ", input$RT2Dstop,
           "\n  useXY : TRUE",
           "\n  x1 : ", selected_area$coord[['x1']],
           "\n  x2 : ", selected_area$coord[['x2']],
           "\n  y1 : ", selected_area$coord[['y1']],
           "\n  y2 : ", selected_area$coord[['y2']],
           "\n Peak_alignment :",
           "\n  do_peak_alignment : ", input$do_peak_alignment,
           "\n  reference : ", input$reference_align_gc2d,
           "\n  method : ", input$method_align_gc2d,   
           "\n  mSPAParam : ",
           "\n   method : ", input$mSPA_method,           
           "\n   w : ", input$mSPA_w, 
           "\n   k : ", input$mSPA_k,
           "\n   rho : ", input$mSPA_rho,
           "\n   distance : ", input$mSPA_distance,
           "\n   arraySimilarity : ", input$mSPA_arraySimilarity,
           "\n   createReferenceAlignment : ",input$mSPA_createReferenceAlignment,
           "\n   referenceAlignmentVariant : ",input$mSPA_referenceAlignmentVariant,
           "\n   sd1Thres : ",input$mSPA_sd1Thres,
           "\n   sd2Thres : ",  input$mSPA_sd2Thres,
           "\n  ERahParam : ", 
           "\n   min.spectra.cor : ",  input$ERah_min.spectra.cor,  
           "\n   max.time.dist : ",  input$ERah_max.time.dist,                                                     
           "\nPeak_annotation : ",
           "\n do_peak_annotation : ", input$do_peak_annotation,
           "\n path_MSPepSearch : ", input$path_MSPepSearch_txt,
           "\n MinInt : ", input$MinInt,
           "\n MinMF : ", input$MinMF,
           "\n HITS : ", input$HITS,
           "\n path_MAIN : ", input$path_MAIN_txt,
           "\n path_REPL : ", input$path_REPL_txt,
           "\nPeak_summary : ",
           "\n do_peak_summary : ", input$do_peak_summary,
           "\n reference : ", input$reference_sum,
           "\n select_name : ", input$select_name)
  })
  selected_area <- reactiveValues(coord = list()) #selected coordinate
  out_values <- reactiveValues()
  ## -- Call functions
  runPipeline <- function(){
    if(input$navbar == "gcpage"){
      out_values$outdat = NULL #reset
      out_values$outdat = generateDat(15) #call gc function
      	#setwd('../')
	print(getwd())
	current_path<- getwd()
	source('run_App-GC.R', local=TRUE)#AS
	setwd(current_path)	
	#AS: Added//  
    }
    if(input$navbar == "gc2dpage"){
      out_values$outdat = NULL #reset
      out_values$outdat = generateDat2(selected_area$coord[['x1']]) #call gc2d function
      	#setwd('../')
	print(getwd())
	current_path<- getwd()
	source('run_App-GC.R', local=TRUE)#AS
	setwd(current_path)	
	#AS: Added//        
      
    }
    if(input$navbar == "helppage"){
      out_values$outdat = NULL #reset
    }
  }
  returnResult <- reactive({ #return final output
    runPipeline()
    sum_result = out_values$outdat
    if(input$do_peak_summary==TRUE){
    	sum_result = read.table(paste(input$output_folder_txt,"/SUMMARY/peak_summary.txt",sep=''),header=T,sep='\t')
    }
    sum_result
    
  })
  ## -- Event listeners
  observeEvent(input$Input_path,{#selected Input_path
    inDirPath = choose_directory()
    updateTextInput(session, "Input_path_txt",
                    label = "", value = path(inDirPath))
  })
  observeEvent(input$Output_folder,{#selected Output_folder
    outDirPath = choose_directory()
    updateTextInput(session, "output_folder_txt",
                    label = "", value = path(outDirPath))
  })
  observeEvent(input$reference,{#selected reference
    refDirPath = choose_directory()
    updateTextInput(session, "reference_txt",
                    label = "", value = path(refDirPath))
  })
  observeEvent(input$path_MSPepSearch,{#selected reference
    MSPepDirPath = choose_directory()
    updateTextInput(session, "path_MSPepSearch_txt",
                    label = "", value = path(MSPepDirPath))
  })
  observeEvent(input$path_MAIN,{#selected reference
    MAINDirPath = choose_directory()
    updateTextInput(session, "path_MAIN_txt",
                    label = "", value = path(MAINDirPath))
  })
  observeEvent(input$path_REPL,{#selected reference
    REPLDirPath = choose_directory()
    updateTextInput(session, "path_REPL_txt",
                    label = "", value = path(REPLDirPath))
  })
  ## -- GC
  observeEvent(input$apply_gc,{#generate gc configuration file
    output$summary <- renderText({ #list configuration parameters
      paste0("Configuration parameters:\n"
             ,gc_config_param(),
             "\n\nConfiguration file was generated to : ", getwd())
    })      
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Generating configuration file ", value = 0)
      for (i in 1:10) {
        progress$inc(0.1, detail = "...")
      }
        
      Sys.sleep(1)
      write(gc_config_param(), file = 'running_App-GC_config.yml') #write to file
      enable("submit_gc")
      hide("downloadgc_btn")
    })  
  })
  observeEvent(input$submit_gc,{#run gc pipeline
    output$summary <- renderPrint({
      isolate({
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Running ", value = 0)
        for (i in 1:10) {
          progress$inc(0.1, detail = "...")
        }
 
        Sys.sleep(1)
        str(returnResult())
        cat("\nJob DONE!")
      })
    })
    show("downloadgc_btn")
    show("view_alignment_btn") #call alignment heatmap
    disable("submit_gc")
  })
  output$download_gc <- downloadHandler(#download output
    filename = function() { 
      paste("peak_summary.txt",sep='')
    },
    content = function(file) {
      write.csv(returnResult(), file) #write to file
    }
  )

  ############################### GCxGC ###################################
  observeEvent(input$apply_gc2d,{#generate gc2d configuration file
    output$summary <- renderText({ #list configuration parameters
      paste0("Configuration parameters:\n"
             ,gc2d_config_param(),
             "\n\nMESSAGE: Configuration file will be generated to : ", getwd())
    })
    isolate({
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Running ", value = 0)
      for (i in 1:10) {
        progress$inc(0.1, detail = "...")
      }
      Sys.sleep(1)
      write(gc2d_config_param(), file = 'running_App-GC_config.yml') #write to file
      enable("proceed_gc2d")
      disable("submit_gc2d")
      hide("downloadgc2d_btn")
      hide("graph_box")
      selected_area$coord <- NULL #reset
    })
  })
  observeEvent(input$proceed_gc2d,{#render graph
    show("graph_box")
    output$graph <- renderPlotly({
      isolate({
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Running ", value = 0)
        for (i in 1:10) {
          progress$inc(0.1, detail = "...")
        }
        ref_path=generatePlot(sampling_rate=input$sampling_rate, modulation_time=input$modulation_time, reference=FALSE, Input_file=input$Input_file, Input_path=input$Input_path_txt,profMat_method=input$profMat_method, profMat_step=input$profMat_step, rm_prev_tic=input$rm_prev_tic,ref.dir=input$output_folder_txt) #call and pass params to generatePlot function
      })
    })
    disable("proceed_gc2d")
    output$summary <- renderText({#print coordinate
      paste0("Use 'Box Select' to select peak area")
    })
  })
  observeEvent(event_data("plotly_brushed", source="graph"), {#selected area
    sel <- event_data("plotly_brushed", source = "graph")
    selected_area$coord[['x1']] <- as.integer(sel$x[1])
    selected_area$coord[['x2']] <- as.integer(sel$x[2])
    selected_area$coord[['y1']] <- as.integer(sel$y[1])
    selected_area$coord[['y2']] <- as.integer(sel$y[2])
    enable("submit_gc2d")
    output$summary <- renderText({#print coordinate
      paste0("Selected area: ",
             "\nX1 : ", selected_area$coord[['x1']],
             "\nX2 : ", selected_area$coord[['x2']],
             "\nY1 : ", selected_area$coord[['y1']],
             "\nY2 : ", selected_area$coord[['y2']])
    })
  })
  observeEvent(input$submit_gc2d,{#run gc2d pipeline
    hide("graph_box")
    write(gc2d_config_param(), file = 'running_App-GC_config.yml') #write to file #AS: Edited (../)
    output$summary <- renderPrint({
      isolate({
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Running ", value = 0)
        for (i in 1:10) {
          progress$inc(0.1, detail = "...")
        }

        Sys.sleep(1)
        cat(c("Configuration parameters:\n",gc2d_config_param(),"\nMESSAGE: Configuration file was generated to : ", getwd(),
              "\n\nOUTPUT:"))
        str(returnResult())

        cat("\nJob DONE!")
      })
    })
    show("downloadgc2d_btn")
    show("view_alignment_btn") #call alignment heatmap
    disable("submit_gc2d")
  })  
  output$download_gc2d <- downloadHandler(#download output
    filename = function() { 
      paste("peak_summary.txt",sep='')
    },
    content = function(file) {
      write.csv(returnResult(), file) #write to file
    }
  )
  
  #render alignment heatmap
  aligned_data <- reactive({#get aligned data
    ali_folder=paste(input$output_folder_txt,"/PEAK_ALIGNMENT/save_ali_peak.Rdata",sep='')
    annot_folder=paste(input$output_folder_txt,"/PEAK_ANNOTATION/save_annotation_peak.Rdata",sep='')
    formatOutputData(ali_folder,annot_folder)
  })
  observeEvent(input$view_alignment,{#render alignment heatmap
    show("heatmap_box")
    #show(aligned_data()) #AP
    output$heatmap <- renderPlotly({
      isolate({
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Running ", value = 0)
        for (i in 1:10) {
          progress$inc(0.1, detail = "...")
        }
        generateHeatmap(aligned_data())
      })
    })
    hide("summary_box")
    #disable("submit_gc2d")
  })
  hmselection <- reactive({
    sel_ind <- event_data("plotly_click",source = "heatmap")
    sel_ind$key
  })
  observeEvent(event_data("plotly_click", source="heatmap"), {#render table of selected alignment
    sel_dt <- dplyr::filter(aligned_data(), mean_RT %in% hmselection())
    output$heatmaptb = renderDataTable(sel_dt, options = list(dom = 't'))
  })

  observeEvent(input$reset, {#reset
    reset("gcapp")
    output$summary <- renderText({paste("Ready to run a job ...")}) 
    disable("submit_gc")
    hide("downloadgc_btn")
    reset("gc2dapp")
    disable("proceed_gc2d")
    disable("submit_gc2d")
    hide("downloadgc2d_btn")
    output$graph <- renderPlotly({plotly_empty()})
    show("summary_box")
    hide("graph_box")
    selected_area$coord <- NULL
    output$heatmap <- renderPlotly({plotly_empty()})
    output$heatmaptb <- renderDataTable(data.frame(), options = list(dom = 't'))
    hide("heatmap_box")    
  })
}
