#'Start GUI
#'@description start GUI as shiny app.
#'@usage launchApp()
#'@return shiny app
#'@export
launchApp<-function(){
  shiny::runApp("App-GC_1.0",display.mode = "normal",launch.browser = TRUE)
}
