library(xcms) #AS: added
source('App-GC/App-GC_GCxGC.R') #AS: added
library(shiny)
library(shinythemes)
library(shinyFiles)
library(shinyjs)
library(plotly)
library(fs)
library(filesstrings)

#options(browser="/usr/bin/google-chrome")

ui <- fluidPage(
    useShinyjs(),
    titlePanel(title="", windowTitle="App-GC pipeline"),
    ##-- head
    tags$head(
      tags$link(rel = "shortcut icon", href = "img/favicon.png"),
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    ##-- top menu ----
    navbarPage(title = div(img(src="img/logo.png", height = "50px")), theme = shinytheme("flatly"),
               id = "navbar", selected = "homepage", fluid = T,
               tabPanel(title = "Home", value = "homepage"),
               tabPanel(title = "GC", value = "gcpage"),
               tabPanel(title = "GCxGC", value = "gc2dpage"),
               tabPanel(title = "Help", value = "helppage")


    ),
    ##-- main content ----
    fluidRow(
      column(12,
             # This outputs the dynamic UI component
             uiOutput("ui")
      )
    ),
    ##-- footer ----
    div(class = "text-center mb-xl-5", hr(), 
        HTML("<span>©2021. All rights reserved.</span>")
    )
)

