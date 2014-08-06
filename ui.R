
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)             # Shiny Framework
library(shinyIncubator)    # Progress indicator
library(shinythings)       # Password Input/Better Action buttons
options(shiny.maxRequestSize=45*1024^2)

shinyUI(fluidPage(
  progressInit(),     # Initialize the Prgoress Indicator

  # Application title
  titlePanel("InCell Data Analysis"),

  # Main Application
  fluidRow(
    column(4,
           tabsetPanel(
             tabPanel("Data",
                      div(
                        conditionalPanel(
                          condition="output.fileUploaded",
                          fileInput('InCell', label=h4("Import InCell Data File"))
                          ),
                        conditionalPanel(
                          condition="!output.fileUploaded",
                          h4('Current Dataset:'),
                          verbatimTextOutput('dataset'),
                          hr(),
                          h4("Download Individual Tables"),
                          fluidRow(
                            column(6, p("Well Level Data")),
                            column(6,  downloadButton('downWell', label="Download"))
                            ),
                          br(),
                          fluidRow(
                            column(6, p("Field Level Data")),
                            column(6,  downloadButton('downField', label="Download"))
                            ),
                          br(),
                          fluidRow(
                            column(6, p("Cell Level Data")),
                            column(6,  downloadButton('downCell', label="Download"))
                            )
                        ),
                        actionButton('console', label="Show Console"),
                        class="well tab-well"
                      )),
             tabPanel("Well",
                      div(
                        selectizeInput('wellColumn', label=h4("Heatmap Data"),
                                       choices=c("Cell Count"),
                                       selected="Cell Count", width="100%"),
                        hr(),
                        verbatimTextOutput('wellData'),
                        hr(),
                        plotOutput('miniHist', height="250px"),
                        class="well tab-well")),
             tabPanel("Field",
                      div(
                        selectizeInput('fieldColumn', label=h4("Data to Plot"),
                                       choices=c("Cell Count"),
                                       selected="Cell Count", width="100%"),
                        selectizeInput('fieldWell', label=h4("Well of Interest"),
                                       choices=c("A1"), selected="A1", width="100%"),
                        hr(),
                        h5("Cell Level Detail"),
                        plotOutput('miniFieldData', height="350px"),
                        verbatimTextOutput('FieldSummary'),
                        class="well tab-well")),
             tabPanel("Cell",
                      div(

                        class="well tab-well")),
             position='left'
           )),
    column(8,
           carouselPanel(
              plotOutput('heatmap', height="500px"),
              plotOutput('fieldwisePlot', height="500px")
           ))
  ),
  tags$head(
    tags$style("
                 .nav-tabs{
                  border-right: none !important;
                  margin-right: 0px !important;
                 }
                 .nav-tabs .active a, .nav-tabs .active a:focus {
                  background-color: #F5F5F5;
                 }
                 .tab-well{
                  border-radius: 0px 4px 4px 4px;
                  height:630px;
                 }
                 .color{
                  width: 60px;
                 }
                 "
    ),
    tags$script(type="text/javascript", src="iframeResizer.contentWindow.min.js")
  )
))
