
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
                          fileInput('InCell', label=h4("Import InCell Data File")),
                          fileInput('annotation', label=h4("REMP Plate Lookup File"))
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
                        br(),
                        actionButton('console', label="Show Console"),
                        class="well tab-well"
                      )),
             tabPanel("Well",
                      div(
                        verbatimTextOutput('wellData'),
                        hr(),
                        plotOutput('miniHist', height="250px"),
                        class="well tab-well")),
             tabPanel("Field",
                      div(
                        h5("Cell Level Detail"),
                        plotOutput('miniFieldData', height="350px"),
                        verbatimTextOutput('FieldSummary'),
                        class="well tab-well")),
             tabPanel("Cell",
                      div(
                        class="well tab-well")),
             tabPanel("Feature",
                      div(
                        textInput('ctlCols', label=h4("Control Columns"), value="example: 21-24 or 21,22,23,24"),
                        fluidRow(
                          column(8, uiOutput('concSlide')),
                          column(4, div(verbatimTextOutput('selConc'), style="margin-top: 21px"))
                        ),
                        hr(),
                        selectizeInput('cmpdSelect', label=h4("Select a Compound"),
                                       choices=c(""), selected=""),
                        fluidRow(
                          column(6, div( actionButton('prv', label='', styleclass='primary',
                                                      icon='chevron-left', icon.library='font awesome'),
                                         class="pull-left")),
                          column(6, div( actionButton('nxt', label='', styleclass='primary',
                                                      icon='chevron-right', icon.library='font awesome'),
                                         class="pull-right"))
                        ),
                        hr(),
                        h5("Negative Controls"),
                        fluidRow(
                                 column(6,verbatimTextOutput('ksStatSummary')),
                                 column(6,verbatimTextOutput('cvmtsStatSummary'))
                        ),
                        fluidRow(
                          column(6, h5("Selected Well"),
                                    verbatimTextOutput('cmpdStatSummary')),
                          column(6, h5("Selected Compound"),
                                    plotOutput('cmpdZplot', height="161px"))
                        ),
                        class="well tab-well")),
             position='left'
           )),
    column(8,
           fluidRow(
             carouselPanel(
                plotOutput('heatmap', height="500px"),
                plotOutput('fieldwisePlot', height="500px"),
                plotOutput('featureDistPlot', height="500px"),
                plotOutput('featureDistGrid', height="500px")
             ))
           ),
           fluidRow(
             column(3,
                    wellPanel(div(
                      selectizeInput('featureCol', label=h4("Select a Feature"),
                                      choices=c("Cell Count"),
                                      selected="Cell Count", width="100%"),
                      selectizeInput('fieldWell', label=h4("Choose a Well"),
                                     choices=c("A01"), selected="A01"),
                        style="height: 169px;")))
             )
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
                  height:700px;
                 }
                 .color{
                  width: 60px;
                 }
                 "
    ),
    tags$script(type="text/javascript", src="iframeResizer.contentWindow.min.js")
  )
))
