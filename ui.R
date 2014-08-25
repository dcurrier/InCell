
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

require(shiny)             # Shiny Framework
require(shinyIncubator)    # Progress indicator
require(shinythings)       # Password Input/Better Action buttons
require(ShinyHighCharts)   # Javascript charting
options(shiny.maxRequestSize=45*1024^2)

shinyUI(fluidPage(
  progressInit(),     # Initialize the Prgoress Indicator

  # Application title
  titlePanel("InCell Data Analysis"),

  # Main Application
  fluidRow(
    column(4,
           tabsetPanel(id="tabs",
             tabPanel("Data",
                      div(
                        conditionalPanel(
                          condition="output.fileUploaded",
                          fileInput('InCell', label=h4("Import InCell Data File")),
                          fileInput('annotation', label=h4("REMP Plate Lookup File")),
                          hr(),
                          h4("Choose from Dropbox"),
                          selectizeInput('InCellDB', label="",
                                         choices=c("InCell File", "No Files Found"), selected="InCell File"),
                          selectizeInput('annotationDB', label="",
                                         choices=c("REMP File", "No Files Found"), selected="REMP File"),
                          textInput('ctlCols', label=h4("Control Columns"), value="example: 21-24 or 21,22,23,24"),
                          actionButton("updateDropBox", label="Check for new files"),
                          hr()
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
             tabPanel("QC",
                      div(
                        selectizeInput('featureCol', label=h4("Select a Feature"),
                                       choices=c("Cell Count"),
                                       selected="Cell Count", width="100%"),
                        selectizeInput('fieldWell', label=h4("Choose a Well"),
                                       choices=c("A01"), selected="A01"),
                        hr(),
                        verbatimTextOutput('wellData'),
                        class="well tab-well")),
             tabPanel("Feature",
                      div(
                        selectizeInput('featureColDist', label=h4("Select a Feature"),
                                       choices=c("Cell Count"),
                                       selected="Cell Count", width="100%"),
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
                        checkboxInput('logFeatureValues', label="Log Transform Feature Values", value=F),
                        hr(),
                        h5("Negative Controls"),
                        fluidRow(
                                 column(6,verbatimTextOutput('ksStatSummary')),
                                 column(6,verbatimTextOutput('cvmtsStatSummary'))
                        ),
                        fluidRow(
                          column(6, h5("Selected Well"),
                                    verbatimTextOutput('cmpdStatSummary'))
                        ),
                        class="well tab-well")),
             position='left'
           )),
    column(8,
           conditionalPanel(
             condition="input.tabs == 'Data'",
             div( class="padding-bottom: 30px")
           ),
           conditionalPanel(
             condition="input.tabs == 'QC'",
             div( highchartsOutput("highHeat", height="470px", include=c("base", "more", "heatmap", "map", "no-data")),
                                   class="padding-bottom: 30px"),
             br(),
             fluidRow(
               column(5, offset=1, highchartsOutput('miniHist', height="250px", include=c("base", "more", "export", "no-data") )),
               column(5, highchartsOutput('miniFieldData', height="250px", include=c("base", "more", "export", "no-data") ))
               )
             ),
           conditionalPanel(
             condition="input.tabs == 'Feature'",
             div( highchartsOutput("distribution", height="470px", include=c("base", "more", "no-data")),
                  class="padding-bottom: 30px"),
             br(),
             fluidRow(
               column(5, offset=1, highchartsOutput('miniZStat', height="250px", include=c("base", "more", "export", "no-data") )),
               column(5, highchartsOutput('miniAUC', height="250px", include=c("base", "more", "export", "no-data") ))
             )
           )
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
