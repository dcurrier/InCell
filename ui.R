
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
    ##################### Inputs Pane #####################
    div(tabsetPanel(id="tabs",
                    ##################### Data Tab #####################
                    tabPanel("Data",
                             # Help Button
                             div(class="pull-right",
                                 modalButton(label="", styleclass="link", size="large", icon="question-circle", modal.id="Datahelp",
                                             css.class="help fa-lg")),
                             div(
                               conditionalPanel(
                                 condition="output.fileUploaded || output.fileUploaded == null",
                                 fileInput('InCell', label=h4("Import InCell Data File")),
                                 selectizeInput('InCellDB', label="",
                                                choices=c("Choose from Dropbox", "No Files Found"), selected="Choose from Dropbox "),
                                 hr(),
                                 fileInput('annotation', label=h4("REMP Plate Lookup File")),
                                 selectizeInput('annotationDB', label="",
                                                choices=c("Choose from Dropbox", "No Files Found"), selected="Choose from Dropbox "),
                                 checkboxInput('noCtlPlate', label="All control wells contain DMSO", value=FALSE),
                                 hr(),
                                 actionButton("updateDropBox", label="Check for new files")
                               ),
                               conditionalPanel(
                                 condition="!output.fileUploaded && output.fileUploaded != null",
                                 h4('Current Dataset:'),
                                 verbatimTextOutput('dataset'),
                                 hr(),
                                 textInput('ctlCols', label=h4("Control Columns"), value="example: 21-24 or 21,22,23,24"),
                                 selectizeInput('setPos', label=h4("Set Positive Control"), choices=c("Select a Compound")),
                                 selectizeInput('posConc', label="", choices="Choose a Concentration"),
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
                               #actionButton('console', label="Show Console"),
                               class="well tab-well"
                             )),
                    ##################### QC Tab #####################
                    tabPanel("QC",
                             # Help Button
                             div(class="pull-right",
                                 modalButton(label="", styleclass="link", size="large", icon="question-circle", modal.id="QChelp",
                                             css.class="help fa-lg")),
                             div(
                               selectizeInput('featureCol', label=h4("Select a Feature"),
                                              choices=c("Cell Count"),
                                              selected="Cell Count", width="100%"),
                               hr(),
                               verbatimTextOutput('wellData'),
                               class="well tab-well")),
                    ##################### Any|Any Tab #####################
                    tabPanel("Any|Any",
                             # Help Button
                             div(class="pull-right",
                                 modalButton(label="", styleclass="link", size="large", icon="question-circle", modal.id="AAhelp",
                                             css.class="help fa-lg")),
                             div(
                               h4("Y Axis Feature"),
                               selectizeInput('yFeat', label="",
                                              choices=c("Cell Count"),
                                              selected="Cell Count", width="100%"),
                               fluidRow(
                                 column(6, radioButtons('yTrans', label="Y Axis Transform",
                                                        choices=c("None", "Log10", "Log2"))   ),
                                 column(6,radioButtons('yThresh', label="Y Axis Threshold",
                                                       choices=c("None", "% Above", "% Below"))   )
                               ),
                               uiOutput('ySlider'),
                               hr(),
                               h4("X Axis Feature"),
                               selectizeInput('xFeat', label="",
                                              choices=c("Cell Count"),
                                              selected="Cell Count", width="100%"),
                               fluidRow(
                                 column(6, radioButtons('xTrans', label="X Axis Transform",
                                                        choices=c("None", "Log10", "Log2"))   ),
                                 column(6,radioButtons('xThresh', label="X Axis Threshold",
                                                       choices=c("None", "% Above", "% Below"))   )
                               ),
                               uiOutput('xSlider'),
                               hr(),
                               h4("Add Lines"),
                               checkboxInput('hLineShow', label="Horizontal", value=F),
                               fluidRow(
                                 column(7, textInput('hLineName', label="", value="Label")),
                                 column(4, offest=1, numericInput('hLine', label="", value=0))
                               ),
                               br(),
                               checkboxInput('vLineShow', label="Vertical", value=F),
                               fluidRow(
                                 column(7, textInput('vLineName', label="", value="Label")),
                                 column(4, offest=1, numericInput('vLine', label="", value=0))
                               ),
                               class="well tab-well")),
                    ##################### Any|Conc Tab #####################
                    tabPanel("Any|Conc",
                             # Help Button
                             div(class="pull-right",
                                 modalButton(label="", styleclass="link", size="large", icon="question-circle", modal.id="AChelp",
                                             css.class="help fa-lg")),
                             div(
                               h4("Y Axis Feature"),
                               selectizeInput('feat', label="",
                                              choices=c("Cell Count"),
                                              selected="Cell Count", width="100%"),
                               fluidRow(
                                 column(6, radioButtons('featTrans', label="Y Axis Transform",
                                                        choices=c("None", "Log10", "Log2"))   ),
                                 column(6,radioButtons('featThresh', label="Y Axis Threshold",
                                                       choices=c("None", "% Above", "% Below"))   )
                               ),
                               uiOutput('featSlider'),

                               class="well tab-well")),
                    ##################### Feature Tab #####################
                    tabPanel("Feature",
                             # Help Button
                             div(class="pull-right",
                                 modalButton(label="", styleclass="link", size="large", icon="question-circle", modal.id="Feathelp",
                                             css.class="help fa-lg")),
                             div(
                               selectizeInput('featureColDist', label=h4("Select a Feature"),
                                              choices=c("Cell Count"),
                                              selected="Cell Count", width="100%"),
                               selectizeInput('cmpdSelect', label=h4("Select a Compound"),
                                              choices=c(""), selected="", width="100%"),
                               fluidRow(
                                 column(6, div( actionButton('prv', label='', styleclass='primary',
                                                             icon='long-arrow-left', icon.library='font awesome'),
                                                class="pull-left")),
                                 column(6, div( actionButton('nxt', label='', styleclass='primary',
                                                             icon='long-arrow-right', icon.library='font awesome'),
                                                class="pull-right"))
                               ),
                               br(),
                               radioButtons('logFeatureValues', label="Log Transform Feature Values",
                                            choices=c("None", "Log10", "Log2"), selected="None"),
                               hr(),
                               h5("Negative Controls"),
                               tableOutput('statSummary'),
                               class="well tab-well")),
                    ##################### Activity Tab #####################
                    tabPanel("Activity",
                             # Help Button
                             div(class="pull-right",
                                 modalButton(label="", styleclass="link", size="large", icon="question-circle", modal.id="Acthelp",
                                             css.class="help fa-lg")),
                             div(
                               selectizeInput('cmpdAct', label=h4("Select a Compound"),
                                              choices=c(""), selected="", width="100%"),
                               fluidRow(
                                 column(6, div( actionButton('prvAct', label='', styleclass='primary',
                                                             icon='long-arrow-left', icon.library='font awesome'),
                                                class="pull-left")),
                                 column(6, div( actionButton('nxtAct', label='', styleclass='primary',
                                                             icon='long-arrow-right', icon.library='font awesome'),
                                                class="pull-right"))
                               ),
                               class="well tab-well")),
                    position='left'
    ),
    class="span4",
    style="max-width: 400px;"),
    ##################### Content Pane #####################
    column(8,
           ##################### Data Pane #####################
           # Data Display Panel
           conditionalPanel(
             condition="input.tabs == 'Data'",
             div( conditionalPanel(
               condition="output.fileUploaded || output.fileUploaded == null",
               img(src="Instructions/Slide1.png")),
               conditionalPanel(
                 condition="!output.fileUploaded && output.fileUploaded != null",
                 img(src="Instructions/Slide2.png")),
               class="padding-bottom: 30px")
           ),
           ##################### QC Pane #####################
           # QC Display Panel - Heatmap, Mini Histogram, Mini Field Scatter Plot
           conditionalPanel(
             condition="input.tabs == 'QC'",
             div(
               # No files have been uploaded
               conditionalPanel(
                 condition="output.fileUploaded || output.fileUploaded == null",
                 img(src="Instructions/NoPlot.png")),
               # Files are ready
               conditionalPanel(
                 condition="!output.fileUploaded && output.fileUploaded != null",
                 div( highchartsOutput("highHeat", height="470px", include=c("base", "more", "heatmap", "map", "no-data")),
                      class="padding-bottom: 30px"),
                 br(),
                 fluidRow(
                   column(5, offset=1, highchartsOutput('miniHist', height="250px", include=c("base", "more", "export", "no-data") )),
                   column(5, highchartsOutput('miniFieldData', height="250px", include=c("base", "more", "export", "no-data") ))
                 )
               )
             )
           ),
           ##################### Any|Any Pane #####################
           # Any over Any Display Panel - Plot any feature over any other
           conditionalPanel(
             condition="input.tabs == 'Any|Any'",
             # No files have been uploaded
             div(
               conditionalPanel(
                 condition="output.fileUploaded || output.fileUploaded == null",
                 img(src="Instructions/NoPlot.png")),
               # Files are ready
               conditionalPanel(
                 condition="!output.fileUploaded && output.fileUploaded != null",
                 div( highchartsOutput('AnyAny', height="470px", include=c("base", "more", "export", "no-data", "dim-on-hover") ),
                      class="padding-bottom: 30px"),
                 br(),
                 fluidRow(
                   column(6, highchartsOutput('yThreshPlot', height="250px", include=c("base", "more", "export", "no-data") )),
                   column(6, highchartsOutput('xThreshPlot', height="250px", include=c("base", "more", "export", "no-data") ))
                 )
               )
             )
           ),
           ##################### Any|Conc Pane #####################
           # Any over Any Display Panel - Plot any feature over any other
           conditionalPanel(
             condition="input.tabs == 'Any|Conc'",
             # No files have been uploaded
             div(
               conditionalPanel(
                 condition="output.fileUploaded || output.fileUploaded == null",
                 img(src="Instructions/NoPlot.png")),
               # Files are ready
               conditionalPanel(
                 condition="!output.fileUploaded && output.fileUploaded != null",
                 div( highchartsOutput('AnyConc', height="470px", include=c("base", "more", "export", "no-data", "dim-on-hover") ),
                      class="padding-bottom: 30px"),
                 br(),
                 fluidRow(
                   column(6, highchartsOutput('featThreshPlot', height="250px", include=c("base", "more", "export", "no-data") ))
                 )
               )
             )
           ),
           ##################### Feature Pane #####################
           # Feature Display Panel - Kernel Density Plot, Z-Stat Plot, AUC Plot
           conditionalPanel(
             condition="input.tabs == 'Feature'",
             # No files have been uploaded
             div(
               conditionalPanel(
                 condition="output.fileUploaded || output.fileUploaded == null",
                 img(src="Instructions/NoPlot.png")),
               # Files are ready
               conditionalPanel(
                 condition="!output.fileUploaded && output.fileUploaded != null",
                 div( highchartsOutput("distribution", height="470px", include=c("base", "more", "no-data")),
                      class="padding-bottom: 30px"),
                 br(),
                 fluidRow(
                   column(6, offset=1, highchartsOutput('miniZStat', height="250px", include=c("base", "more", "export", "no-data") )),
                   column(4, highchartsOutput('miniAUC', height="250px", include=c("base", "more", "export", "no-data") ))
                 )
               )
             )
           ),
           ##################### Action Pane #####################
           # Mechanism of Action Panel -
           conditionalPanel(
             condition="input.tabs == 'Activity'",
             # No files have been uploaded
             div(
               conditionalPanel(
                 condition="output.fileUploaded || output.fileUploaded == null",
                 img(src="Instructions/NoPlot.png")),
               # Files are ready
               conditionalPanel(
                 condition="!output.fileUploaded && output.fileUploaded != null",
                 div( highchartsOutput("allAUC", height="270px", include=c("base", "more", "no-data")),
                      class="padding-bottom: 30px"),
                 br(),
                 fluidRow(
                 column(6, div( highchartsOutput("clusterAUC", height="370px", include=c("base", "more", "heatmap", "map", "no-data")),
                      class="padding-bottom: 30px")),
                 column(6, div())
                 )
               )
             )
           )
    )
  ),
  ##################### Help Modals #####################
  modalPanel("Datahelp", header=h2("Data Import"), body=Data),
  modalPanel("QChelp", header=h2("Quality Control"), body=QC),
  modalPanel("AAhelp", header=h2("Any Feature versus Any Feature"), body=AA),
  modalPanel("AChelp", header=h2("Any Feature versus Concentration"), body=AC),
  modalPanel("Feathelp", header=h2("Feature Value Distribution"), body=Feat),
  modalPanel("Acthelp", header=h2("Mechanism of Action Prediction"), body=Act),

  ##################### Header Tags #####################
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
                .help, .help:focus{
                 color: #F89406;
                }
                .help:hover{
                 color: #F88300;
                }
                 "
    )
    #tags$script(type="text/javascript", src="iframeResizer.contentWindow.min.js")
  )
))
