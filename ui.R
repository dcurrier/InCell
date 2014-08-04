
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
    column(5,
           tabsetPanel(
             tabPanel("Data",
                      div(
                        p("Import Data File"),
                        fileInput('InCell', label="InCell Output File"),
                        br(),
                        actionButton('console', label="Show Console"),
                        class="well tab-well"
                      )),
             tabPanel("Next",
                      div(
                        p("Next Tab"),
                        class="well tab-well")),
             position='left'
           )),
    column(7,
           h5("Cell Number"),
           plotOutput('heatmap'))
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
