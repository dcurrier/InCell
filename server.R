
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)             # Shiny Framework
library(shinyIncubator)    # Progress indicator
library(shinythings)       # Password Input/Better Action buttons
library(plotrix)           # Heatmap
source("Parse_InCell.R")

shinyServer(function(input, output, session) {

  # Store InCell Data as a list
  InCell = reactive({
    if(!(is.null(input$InCell)) ){
      withProgress(session, min=0, max=1, {
        setProgress(message='Parsing InCell File')
        ReadInCell(input$InCell[1,'datapath'], progressBar=TRUE )
      })
    }
  })

  output$heatmap = renderPlot({
    if( !(is.null(input$InCell)) && !(is.null(InCell()$well)) ){
      #browser()
      d = matrix(InCell()$well$`Cell Count`, nrow=16, ncol=24, byrow=T)
      par(mar=c(4,0,0,0))
      color2D.matplot(d,
                      show.values = FALSE,
                      show.legend = TRUE,
                      axes = FALSE,
                      xlab = "",
                      ylab = "",
                      vcex = 1,
                      vcol = "#FF3526",
                      extremes = c("#FFFFFF", "#00436B") # http://adobe.ly/1nMjd54
        )

    }



  })



  # Debug
  observe(label="console",{
    if(input$console != 0) {
      options(browserNLdisabled=TRUE)
      saved_console<-".RDuetConsole"
      if (file.exists(saved_console)) load(saved_console)
      isolate(browser())
      save(file=saved_console,list=ls(environment()))
    }
  })

})
