
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

  ############### Reactives ###############

  # Store InCell Data as a list
  InCell = reactive({
    if(!(is.null(input$InCell)) ){
      withProgress(session, min=0, max=1, {
        setProgress(message='Parsing InCell File')
        data = ReadInCell(input$InCell[1,'datapath'], progressBar=TRUE )

        # Update Selectize Inputs
        updateSelectizeInput(session, 'wellColumn',
                             choices = names(data$well)[which(names(data$well) != "Well")],
                             selected="Cell Count")

        updateSelectizeInput(session, 'fieldColumn',
                             choices = names(data$well)[which(names(data$well) != "Well")],
                             selected="Cell Count")

        updateSelectizeInput(session, 'fieldWell',
                             choices = as.character(data$well$Well),
                             selected=as.character(data$well$Well[1]))

        return(data)
      })
    }else{
      return(NULL)
    }
  })





  ############### Outputs ###############

  ## Data Tab ##

  # Tests for a file upload
  output$fileUploaded <- reactive({
    return(is.null(InCell()))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)

  # Statistics secton of Data tab
  output$dataset = renderText({
    paste0(
      "Data File:\n",
      " ",input$InCell[1,'name'], "\n\n",
      "Total Wells: ", dim(InCell()$well)[1], "\n",
      "Total Cells Counted: ", format(dim(InCell()$cell)[1], big.mark=","), "\n",
      "Total Fields: ", dim(InCell()$field)[1], "\n",
      "Mean Fields Per Well: ", dim(InCell()$field)[1]/dim(InCell()$well)[1], "\n\n"
    )
  })

  # Download Handler for Well Level Data
  output$downWell = downloadHandler(
    filename = function(){
      name = strsplit(input$InCell[1,'name'], "[.]")[[1]]
      paste(paste(name[1:length(name)-1], collapse="."), "-Well.csv", sep = "")
    },
    content = function(file){
      write.csv(InCell()$well, file=file, row.names=F)
    }
    )

  # Download Handler for Field Level Data
  output$downField = downloadHandler(
    filename = function(){
      name = strsplit(input$InCell[1,'name'], "[.]")[[1]]
      paste(paste(name[1:length(name)-1], collapse="."), "-Field.csv", sep = "")
    },
    content = function(file){
      write.csv(InCell()$field, file=file, row.names=F)
    }
  )

  # Download Handler for Cell Level Data
  output$downCell = downloadHandler(
    filename = function(){
      name = strsplit(input$InCell[1,'name'], "[.]")[[1]]
      paste(paste(name[1:length(name)-1], collapse="."), "-Cell.csv", sep = "")
    },
    content = function(file){
      write.csv(InCell()$cell, file=file, row.names=F)
    }
  )


  ## Well Tab ##

  # Statistics section of Well tab
  output$wellData = renderText({
    # Generate parameter distribution statistics
    v = eval(parse(text=paste0('InCell()$well$`',input$wellColumn, '`')))
    stats = as.list(summary(v))

    # make output
    paste0(
      "Summary Statistics for\n",
      input$wellColumn, "\n",
      " Min: ", stats$Min, "\n",
      " 1st Quartile: ", stats$`1st Qu.`, "\n",
      " Median: ", stats$Median, "\n",
      " Mean: ", stats$Mean, "\n",
      " 3rd Quartile: ", stats$`3rd Qu.`, "\n",
      " Max: ", stats$Max, "\n"
    )
  })

  # Mini Histogram in Well tab
  output$miniHist = renderPlot({
    if( !is.null(InCell()$well) ){
      v = eval(parse(text=paste0('InCell()$well$`',input$wellColumn, '`')))
      par(mar=c(3, 4, 6, 2)+0.1, bg=rgb(0,0,0,0), fg="#333333")
      hist(v,
           col="#58849e",
           main=paste0("Histogram of\n", input$wellColumn),
           xlab="",
           ylab="Freq",
           las=1,
           border="#BBBBBB",
           col.axis="#333333",
           col.lab="#333333",
           col.main="#000000")
    }
  })



  ## Main Panel ##

  # Main Heatmap data output
  output$heatmap = renderPlot({
    if( !(is.null(input$InCell)) && !(is.null(InCell()$well)) ){
      #browser()
      d = matrix(eval(parse(text=paste("InCell()$well$`", input$wellColumn,"`", sep=""))),
                 nrow=16, ncol=24, byrow=T)
      par(mar=c(6,2,6,0), bg=rgb(0,0,0,0))
      color2D.matplot(d,
                      show.values = FALSE,
                      show.legend = TRUE,
                      axes = FALSE,
                      main = input$wellColumn,
                      xlab = "",
                      ylab = "",
                      vcex = 1,
                      vcol = "#FF3526",
                      extremes = c("#FFFFFF", "#00436B") # http://adobe.ly/1nMjd54
      )
      axis(3, at = seq_len(ncol(d)) - 0.5,
           labels = as.character(1:24), tick = FALSE)
      axis(2, at = seq_len(nrow(d)) -0.5,
           labels = LETTERS[16:1], tick = FALSE, las = 1)
    }
  })

  output$fieldwisePlot = renderPlot({
    if( !(is.null(input$InCell)) && !(is.null(InCell()$field)) ){
      #browser()
      d = eval(parse(text=paste("InCell()$field$`", input$fieldColumn,"`", sep="")))

      l = mapply(function(well){
        index = grep(paste0(well, "[[:punct:]]"), InCell()$field$Well)
        return(d[index])
      }, as.character(InCell()$well$Well), SIMPLIFY=F, USE.NAMES=F)
      names(l) = as.character(InCell()$well$Well)

      chosen = grep(paste0(input$fieldWell, "$"), names(l))

      if(max(d) > 1) par(mar=c(6,2,6,0), bg=rgb(0,0,0,0))
      if(max(d) > 100) par(mar=c(6,3,6,0), bg=rgb(0,0,0,0))
      if(max(d) > 1000) par(mar=c(6,4,6,0), bg=rgb(0,0,0,0))
      stripchart(l, vertical = T, col="#58849e",
                 pch=16, las=1, axes=F,
                 main=input$fieldColumn)
      points(chosen, l[chosen], pch=16, cex=2, col="#bd0000")
      axis(2, las=1)
      axis(1, at=c(seq(1, 384, by=24), 384), labels=names(l)[c(seq(1, 384, by=24), 384)])
      box()
    }
  })
  outputOptions(output, 'fieldwisePlot', suspendWhenHidden=FALSE)

  output$miniFieldData = renderPlot({
    if( !(is.null(input$InCell)) && !(is.null(InCell()$field)) ){
      if(input$fieldColumn != "Cell Count"){
        d = eval(parse(text=paste("InCell()$cell$`", input$fieldColumn,"`", sep="")))
        fieldNames = as.character(unique( InCell()$field$Well[grep( paste0(input$fieldWell, "[[:punct:]]"), InCell()$field$Well )] ))

        l = mapply(function(field){
          field = gsub("\\(", "\\\\(", field)
          field = gsub("\\)", "\\\\)", field)
          index = grep(field, InCell()$cell$Well)
          eval(parse(text=paste("InCell()$cell$`", input$fieldColumn,"`[index]", sep="")))
        }, fieldNames, SIMPLIFY=F, USE.NAMES=F)
        names(l) = fieldNames



        par(mar=c(5, 4, 4, 2)+0.1, bg=rgb(0,0,0,0))
        stripchart(l, vertical=T, method="jitter", jitter=1, pch=16, col="#58849e",
                   las=1, main=input$fieldColumn)
        axis(1, at=1:length(l), label=names(l))
      }else{
        par(mar=c(5, 4, 4, 2)+0.1, bg=rgb(0,0,0,0))
      }
    }
  })

  output$fieldSummary = renderText({
    paste0(
      "Number of Cells: "

      )
  })

  ############### Observers ###############

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

  observe({
    updateSelectizeInput(session, 'fieldColumn', selected=input$wellColumn)
  })

  observe({
    updateSelectizeInput(session, 'wellColumn', selected=input$fieldColumn)
  })

})
