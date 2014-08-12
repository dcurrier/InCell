
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
source("Parse_REMP_PlateLookup.R")
source("Analysis_Helpers.R")

noDataPlot = function(){
  par(mar=c(0,0,0,0))
  plot(1,1, type="n", bty="n", axes=F, main = "", ylab="", xlab="")
  text(1,1 , labels="No Data to Display", cex=2, col="#919191")
}

shinyServer(function(input, output, session) {

  ############### Reactives ###############

  # Store InCell Data as a list
  InCell = reactive({
    if( !is.null(input$InCell) ){
      withProgress(session, min=0, max=1, {
        setProgress(message='Parsing InCell File')
        data = ReadInCell(input$InCell[1,'datapath'], progressBar=TRUE )

        # Update Selectize Inputs
        updateSelectizeInput(session, 'featureCol',
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

  # Store compound annotations as a table
  REMP = reactive({
    if( !is.null(input$annotation) ){
      # Get the datapath
      path = input$annotation[1, 'datapath']

      # Parse the control column input
      ctlString = input$ctlCols
      if( ctlString == "example: 21-24 or 21,22,23,24" ){
        ctl = 21:24
      }else{
        ctlString = gsub("[[:lower:]]", "", ctlString)
        ctlString = gsub("[[:upper:]]", "", ctlString)
        ctlString = gsub("[[:digit:]] [[:digit:]]", ",", ctlString)
        ctlString = gsub("-", ":", ctlString)
        ctlString = gsub(" ", "", ctlString)
        ctl = eval(parse( text=paste0("c(", ctlString, ")") ))
      }

      # Parse the annotation file
      t = Parse_REMP(path, ctl)

      # Update selectize input for compound names
      updateSelectizeInput(session, 'cmpdSelect',
                           choices = unique(as.character(t$SAMPLE[which(t$SAMPLE != "DMSO")])),
                           selected = unique(as.character(t$SAMPLE[which(t$SAMPLE != "DMSO")]))[1] )

      return(t)
    }
  })

  # Neg Control Distributions
  negCtrlDist = reactive({
    if( is.null(InCell()) || is.null(REMP()) ) return()

    withProgress(session, min=0, max=1, {
      setProgress(message='Calculating Negative Control Distributions')

      # Get the list of negative control wells
      negCtrlWells = REMP()$well[which(REMP()$SAMPLE == "DMSO")]

      setProgress(value=0.01)
      negStats=mapply(function(feature, name){

        # Tabulate the ks and cvmts statistics for each well
        t=mapply(function(well){
          x = feature[[well]]
          y = unlist(
            mapply(function(negWell, well){
              if(negWell != well) feature[[negWell]]
            }, as.character(negCtrlWells), well, SIMPLIFY=T, USE.NAMES=F))

          list(ks=suppressWarnings(f_ks(x,y)), cvmts=f_cvmts(x,y))
        }, as.character(negCtrlWells), SIMPLIFY=T, USE.NAMES=T)

        # Set progress bar
        setProgress(value=which(names(InCell()$cellList) == name)/length(names(InCell()$cellList)), detail=name)

        # Make list of negCtl Statistics for each feature
        list( table=t,
              ks.median=median( unlist(t["ks", ]) ),
              ks.mad=mad( unlist(t["ks", ]) ),
              ks.mean=mean( unlist(t["ks", ]) ),
              ks.sd=sd( unlist(t["ks", ]) ),
              cvmts.median=median( unlist(t["cvmts", ]) ),
              cvmts.mean=mean( unlist(t["cvmts", ]) ),
              cvmts.mad=mad( unlist(t["cvmts", ]) ),
              cvmts.sd=sd( unlist(t["cvmts", ]) )
              )

      }, InCell()$cellList, names(InCell()$cellList), SIMPLIFY=F, USE.NAMES=T)

      setProgress(value=0.99, detail="wrapping up")

      return(negStats)
    })
  })

  doseCurveStats = reactive({
    if( is.null(REMP()) || is.null(input$cmpdSelect) || is.null(InCell()) ||
          is.null(input$featureCol) || is.null(negCtrlDist()) )  return()

    # Get the list of the wells that the selected compound is in
    conc = as.numeric(REMP()$CONC[which(REMP()$SAMPLE == input$cmpdSelect)])
    concOrder = sort.int(conc, index.return=T)
    wells = as.character(REMP()$well[which(REMP()$SAMPLE == input$cmpdSelect)])[concOrder$ix]

    # Get the list of negative control wells
    negCtrlWells = REMP()$well[which(REMP()$SAMPLE == "DMSO")]

    y = unlist(
      mapply(function(negWell){
        InCell()$cellList[[input$featureCol]][[negWell]]
      }, as.character(negCtrlWells), SIMPLIFY=T, USE.NAMES=F))

    d=mapply(function(w, c){
      # Get the values for the stat tests
      x = InCell()$cellList[[input$featureCol]][[w]]
      if( length(x) < 1 ){
        ks=NA
        cvmts=NA
      }else{
        # Calculate Statistics
        ks=suppressWarnings(f_ks(x,y))
        cvmts=f_cvmts(x,y)
      }
      # Z Transform
      ks.z=(ks - negCtrlDist()[[input$featureCol]]$ks.median)/negCtrlDist()[[input$featureCol]]$ks.mad
      cvmts.z=(cvmts - negCtrlDist()[[input$featureCol]]$cvmts.median)/negCtrlDist()[[input$featureCol]]$cvmts.mad

      # Calculate z transformed cell number
      DMSO = which(InCell()$well$Well %in% REMP()$well[which(REMP()$SAMPLE == "DMSO")])
      current = which(InCell()$well$Well %in% w)
      cn.z = (InCell()$well$`Cell Count`[current] - median(InCell()$well$`Cell Count`[DMSO]))/mad(InCell()$well$`Cell Count`[DMSO])

      # store stats
      list(ks=ks, ks.z=ks.z, cvmts=cvmts, cvmts.z=cvmts.z, conc=c, cellNum=cn.z )
    }, wells, concOrder$x, SIMPLIFY=T, USE.NAMES=T)

    return(d)
  })




  ############### Outputs ###############

  ## Data Tab ##

  # Tests for a file upload
  output$fileUploaded <- reactive({
    if( is.null(InCell()) || is.null(REMP()) ){
      return(T)
    }else{
      return(F)
    }

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
    v = eval(parse(text=paste0('InCell()$well$`',input$featureCol, '`')))
    stats = as.list(summary(v))

    # make output
    paste0(
      "Summary Statistics for\n",
      input$featureCol, "\n",
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
      v = eval(parse(text=paste0('InCell()$well$`',input$featureCol, '`')))
      par(mar=c(3, 4, 6, 2)+0.1, bg=rgb(0,0,0,0), fg="#333333")
      hist(v,
           col="#58849e",
           main=paste0("Histogram of\n", input$featureCol),
           xlab="",
           ylab="Freq",
           las=1,
           border="#BBBBBB",
           col.axis="#333333",
           col.lab="#333333",
           col.main="#000000")
    }else{
      par(mar=c(5, 4, 4, 2)+0.1, bg=rgb(0,0,0,0))
    }
  })


  ## Field Tab ##
  output$miniFieldData = renderPlot({
    if( !(is.null(input$InCell)) && !(is.null(InCell()$field)) ){
      if(input$featureCol != "Cell Count"){
        d = eval(parse(text=paste("InCell()$cell$`", input$featureCol,"`", sep="")))
        fieldNames = as.character(unique( InCell()$field$Well[grep(input$fieldWell, InCell()$field$Well )] ))

        l = mapply(function(field){
          index = grep(field, InCell()$cell$Well)
          eval(parse(text=paste("InCell()$cell$`", input$featureCol,"`[index]", sep="")))
        }, fieldNames, SIMPLIFY=F, USE.NAMES=F)
        names(l) = fieldNames



        par(mar=c(5, 4, 4, 2)+0.1, bg=rgb(0,0,0,0))
        stripchart(l, vertical=T, method="jitter", jitter=1, pch=16, col="#58849e",
                   las=1, main=input$fieldColumn)
        axis(1, at=1:length(l), label=names(l))
      }else{
        par(mar=c(5, 4, 4, 2)+0.1, bg=rgb(0,0,0,0))
      }
    }else{
      par(mar=c(5, 4, 4, 2)+0.1, bg=rgb(0,0,0,0))
    }
  })

  output$fieldSummary = renderText({
    paste0(
      "Number of Cells: ",
      "Add calculation"
    )
  })



  ## Feature Tab ##
  output$concSlide = renderUI({
    if( !is.null(REMP()) && !is.null(input$cmpdSelect) && input$cmpdSelect != ""){
      count = length(REMP()$CONC[which(REMP()$SAMPLE == input$cmpdSelect)])
      sliderInput('conc', label=h4("Concentration"),
                  value=ceiling(count/2), min=1, max=count, step=1)
    }
  })

  output$selConc = renderText({
    if( !is.null(input$conc) ){
      value = as.numeric(REMP()$CONC[which(REMP()$SAMPLE == input$cmpdSelect)][input$conc])
      paste0(round(value, digits=0), " nM")
    }
  })

  output$ksStatSummary = renderText({
    if( !is.null(REMP()) && !is.null(input$cmpdSelect) &&
          !is.null(negCtrlDist()) && input$featureCol != "Cell Count" ){
      paste0(
        "K-S Test\n",
        "Median: ", round(negCtrlDist()[[input$featureCol]]$ks.median, digits=2), "\n",
        "MAD: ", round(negCtrlDist()[[input$featureCol]]$ks.mad, digits=2), "\n",
        "Mean: ", round(negCtrlDist()[[input$featureCol]]$ks.mean, digits=2), "\n",
        "SD: ", round(negCtrlDist()[[input$featureCol]]$ks.sd, digits=2)
        )
    }else{
      paste0(  "    Statistics\n\n   not available\n"  )
    }
  })

  output$cvmtsStatSummary = renderText({
    if( !is.null(REMP()) && !is.null(input$cmpdSelect) &&
          !is.null(negCtrlDist()) && input$featureCol != "Cell Count" ){
      paste0(
        "CVMTS Test\n",
        "Median: ", round(negCtrlDist()[[input$featureCol]]$cvmts.median, digits=2), "\n",
        "MAD: ", round(negCtrlDist()[[input$featureCol]]$cvmts.mad, digits=2), "\n",
        "Mean: ", round(negCtrlDist()[[input$featureCol]]$cvmts.mean, digits=2), "\n",
        "SD: ", round(negCtrlDist()[[input$featureCol]]$cvmts.sd, digits=2)
      )
    }else{
      paste0(  "    Statistics\n\n   not available\n"  )
    }
  })

  output$cmpdStatSummary = renderText({
    if( !is.null(InCell()) && !is.null(REMP()) && !is.null(negCtrlDist())
        && input$featureCol != 'Cell Count' ){

      # Get the list of negative control wells
      negCtrlWells = REMP()$well[which(REMP()$SAMPLE == "DMSO")]


      x = InCell()$cellList[[input$featureCol]][[input$fieldWell]]
      if( length(x) < 1 ){
        # Print stats
        paste0(
          "K-S: NA\n",
          "z: NA\n\n",
          "CVMTS: NA\n",
          "z: NA\n\n",
          "n: NA"
        )
      }else{
        y = unlist(
          mapply(function(negWell, well){
            if(negWell != well) InCell()$cellList[[input$featureCol]][[negWell]]
          }, as.character(negCtrlWells), input$fieldWell, SIMPLIFY=T, USE.NAMES=F))

        #Calculate the statistics
        ks=suppressWarnings(f_ks(x,y))
        cvmts=f_cvmts(x,y)

        # Print stats
        paste0(
          "K-S: ", round(ks, 2), "\n",
          "z: ", round((ks - negCtrlDist()[[input$featureCol]]$ks.median)/negCtrlDist()[[input$featureCol]]$ks.mad, 2), "\n\n",
          "CVMTS: ", round(cvmts, 2), "\n",
          "z: ", round((cvmts - negCtrlDist()[[input$featureCol]]$cvmts.median)/negCtrlDist()[[input$featureCol]]$cvmts.mad, 2), "\n\n",
          "n: ", length(x)
          )
      }
    }else{
      paste0(  "    Statistics\n   not available\n"  )
    }
  })

  output$cmpdZplot = renderPlot({
    if( is.null(doseCurveStats()) || input$featureCol == "Cell Count" ) return()

    # get the z stat vectors
    ks.z = unlist(doseCurveStats()['ks.z', ])
    cvmts.z = unlist(doseCurveStats()['cvmts.z', ])
    conc = unlist(doseCurveStats()['conc', ])
    cells = unlist(doseCurveStats()['cellNum', ])

    par(mar=c(4,2,2,0)+0.1, bg=rgb(0,0,0,0))

    plot(x=log10(conc/1000), y=ks.z,
         las=1,
         type="l",
         bty="n",
         lwd=3,
         ylim=c(min(c(ks.z,cvmts.z,cells), na.rm=T),max(c(ks.z,cvmts.z,cells), na.rm=T)),
         ylab="",
         xlab="log10 Conc [nM]",
         col="#2A6373"
    )

    points(x=log10(conc/1000), y=cvmts.z,
           type="l", lwd=3, col="#DD8702"
           )

    points(x=log10(conc/1000), y=cells,
           type="l", lwd=3, col="#C1C1C1" )

  })



  ## Main Panel ##

  # Main Heatmap data output
  output$heatmap = renderPlot({
    if( !(is.null(input$InCell)) && !(is.null(InCell()$well)) ){
      #browser()
      d = matrix(eval(parse(text=paste("InCell()$well$`", input$featureCol,"`", sep=""))),
                 nrow=16, ncol=24, byrow=T)
      par(mar=c(6,2,6,0), bg=rgb(0,0,0,0))
      color2D.matplot(d,
                      show.values = FALSE,
                      show.legend = TRUE,
                      axes = FALSE,
                      main = input$featureCol,
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
    }else{
      noDataPlot()
    }
  })
  outputOptions(output, 'heatmap', suspendWhenHidden=FALSE)

  # Plot of each field
  output$fieldwisePlot = renderPlot({
    if( !(is.null(input$InCell)) && !(is.null(InCell()$field)) ){
      #browser()
      d = eval(parse(text=paste("InCell()$field$`", input$featureCol,"`", sep="")))

      l = mapply(function(well){
        index = grep(well, InCell()$field$Well)
        return(d[index])
      }, as.character(InCell()$well$Well), SIMPLIFY=F, USE.NAMES=F)
      names(l) = as.character(InCell()$well$Well)

      chosen = grep(paste0(input$fieldWell, "$"), names(l))

      if(max(d) > 1) par(mar=c(6,2,6,0), bg=rgb(0,0,0,0))
      if(max(d) > 100) par(mar=c(6,3,6,0), bg=rgb(0,0,0,0))
      if(max(d) > 1000) par(mar=c(6,4,6,0), bg=rgb(0,0,0,0))
      stripchart(l, vertical = T, col="#58849e",
                 pch=16, las=1, axes=F,
                 main=input$featureCol)
      points(chosen, l[chosen], pch=16, cex=2, col="#bd0000")
      axis(2, las=1)
      axis(1, at=c(seq(1, 384, by=24), 384), labels=names(l)[c(seq(1, 384, by=24), 384)])
      box()
    }else{
      noDataPlot()
    }
  })
  outputOptions(output, 'fieldwisePlot', suspendWhenHidden=FALSE)



  # Plot of the feature of interest distribution
  output$featureDistPlot = renderPlot({
    if( !is.null(InCell()) && !is.null(REMP()) && !is.null(negCtrlDist())
        && !is.null(input$conc) && input$featureCol != 'Cell Count' ){

      # Get the list of negative control wells
      negCtrlWells = REMP()$well[which(REMP()$SAMPLE == "DMSO")]

      # Get the current Compount Conc
      value = as.numeric(REMP()$CONC[which(REMP()$SAMPLE == input$cmpdSelect)][input$conc])


      x = InCell()$cellList[[input$featureCol]][[input$fieldWell]]
      if( length(x) < 1 ){
        noDataPlot()
      }else{
        y = unlist(
          mapply(function(negWell, well){
            if(negWell != well) InCell()$cellList[[input$featureCol]][[negWell]]
          }, as.character(negCtrlWells), input$fieldWell, SIMPLIFY=T, USE.NAMES=F))

        #Calculate the statistics
        ks=suppressWarnings(f_ks(x,y))
        cvmts=f_cvmts(x,y)

        # Calculate breaks
        breaks=hist(c(x,y), plot = FALSE, breaks = 100)$breaks

        # Set up a multiplot layout
        layout( matrix(c(1:2), nrow=2, ncol=1), heights=c(7,1) )
        par(mar=c(0,4,4,2)+0.1)

        # Draw Histogram of DMSO Wells
        hist( y,
              freq=F,
              breaks=breaks,
              lwd = 0.5,
              ylab = "", xlab = "",
              main=paste0(input$fieldWell, " [", round(value, digits=0), " nM]"),
              cex.main=1.5,
              axes = FALSE,
              ylim = c(0, max(c(pretty(density(y)$y), pretty(density(x)$y)))),
              col="#FFE8C7",
              border="#FFD599"
        )

        hist( x,
              freq=F,
              breaks=breaks,
              add=T,
              col="#9CCCEB",
              border="#6BB9E3")

        lines(density(y), lwd = 2, col = "#EBAF59")
        lines(density(x), lwd = 2, col = "#256E9E")

        if( length(x) < 30 ){
          box(which="plot", lty="solid", lwd=3, col="red")
          legend("topright", legend=paste0("n=",length(x)),
                 pch="", text.col="red", cex=2, bty="n")
        }

        # Plot the legend
        par(mar=c(0,0,0,0))
        frame()

        legend("center", legend=c("DMSO", input$cmpdSelect),
               col=c("#EBAF59", "#256E9E"), lwd=2,
               cex=1.5, ncol=2, bty="n")
      }
    }else{
      noDataPlot()
    }
  })
  outputOptions(output, 'featureDistPlot', suspendWhenHidden=FALSE)


  output$featureDistGrid = renderPlot({
    if( is.null(InCell()) && is.null(REMP()) && is.null(negCtrlDist())
        && is.null(input$conc) && input$featureCol == 'Cell Count' ) noDataPlot()

    # Get the list of negative control wells
    negCtrlWells = REMP()$well[which(REMP()$SAMPLE == "DMSO")]

    # Get the list of the wells that the selected compound is in
    conc = as.numeric(REMP()$CONC[which(REMP()$SAMPLE == input$cmpdSelect)])
    concOrder = sort.int(conc, index.return=T)
    wells = as.character(REMP()$well[which(REMP()$SAMPLE == input$cmpdSelect)])[concOrder$ix]

    # Get the negative control well values
    y = unlist(
      mapply(function(negWell){
        InCell()$cellList[[input$featureCol]][[negWell]]
      }, as.character(negCtrlWells), SIMPLIFY=T, USE.NAMES=F))

    # Setup the layout
    n = length(wells)
    if(n<6){
      # Calculate Matrix Dimensions
      cols = 3
      rows = 1
    }else if(n>6 && n<12){
      # Calculate Matrix Dimensions
      cols = ceiling(n/2)
      rows = ceiling(n/cols)
    }else if(n>12){
      # Calculate Matrix Dimensions
      rows = ceiling(n/6)
      cols = 6
    }

    # Construct Matrix
    mat = matrix( c(1:(rows*cols)), nrow=rows, ncol=cols, byrow=T )
    mat = rbind(mat, rep(max(mat)+1, cols))

    # Setup the heights
    heights = c(rep(3, rows), 1)

    # Make the layout
    layout(mat, heights=heights)

    # Generate the plots
    for( i in 1:length(wells) ){
      # Get the values for the stat tests
      x = InCell()$cellList[[ input$featureCol ]][[ wells[i] ]]

      if( length(x) < 1 ){
        frame()
      }else{
        #Calculate the statistics
        ks=suppressWarnings(f_ks(x,y))
        cvmts=f_cvmts(x,y)

        # Calculate breaks
        breaks=hist(c(x,y), plot = FALSE, breaks = 100)$breaks

        # Set margins
        if( wells[i] == input$fieldWell ) {
          #browser()
          par(mar=c(0,4,8,2)+0.1, bg=rgb(223, 238, 255, 100, maxColorValue = 255) )
        }else{
          par(mar=c(0,4,8,2)+0.1)
        }

        # Draw Histogram of DMSO Wells
        hist( y,
              freq=F,
              breaks=breaks,
              lwd = 0.5,
              ylab = "", xlab = "",
              main=paste0(wells[i], " [", round(concOrder$x[i], digits=0), " nM]"),
              cex.main=2,
              axes = FALSE,
              ylim = c(0, max(c(pretty(density(y)$y), pretty(density(x)$y)))),
              col="#FFE8C7",
              border="#FFD599"
        )

        hist( x,
              freq=F,
              breaks=breaks,
              add=T,
              col="#9CCCEB",
              border="#6BB9E3")

        lines(density(y), lwd = 2, col = "#EBAF59")
        lines(density(x), lwd = 2, col = "#256E9E")

        if( length(x) < 30 ){
          box(which="plot", lty="solid", lwd=2, col="red")
          legend("topright", legend=paste0("n=",length(x)),
                 pch="", text.col="red", cex=1.5, bty="n")
        }
      }
    }

      if( length(wells) < (max(mat)-1) ){
        for( i in (length(wells)+1):(max(mat)-1) ){
          frame()
        }
    }

    # Generate legend
    par(mar=c(0,0,0,0))
    frame()

    legend("center", legend=c("DMSO", input$cmpdSelect),
           col=c("#EBAF59", "#256E9E"), lwd=2,
           cex=1.5, ncol=2, bty="n")
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

  observe({
    updateSelectizeInput(session, 'featureColumn', selected=input$wellColumn)
  })

   observe({
     if( !is.null(InCell()) && !is.null(REMP()) && !is.null(input$conc) &&
           input$cmpdSelect != ""){
       # Get the well coordinates
       cmpd = input$cmpdSelect
       well = as.character(REMP()$well[which(REMP()$SAMPLE == cmpd)[input$conc]])

       updateSelectizeInput(session, 'fieldWell', selected=well)
     }
   })

  if( TRUE ){
    observe({
      if( input$prv > 0 ){
        isolate({
          cList = unique(as.character(REMP()$SAMPLE[which(REMP()$SAMPLE != "DMSO")]))
          n = which(cList == input$cmpdSelect[1])-1
          if( n <= 0 ) {n = 1}
          newName = cList[n]
        })
        #browser()
        updateSelectizeInput(session, 'cmpdSelect', selected=newName )
      }
    }) # Previous

    observe({
      if( input$nxt > 0 ){
        isolate({
          cList = unique(as.character(REMP()$SAMPLE[which(REMP()$SAMPLE != "DMSO")]))
          n = which(cList == input$cmpdSelect[1])+1
          if( n > length(cList) ) {n = length(cList)}
          newName = cList[n]
        })
        #browser()
        updateSelectizeInput(session, 'cmpdSelect', selected=newName )
      }
    }) # Next

  }  # Previous/Next compound navigation

})
