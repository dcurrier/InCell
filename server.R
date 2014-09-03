
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

require(shiny)             # Shiny Framework
require(shinyIncubator)    # Progress indicator
require(shinythings)       # Password Input/Better Action buttons
require(MESS)               # AUC Analysis
require(ShinyHighCharts)   # Javascript charting


noDataPlot = function(){
  par(mar=c(0,0,0,0))
  plot(1,1, type="n", bty="n", axes=F, main = "", ylab="", xlab="")
  text(1,1 , labels="No Data to Display", cex=2, col="#919191")
}

getDropboxFiles = function(){
  if( file.exists("Y:/instrument/HTS Instruments QC and Study/InCell 6000/ShinyDropbox") ){
    path = "Y:/instrument/HTS Instruments QC and Study/InCell 6000/ShinyDropbox"
  }

  if( file.exists("Z:/instrument/HTS Instruments QC and Study/InCell 6000/ShinyDropbox") ){
    path = "Z:/instrument/HTS Instruments QC and Study/InCell 6000/ShinyDropbox"
  }

  if( file.exists( "/media/sf_ShinyDropbox" ) ){
    path = "/media/sf_ShinyDropbox"
  }

  # Get list of files in the folder
  allFiles = list.files(path)

  # Get list of csv Files
  csvFiles = mapply(function(f){
    splt = strsplit(f, "[.]")[[1]]
    type = splt[length(splt)]
    if( type %in% c("csv", "CSV") ){
      info = file.info(paste0(path, "/", f))
      list(path=paste0(path, "/", f), size=info[1, 'size']/1000000)
    }
  },allFiles, SIMPLIFY=F, USE.NAMES=T)
  csvFiles[sapply(csvFiles, is.null)] <- NULL

  #get List of txt Files
  txtFiles = mapply(function(f){
    splt = strsplit(f, "[.]")[[1]]
    type = splt[length(splt)]
    if( type %in% c("txt", "TXT") ){
      info = file.info(paste0(path, "/", f))
      list(path=paste0(path, "/", f), size=info[1, 'size']/1000000)
    }
  },allFiles, SIMPLIFY=F, USE.NAMES=T)
  txtFiles[sapply(txtFiles, is.null)] <- NULL

  # Return the lists of files
  list(csv=csvFiles, txt=txtFiles)
}

makeColorCode = function(min, max, data){
  # Test for min, max types
  if( typeof(min) != "character" ) stop("'min' must be supplied as a character vector")
  if( typeof(max) != "character" ) stop("'max' must be supplied as a character vector")

  if( nchar(min) %in% c(7,9) && strsplit(min, "")[[1]][1] == "#" &&
      nchar(max) %in% c(7,9) && strsplit(max, "")[[1]][1] == "#" ){
    # Scale data range onto 0-1 scale
    range01 <- function(x)(x-min(x))/diff(range(x))
    # generate colors
    cols <- colorRamp(c(min, max))(range01(data))
    # convert colors back to hex
    apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
  }else{
    stop("'min' and 'max' must be supplied as hex color codes with preceeding '#'")
  }
}

shinyServer(function(input, output, session) {



  ############### Reactives ###############

  # Get Dropbox File List
  Dropbox = reactive({
    input$updateDropbox
    getDropboxFiles()
  })





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
        updateSelectizeInput(session, 'featureColDist',
                             choices = names(data$well)[which(names(data$well) != "Well")],
                             selected="Cell Count")
        updateSelectizeInput(session, 'yFeat',
                             choices = names(data$well)[which(names(data$well) != "Well")],
                             selected="Cell Count")
        updateSelectizeInput(session, 'xFeat',
                             choices = names(data$well)[which(names(data$well) != "Well")],
                             selected="Cell Count")
        updateSelectizeInput(session, 'fieldWell',
                             choices = as.character(data$well$Well),
                             selected=as.character(data$well$Well[1]))
        # Return
        return(data)
      })
    }else if( !is.null(input$InCellDB) && input$InCellDB != "Choose from Dropbox" ){
      withProgress(session, min=0, max=1, {
        setProgress(message='Parsing InCell File')
        data = ReadInCell(Dropbox()$csv[[input$InCellDB]]$path, progressBar=TRUE )

        # Update Selectize Inputs
        updateSelectizeInput(session, 'featureCol',
                             choices = names(data$well)[which(names(data$well) != "Well")],
                             selected="Cell Count")
        updateSelectizeInput(session, 'featureColDist',
                             choices = names(data$well)[which(names(data$well) != "Well")],
                             selected="Cell Count")
        updateSelectizeInput(session, 'yFeat',
                             choices = names(data$well)[which(names(data$well) != "Well")],
                             selected="Cell Count")
        updateSelectizeInput(session, 'xFeat',
                             choices = names(data$well)[which(names(data$well) != "Well")],
                             selected="Cell Count")
        updateSelectizeInput(session, 'fieldWell',
                             choices = as.character(data$well$Well),
                             selected=as.character(data$well$Well[1]))
        # Return
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
    }else if( !is.null(input$annotationDB) && input$annotationDB != "Choose from Dropbox" ){
      # Get the datapath
      path = Dropbox()$txt[[input$annotationDB]]$path

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
    }else{
      NULL
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
      negStats=mapply(function(name){
        feature=InCell()$cellList[[name]]

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

      }, names(InCell()$cellList), SIMPLIFY=F, USE.NAMES=T)

      setProgress(value=0.99, detail="wrapping up")

      return(negStats)
    })
  })

  doseCurveStats = reactive({
    if( is.null(REMP()) || is.null(input$cmpdSelect) || is.null(InCell()) ||
          (is.null(input$featureColDist) && !(input$featureColDist %in% names(InCell()$cellList)))
          || is.null(negCtrlDist()) )  return()

    # Get the list of the wells that the selected compound is in
    conc = as.numeric(REMP()$CONC[which(REMP()$SAMPLE == input$cmpdSelect)])
    concOrder = sort.int(conc, index.return=T)
    wells = as.character(REMP()$well[which(REMP()$SAMPLE == input$cmpdSelect)])[concOrder$ix]

    # Get the list of negative control wells
    negCtrlWells = REMP()$well[which(REMP()$SAMPLE == "DMSO")]

    y = unlist(
      mapply(function(negWell){
        InCell()$cellList[[input$featureColDist]][[negWell]]
      }, as.character(negCtrlWells), SIMPLIFY=T, USE.NAMES=F))

    d=mapply(function(w, c){
      # Get the values for the stat tests
      x = InCell()$cellList[[input$featureColDist]][[w]]
      if( length(x) < 1 ){
        ks=NA
        cvmts=NA
      }else{
        # Calculate Statistics
        ks=suppressWarnings(f_ks(x,y))
        cvmts=f_cvmts(x,y)
      }
      # Z Transform
      ks.z=(ks - negCtrlDist()[[input$featureColDist]]$ks.median)/negCtrlDist()[[input$featureColDist]]$ks.mad
      cvmts.z=(cvmts - negCtrlDist()[[input$featureColDist]]$cvmts.median)/negCtrlDist()[[input$featureColDist]]$cvmts.mad

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
    if( !is.null(input$InCell) ) ICname = input$InCell[1,'name']
    if( input$InCellDB != "InCell File" ) ICname = input$InCellDB

    if( !is.null(input$annotation) ) REMPname = input$annotation[1,'name']
    if( input$annotationDB != "InCell File" ) REMPname = input$annotationDB

    if( exists('ICname') && exists('REMPname') && !is.null(ICname) && !is.null(REMPname) ){
      paste0(
        "Data Files:\n",
        " ",ICname, "\n",
        " ",REMPname, "\n\n",
        "Total Wells: ", dim(InCell()$well)[1], "\n",
        "Total Cells Counted: ", format(dim(InCell()$cell)[1], big.mark=","), "\n",
        "Total Fields: ", dim(InCell()$field)[1], "\n",
        "Mean Fields Per Well: ", dim(InCell()$field)[1]/dim(InCell()$well)[1], "\n\n"
      )
    }
  })

  # Download Handler for Well Level Data
  output$downWell = downloadHandler(
    filename = function(){
      if( !is.null(input$InCell) ) name = strsplit(input$InCell[1,'name'], "[.]")[[1]]
      if( input$InCellDB != "InCell File" ) name = strsplit(input$InCellDB, "[.]")[[1]]
      paste(paste(name[1:length(name)-1], collapse="."), "-Well.csv", sep = "")
    },
    content = function(file){
      write.csv(InCell()$well, file=file, row.names=F)
    }
    )

  # Download Handler for Field Level Data
  output$downField = downloadHandler(
    filename = function(){
      if( !is.null(input$InCell) ) name = strsplit(input$InCell[1,'name'], "[.]")[[1]]
      if( input$InCellDB != "InCell File" ) name = strsplit(input$InCellDB, "[.]")[[1]]
      paste(paste(name[1:length(name)-1], collapse="."), "-Field.csv", sep = "")
    },
    content = function(file){
      write.csv(InCell()$field, file=file, row.names=F)
    }
  )

  # Download Handler for Cell Level Data
  output$downCell = downloadHandler(
    filename = function(){
      if( !is.null(input$InCell) ) name = strsplit(input$InCell[1,'name'], "[.]")[[1]]
      if( input$InCellDB != "InCell File" ) name = strsplit(input$InCellDB, "[.]")[[1]]
      paste(paste(name[1:length(name)-1], collapse="."), "-Cell.csv", sep = "")
    },
    content = function(file){
      write.csv(InCell()$cell, file=file, row.names=F)
    }
  )


  ## QC Tab ##

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



  ## Feature Tab ##
  output$statSummary = renderTable({
    if( !is.null(REMP()) && !is.null(input$cmpdSelect) &&
          !is.null(negCtrlDist()) && !is.null(input$featureColDist) &&
          input$featureColDist %in% names(InCell()$cellList) ){

      table = data.frame(KS=c(
                              round(negCtrlDist()[[input$featureColDist]]$ks.median, digits=2),
                              round(negCtrlDist()[[input$featureColDist]]$ks.mad, digits=2),
                              round(negCtrlDist()[[input$featureColDist]]$ks.mean, digits=2),
                              round(negCtrlDist()[[input$featureColDist]]$ks.sd, digits=2)
                              ),
                         CVMTS=c(
                              round(negCtrlDist()[[input$featureColDist]]$cvmts.median, digits=2),
                              round(negCtrlDist()[[input$featureColDist]]$cvmts.mad, digits=2),
                              round(negCtrlDist()[[input$featureColDist]]$cvmts.mean, digits=2),
                              round(negCtrlDist()[[input$featureColDist]]$cvmts.sd, digits=2)
                           ))
      row.names(table) = c("Median", "MAD", "Mean", "SD")

    }else{
      table = data.frame(KS=rep(NA, 4), CVMTS=rep(NA, 4))
      row.names(table) = c("Median", "MAD", "Mean", "SD")
    }

    return(table)
  })



  ## Any|Any Tab ##
  output$ySlider = renderUI({
    if( !is.null(InCell()) && !is.null(REMP()) && input$yThresh != "None" && input$yFeat != "Cell Count" ){

      # Calculate the range of data for this feature
      rng = range(unlist(InCell()$cellList[[input$yFeat]]), na.rm=T)

      # Generate slider
      sliderInput('ySlide', label="Threshold", min=floor(rng[1]), max=ceiling(rng[2]), value=mean(rng), step=1)
    }
  })

  output$xSlider = renderUI({
    if( !is.null(InCell()) && !is.null(REMP()) && input$xThresh != "None" && input$xFeat != "Cell Count" ){

      # Calculate the range of data for this feature
      rng = range(unlist(InCell()$cellList[[input$xFeat]]), na.rm=T)

      # Generate slider
      sliderInput('xSlide', label="Threshold", min=floor(rng[1]), max=ceiling(rng[2]), value=mean(rng), step=1)
    }
  })


  ############### Main Panel ###############
  ## QC Tab ##
  # Large Heatmap
  output$highHeat = renderHighcharts({
    if( !is.null(InCell()$well) ){
      rows = rep(c(0:15), 24)

      cols = vector()
      for(i in 0:23){
        cols = c(cols, rep(i, 16))
      }

      name = paste0(rep(LETTERS[1:16], 24), unlist(mapply(function(i){rep(i,16)},1:24, SIMPLIFY=F, USE.NAMES=F)))

      values = eval(parse(text=paste("InCell()$well$`", input$featureCol,"`", sep="")))

      data = data.frame(cols,rows,values,name)

      # Highcharts options
      myChart=list(
        credits=list(
          enabled=FALSE
        ),

        chart=list(
          type='heatmap'
        ),

        title=list(
          text=input$featureCol,
          align='left'
        ),

        subtitle=list(
          text="Well-level averages",
          align='left'
        ),

        xAxis=list(
          categories=as.character(1:24),
          opposite=TRUE,
          lineWidth=0,
          minorGridLineWidth=0,
          minorTickLength=0,
          tickLength=0,
          lineColor='transparent'
        ),

        yAxis=list(
          reversed=TRUE,
          categories=LETTERS[1:16],
          lineWidth=0,
          minorGridLineWidth=0,
          minorTickLength=0,
          tickLength=0,
          lineColor='transparent',
          title=list(
            text=""
          )
        ),

        colorAxis=list(
          min=min(data$values),
          minColor='#ffffff',
          maxColor=getHighchartsColors()[1]
        ),

        legend=list(
          align='right',
          layout='vertical',
          margin=0,
          symbolHeight=320
        ),

        tooltip=list(
          headerFormat='{series.name} <br/>',
          pointFormat='{point.name}: <b>{point.value}</b><br/>'
        ),

        series=list(
          list(
            name=input$featureCol,
            borderWidth=1,
            data=JSONify(data, element.names=c("x", "y", "value", "name"))
          )
        )
      )

      return( list(chart=myChart) )
    }
  } )

  # Mini Histogram
  output$miniHist = renderHighcharts({
    if( !is.null(InCell()$well) ){
      v = eval(parse(text=paste0('InCell()$well$`',input$featureCol, '`')))
      h=hist(v, plot=F)

      # Highcharts options
      myChart=list(
        credits=list(
          enabled=FALSE
        ),
        chart=list(
          type='column'
        ),

        title=list(
          text=""
        ),

        subtitle=list(
          text="Histogram of Well Mean",
          align='left'
        ),

        xAxis=list(
          categories=as.character(h$breaks),
          labels=list(
            rotation=270
            )
        ),

        yAxis=list(
          title=list(
            text="Frequency"
          )
        ),

        legend=list(
          enabled=FALSE
        ),

        plotOptions=list(
          series=list(
            groupPadding=0.1,
            pointPadding=0
          )
        ),

        tooltip=list(
          headerFormat="{series.name}<br/>",
          pointFormat="<b>{point.y}</b><br/>"
        ),

        series=list(
          list(
            name=input$featureCol,
            data=h$counts
          )
        )
      )

      return( list(chart=myChart) )
    }
  })

  # Mini Stripchart
  output$miniFieldData = renderHighcharts({
    if( !is.null(input$highHeat) ){
      if( nchar(input$highHeat$x) == 1 ) clickX = paste0("0", input$highHeat$x+1) else clickX = input$highHeat$x+1
      well = paste0(LETTERS[input$highHeat$y+1], clickX)
      if(input$featureCol != "Cell Count"){
        d = eval(parse(text=paste("InCell()$cell$`", input$featureCol,"`", sep="")))
        fieldNames = unique( InCell()$field$Field[grep(well, InCell()$field$Well )] )

        l = mapply(function(field, n){
          wellIdx = grep(well, InCell()$cell$Well)
          fldIdx = grep(field, InCell()$cell$Field)
          index = fldIdx[which(fldIdx %in% wellIdx)]

          values = eval(parse(text=paste("InCell()$cell$`", input$featureCol,"`[index]", sep="")))

          mapply(function(v, n){
            c( (jitter(0, factor=5)+n), v )
          }, values, n-1, SIMPLIFY=F, USE.NAMES=F)


        }, fieldNames, 1:length(fieldNames), SIMPLIFY=F, USE.NAMES=F)


        # Highcharts Options
        myChart=list(
          credits=list(
            enabled=FALSE
          ),
          chart=list(
            type='boxplot'
          ),
          title=list(
            text=""
          ),
          subtitle=list(
            text=paste0("Individual Cell Values (", well, ")"),
            align='left'
          ),
          xAxis=list(
            categories=as.list(paste0("Fld ", as.character(fieldNames))),
            min=-0.2,
            max=0.2+length(fieldNames)-1
          ),
          yAxis=list(
            title=list(
              text="Value"
            )
          ),
          tooltip=list(
            headerFormat="{series.name}<br/>",
            pointFormat="<b>{point.y}</b><br/>"
          ),
          legend=list(
            enabled=FALSE
          ),

          series=list(
            list(
              name=input$featureCol,
              color=getHighchartsColors()[1],
              type="scatter",
              data=unlist(l, recursive=F),
              marker=list(
                fillColor = '#ffffff',
                lineWidth = 1,
                lineColor = getHighchartsColors()[1]
              )
            )
          )
        )


        return(list(chart=myChart))
      }
    }
  })




  ## Any|Any Tab ##
  # Large Plot
  output$AnyAny = renderHighcharts({
    noPlot = is.null(input$yFeat) || is.null(input$xFeat) || is.null(input$yTrans) || is.null(input$xTrans) ||
             is.null(input$yThresh) || is.null(input$xThresh) || is.null(InCell()) || is.null(REMP()) ||
             (input$xFeat == "Cell Count") || (input$yFeat == "Cell Count")

    if( !noPlot ){
      compounds = unique(REMP()$SAMPLE[which(REMP()$SAMPLE != "DMSO")])
      seriesData = mapply(function(cmpd, n){
        # Pull Values for the x and y data
        xF = InCell()$cellList[[input$xFeat]]
        yF = InCell()$cellList[[input$yFeat]]

        # Get the list of the wells that the selected compound is in
        conc = as.numeric(REMP()$CONC[which(REMP()$SAMPLE == cmpd)])
        concOrder = sort.int(conc, index.return=T)
        wells = as.character(REMP()$well[which(REMP()$SAMPLE == cmpd)])[concOrder$ix]

        # Generate X Data
        x = unlist(mapply(function(w){
          # Get data
          t = xF[[w]]
          c = if(input$xThresh != "None") input$xSlide else 1

          # Transform if needed
          if( input$xTrans == "Log10" ){
            t = log10(t)
            c = log10(c)
          }
          if( input$xTrans == "Log2" ){
            t = log2(t)
            c = log2(c)
          }

          # Threshold
          switch( input$xThresh,
                  "% Above" = prettyNum(sum(t > c)/length(t)*100, digits=3, width=4),
                  "% Below" = prettyNum(sum(t < c)/length(t)*100, digits=3, width=4),
                  "None" = prettyNum(mean(t), digits=3, width=4) )

        }, wells, SIMPLIFY=F, USE.NAMES=F))

        # Generate Y Data
        y = unlist(mapply(function(w){
          # Get data
          t = yF[[w]]
          c = if(input$yThresh != "None") input$ySlide else 1

          # Transform if needed
          if( input$yTrans == "Log10" ){
            t = log10(t)
            c = log10(c)
          }
          if( input$yTrans == "Log2" ){
            t = log2(t)
            c = log2(c)
          }

          # Threshold
          switch( input$yThresh,
                  "% Above" = prettyNum(sum(t > c)/length(t)*100, digits=3, width=4),
                  "% Below" = prettyNum(sum(t < c)/length(t)*100, digits=3, width=4),
                  "None" = prettyNum(mean(t), digits=3, width=4) )

        }, wells, SIMPLIFY=F, USE.NAMES=F))

        list(
          events=list(
            mouseOver = dimOtherSeries,
            mouseOut = resetAllSeries
          ),
          animation=FALSE,
          name=cmpd,
          type="scatter",
          data=JSONify(data.frame(x=as.numeric(x), y=as.numeric(y),
                                  name=paste0(prettyNum(concOrder$x, digits=3, width=4), "nM"),
                                  radius=log10(concOrder$x)*3),
                       element.names=c("x", "y", "name", "radius"))
        )

      }, compounds, length(compounds), SIMPLIFY=F, USE.NAMES=F)


      ## Highcarts Options ##
      myChart=list(
        credits=list(
          enabled=FALSE
        ),
        title=list(
          text=paste0(input$yFeat, " Versus ", input$xFeat),
          align='left'
        ),
        subtitle=list(
          text=input$cmpdAnyAny,
          align='left'
        ),
        xAxis=list(
          title=list(
            text= switch( input$xThresh,
                         "% Above" = paste0("% Cells with ", input$xFeat, " > ", input$xSlide),
                         "% Below" = paste0("% Cells with ", input$xFeat, " < ", input$xSlide),
                         "None" = input$xFeat)
          ),
          labels=list(
            format=if(input$xThresh == "None") '{value}' else '{value}%'
          ),
          min=if(input$xThresh == "None") NULL else 0,
          max=if(input$xThresh == "None") NULL else 100
        ),
        yAxis=list(
          title=list(
            text=switch( input$yThresh,
                         "% Above" = paste0("% Cells with ", input$yFeat, " > ", input$ySlide),
                         "% Below" = paste0("% Cells with ", input$yFeat, " < ", input$ySlide),
                         "None" = input$yFeat)
          ),
          labels=list(
            format=if(input$yThresh == "None") '{value}' else '{value}%'
          ),
          min=if(input$yThresh == "None") NULL else 0,
          max=if(input$yThresh == "None") NULL else 100
        ),
        legend=list(
          enabled = TRUE
        ),
        tooltip=list(
          enabled=FALSE,
          formatter="return this.series.name + ' [' + this.point.name + ']<br/>X: <b>' + this.x + '</b>   Y: <b>' + this.y + '</b>'"
        ),
        series=seriesData
      )

      return(list(chart=myChart))



    }
  })




  ## Feature Tab ##
  # Large Kernel Density Plot
  output$distribution = renderHighcharts({
    noPlot = is.null(negCtrlDist()) || (is.null(input$featureColDist) || !(input$featureColDist %in% names(InCell()$cellList)) ||
             is.null(input$logFeatureValues))

    if( !noPlot ){


      # Get the list of negative control wells
      negCtrlWells = REMP()$well[which(REMP()$SAMPLE == "DMSO")]

      # Get the list of the wells that the selected compound is in
      conc = as.numeric(REMP()$CONC[which(REMP()$SAMPLE == input$cmpdSelect)])
      concOrder = sort.int(conc, index.return=T)
      wells = as.character(REMP()$well[which(REMP()$SAMPLE == input$cmpdSelect)])[concOrder$ix]

      # Get the negative control well values
      if( input$logFeatureValues ){
        y = log10(unlist(
          mapply(function(negWell){
            InCell()$cellList[[ input$featureColDist ]][[ negWell ]]
          }, as.character(negCtrlWells), SIMPLIFY=T, USE.NAMES=F)))
      }else{
        y = unlist(
          mapply(function(negWell){
            InCell()$cellList[[ input$featureColDist ]][[ negWell ]]
          }, as.character(negCtrlWells), SIMPLIFY=T, USE.NAMES=F))
      }

      data=list()
      data$neg = JSONify(data.frame(x=density(y, na.rm=T)$x, y=density(y, na.rm=T)$y))

      # Generate the plots
      n = length(wells)
      for( i in 1:length(wells) ){
        # Get the values for the stat tests
        if( input$logFeatureValues ){
          x = log10(InCell()$cellList[[ input$featureColDist ]][[ wells[i] ]])
        }else{
          x = InCell()$cellList[[ input$featureColDist ]][[ wells[i] ]]
        }

        if( length(x) < 1 ){
          data[[paste0("c",i)]] = NULL
        }else{
          data[[paste0("c",i)]] = JSONify(data.frame(x=density(x, na.rm=T)$x, y=density(x, na.rm=T)$y))
        }
      }

      # Get the values for
      if( input$logFeatureValues ){
        x = log10(unlist(InCell()$cellList[[ input$featureColDist ]][ wells ]))
        label = "Log10 Feature Value"
      }else{
        x = unlist(InCell()$cellList[[ input$featureColDist ]][ wells ])
        label = "Feature Value"
      }
      data$all = JSONify(data.frame(x=density(x, na.rm=T)$x, y=density(x, na.rm=T)$y))

      seriesData=mapply(function(n){
        list(
          name=paste0(prettyNum(concOrder$x[n], digits=3, width=4)," nM"),
          type="line",
          linewidth=1,
          zIndex=1,
          data=data[[paste0("c",n)]],
          visible=is.even(n)
        )
      }, 1:length(wells), SIMPLIFY=F, USE.NAMES=F)

      seriesData[[length(wells)+1]] = list(
        name="DMSO",
        type="area",
        linewidth=1,
        zIndex=0,
        data=data$neg

        )



      ## Highcarts Options ##
      myChart=list(
        credits=list(
          enabled=FALSE
        ),
        title=list(
          text=input$cmpdSelect,
          align='left'
        ),
        subtitle=list(
          text="Kernel Density Distribution Estimate",
          align='left'
        ),
        xAxis=list(
          title=list(
            text=label
          )
        ),
        yAxis=list(
          title=list(
            text="Density"
          )
        ),
        tooltip=list(
          shared=TRUE,
          crosshairs=TRUE
        ),
        legend=list(
          enabled=TRUE
        ),
        series=seriesData
      )

      return(list(chart=myChart))
    }
  })

  # Mini Stat Z Plot
  output$miniZStat = renderHighcharts({
    if( !is.null(doseCurveStats()) && !is.null(input$featureColDist) &&
          (input$featureColDist %in% names(InCell()$cellList)) ){

      # get the z stat vectors
      ks.z = round(as.numeric(unlist(doseCurveStats()['ks.z', ])), 3)
      cvmts.z = round(as.numeric(unlist(doseCurveStats()['cvmts.z', ])), 3)
      conc = round(as.numeric(log10(unlist(doseCurveStats()['conc', ]))), 2)
      cells = round(as.numeric(unlist(doseCurveStats()['cellNum', ])), 3)




      ## Highcarts Options ##
      myChart=list(
        credits=list(
          enabled=FALSE
        ),
        title=list(
          text="",
          align='left'
        ),
        subtitle=list(
          text="Z-Score Normalized Test Statistics",
          align='left'
        ),
        xAxis=list(
          title=list(
            text="Log10 Concentration [nM]"
          )
        ),
        yAxis=list(
          title=list(
            text="Z-Score"
          ),
          tickPositions=sort(c(min(c(ks.z, cvmts.z, cells)), max(c(ks.z, cvmts.z, cells)), 0))
        ),
        tooltip=list(
          shared=TRUE,
          crosshairs=TRUE
        ),
        legend=list(
          enabled=TRUE
        ),
        series=list(
          list(
            name="K-S Test",
            type="line",
            linewidth=1,
            zIndex=1,
            data=JSONify(data.frame(x=conc, y=ks.z))
          ),
          list(
            name="CVMTS Test",
            type="line",
            linewidth=1,
            zIndex=1,
            data=JSONify(data.frame(x=conc, y=cvmts.z))
          ),
          list(
            name="Cell Number",
            type="area",
            linewidth=1,
            zIndex=0,
            data=JSONify(data.frame(x=conc, y=cells))
          )
        )
      )

      return(list(chart=myChart))
    }

  })


  # Mini Area Under the Curve Plot
  output$miniAUC = renderHighcharts({
    if( !is.null(doseCurveStats()) && !is.null(input$featureColDist) &&
          (input$featureColDist %in% names(InCell()$cellList)) ){

      # get the z stat vectors
      ks.z = round(as.numeric(unlist(doseCurveStats()['ks.z', ])), 3)
      cvmts.z = round(as.numeric(unlist(doseCurveStats()['cvmts.z', ])), 3)
      conc = round(as.numeric(log10(unlist(doseCurveStats()['conc', ]))), 2)
      cells = round(as.numeric(unlist(doseCurveStats()['cellNum', ])), 3)

      # Get indexes of the ordered conc
      id = order(conc)

      # Calculate the Area Under the Curves
      ks.AUC = auc(x=conc[id], y=ks.z[id])
      cvmts.AUC = auc(x=conc[id], y=cvmts.z[id])

      ## Highcarts Options ##
      myChart=list(
        credits=list(
          enabled=FALSE
        ),
        title=list(
          text="",
          align='left'
        ),
        subtitle=list(
          text="Area Under The Curve",
          align='left'
        ),
        xAxis=list(
          title=list(
            text=""
          ),
          categories=c("K-S Test", "CVMTS Test")
        ),
        yAxis=list(
          title=list(
            text="AUC"
          )
        ),
        tooltip=list(
          shared=TRUE,
          crosshairs=TRUE
        ),
        legend=list(
          enabled=FALSE
        ),
        series=list(
          list(
            name="AUC",
            type="column",
            data=round(c(ks.AUC, cvmts.AUC), 3)
          )
        )
      )

      return(list(chart=myChart))
    }

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

  # Update File dropdowns
  observe({
    input$updateDropbox

    updateSelectizeInput(session, 'InCellDB', choices = c("Choose from Dropbox", names(Dropbox()$csv)))
    updateSelectizeInput(session, 'annotationDB', choices = c("Choose from Dropbox", names(Dropbox()$txt)))
  })

  # Sync the Feature dropdowns
  observe({
    updateSelectizeInput(session, 'featureColDist', selected=input$featureCol)
  })
  observe({
    updateSelectizeInput(session, 'featureCol', selected=input$featureColDist)
  })



})
