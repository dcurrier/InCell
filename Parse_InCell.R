#
# Helper scripts for parsing InCell data files
#
# Duane Currier - 08/01/2014
#

ReadInCell = function(file, tables=c("cell", "field", "well"), progressBar=FALSE){

  # Begin Progress Bar
  if( progressBar ){ eval(setProgress(value = 0.01, detail='Scanning File'), parent.frame(n=2)) }

  # Scan the entire file
  temp = scan(file, what="Character", sep="\n", blank.lines.skip=F)

  # Make the output list
  out = vector("list", length(tables))
  names(out) = tables

  # Find the Cell and Field Level data section
  if( any(c("cell", "field") %in% tables) ){
    well.field = grep("[[:upper:]] - [[:digit:]]+\\(fld [[:digit:]]+\\)", temp)

    if( length(well.field) == 0 ) {
      well.field = grep("[[:upper:]] - [[:digit:]]", temp)
      fld.name = 1
    }

    # Find headers
    cellHead = grep("Well,", temp)[which(grep("Well,", temp) < min(well.field))]
    fieldHead = grep("Well,", temp)[which(grep("Well,", temp) %in% min(well.field):max(well.field))][1]
  }

  # Update Progress Bar
  if( progressBar ){ eval(setProgress(value = 0.1, detail='Parsing Cell Level Data'), parent.frame(n=2)) }


  # Get Cell Level Data
  if( "cell" %in% tables ){
    # Find data
    cell = well.field[which(well.field < fieldHead)]

    # Make DataFrame
    cellData = read.csv(file, header=F, skip=cellHead, nrows=length(cell) )
    head1 = strsplit(temp[cellHead-1], ",")[[1]]
    head2 = strsplit(temp[cellHead], ",")[[1]]
    header = unlist(mapply(function(t,b){
                            if(t==""){b}else{paste(t,b,sep="-")}
                          }, head1, head2, SIMPLIFY=F, USE.NAMES=F))
    names(cellData) = header
    cellData = cellData[, !is.na(names(cellData))]

    if( is.null(fld.name) ){
      # Parse Well names
      wellID = gsub("([[:upper:]]) - ([[:digit:]]{2})([[:punct:]]fld [[:digit:]]+[[:punct:]])", paste0("\\1","\\2"), cellData$Well)
      wellID = gsub("([[:upper:]]) - ([[:digit:]]{1})([[:punct:]]fld [[:digit:]]+[[:punct:]])", paste0("\\1","0\\2"), wellID)
      # Parse Field IDs
      fieldID = gsub("([[:upper:]] - [[:digit:]]+)[[:punct:]]fld ([[:digit:]]+)[[:punct:]]", paste0("\\2"), cellData$Well)
    }else{
      # Parse Well names
      wellID = gsub("([[:upper:]]) - ([[:digit:]]{2})", paste0("\\1","\\2"), cellData$Well)
      wellID = gsub("([[:upper:]]) - ([[:digit:]]{1})", paste0("\\1","0\\2"), wellID)

      fieldID = fld.name
    }

    # Append to table
    cellData = data.frame(Well=wellID, Field=fieldID, cellData[, which(names(cellData) != "Well")])
    names(cellData) = c("Well", "Field", header[-1])
    out[["cell"]] = cellData

    # Update Progress Bar
    if( progressBar ){ eval(setProgress(value = 0.2, detail='Parsing Cell Level Data'), parent.frame(n=2)) }


    # Split the Cell Level data into a list of lists
    dataCols = names(out$cell)[!(names(out$cell) %in% c("Well", "Field", "Cell"))]

    # Reformat Cell table to a list of features, each containing a list of the wells, each containing a vector of the values
    feature=mapply(function(col){
      w=mapply(function(well){
        # Get the values for the current column
        eval(parse(text=paste0("out$cell$`",col,"`[which(out$cell$Well == '",well,"')]")))
      }, as.character(unique(out$cell$Well)), SIMPLIFY=T, USE.NAMES=T)
      #names(w) = unique(out$cell$Well)

      # Update Progress Bar
      if( progressBar ){ eval(setProgress(value = (0.2+(which(dataCols == col)/length(dataCols)*0.6)),
                                          detail='Parsing Cell Level Data'), parent.frame(n=3)) }

      list(w)
    }, dataCols, SIMPLIFY=T, USE.NAMES=T)
    #names(feature) = dataCols

    out[["cellList"]] = feature
  }

  # Update Progress Bar
  if( progressBar ){ eval(setProgress(value = 0.8, detail='Parsing Field Level Data'), parent.frame(n=2)) }


  # Get Field Level Data
  if( "field" %in% tables ){
    # Find data
    field = well.field[which(well.field > fieldHead)]

    # Make DataFrame
    fieldData = read.csv(file, header=F, skip=fieldHead, nrows=length(field) )
    head1 = strsplit(temp[fieldHead-1], ",")[[1]]
    head2 = strsplit(temp[fieldHead], ",")[[1]]
    header = unlist(mapply(function(t,b){
                            if(t==""){b}else{paste(t,b,sep="-")}
                          }, head1, head2, SIMPLIFY=F, USE.NAMES=F))
    names(fieldData) = header
    fieldData = fieldData[, !is.na(names(fieldData))]

    if( is.null(fld.name) ){
      # Parse Well names
      wellID = gsub("([[:upper:]]) - ([[:digit:]]{2})([[:punct:]]fld [[:digit:]]+[[:punct:]])", paste0("\\1","\\2"), fieldData$Well)
      wellID = gsub("([[:upper:]]) - ([[:digit:]]{1})([[:punct:]]fld [[:digit:]]+[[:punct:]])", paste0("\\1","0\\2"), wellID)
      # Parse Field IDs
      fieldID = gsub("([[:upper:]] - [[:digit:]]+)[[:punct:]]fld ([[:digit:]]+)[[:punct:]]", paste0("\\2"), fieldData$Well)
    }else{
      # Parse Well names
      wellID = gsub("([[:upper:]]) - ([[:digit:]]{2})", paste0("\\1","\\2"), fieldData$Well)
      wellID = gsub("([[:upper:]]) - ([[:digit:]]{1})", paste0("\\1","0\\2"), wellID)

      fieldID = fld.name
    }


    # Append to table
    fieldData = data.frame(Well=wellID, Field=fieldID, fieldData[, which(names(fieldData) != "Well")])
    names(fieldData) = c("Well", "Field", header[-1])
    out[["field"]] = fieldData
  }

  # Update Progress Bar
  if( progressBar ){ eval(setProgress(value = 0.9, detail='Parsing Well Level Data'), parent.frame(n=2)) }


  # Get Well Level Data
  if( "well" %in% tables ){
    # Find header
    wellHead = grep("Well,,", temp)

    # Find data
    well = grep("[[:upper:]] - [[:digit:]]+,", temp)

    # Make DataFrame
    wellData = read.csv(file, header=F, skip=wellHead, nrows=length(well) )
    head1 = strsplit(temp[wellHead-1], ",")[[1]]
    head2 = strsplit(temp[wellHead], ",")[[1]]
    header = unlist(mapply(function(t,b){
                            if(t==""){b}else{paste(t,b,sep="-")}
                          }, head1, head2, SIMPLIFY=F, USE.NAMES=F))
    names(wellData) = header
    wellData = wellData[, !is.na(names(wellData))]
    wellData$Well = gsub("([[:upper:]]) - ([[:digit:]]{2})", paste0("\\1","\\2"), wellData$Well)
    wellData$Well = gsub("([[:upper:]]) - ([[:digit:]]{1})", paste0("\\1","0\\2"), wellData$Well)

    out[["well"]] = wellData
  }


  # Update Progress Bar
  if( progressBar ){ eval(setProgress(value = 0.99, detail='Wrapping up'), parent.frame(n=2)) }

  return(out)

}