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

    # Find headers
    cellHead = grep("Well,Cell", temp)
    fieldHead = grep("Well,", temp)[which(grep("Well,", temp) %in% min(well.field):max(well.field))]
  }

  # Update Progress Bar
  if( progressBar ){ eval(setProgress(value = 0.3, detail='Parsing Cell Level Data'), parent.frame(n=2)) }


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

    out[["cell"]] = cellData
  }

  # Update Progress Bar
  if( progressBar ){ eval(setProgress(value = 0.6, detail='Parsing Field Level Data'), parent.frame(n=2)) }


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

    out[["field"]] = fieldData
  }

  # Update Progress Bar
  if( progressBar ){ eval(setProgress(value = 0.8, detail='Parsing Well Level Data'), parent.frame(n=2)) }


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

    out[["well"]] = wellData
  }

  # Update Progress Bar
  if( progressBar ){ eval(setProgress(value = 0.99, detail='Wrapping up'), parent.frame(n=2)) }

  return(out)

}