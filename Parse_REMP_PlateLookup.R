#
# Helper scripts for parsing REMP Plate Lookup files
#
# Duane Currier - 08/06/2014
#

Parse_REMP = function(path, ctl, noCtl){
  # Read in the data file
  data = read.table(path, header=T, sep="\t", fill=TRUE)

  nPlates = length(unique(data$currplateid))

  # Check for too few plates in the file
  if( nPlates < 2 && !noCtl ){
    return(NULL)
  }

  # check for too many plates in the file
  if( nPlates > 2 ){
    n=unlist(mapply(function(id){
      sum(data$currplateid == id)
    },unique(data$currplateid), SIMPLIFY=F, USE.NAMES=F))

    ctlPlate = unique(data$currplateid)[ min(( 1:length(n) )[which( n < (length(ctl)*16) )]) ]
    cmpdPlate = unique(data$currplateid)[ min(( 1:length(n) )[which( n > (length(ctl)*16) )]) ]

    data = data[ data$currplateid %in% c(ctlPlate, cmpdPlate) ,]
  }

  # Generate a 384 well list
  wells = unlist(mapply(function(L){
    paste0(L, formatC(1:24, width=2, format='d', flag='0'))
  }, LETTERS[1:16], SIMPLIFY=F, USE.NAMES=F))

  # Find the unannotated wells in the control columns
  DMSO = wells[which(!(wells %in% data$well))][as.integer(gsub("[[:upper:]]", "", wells[which(!(wells %in% data$well))])) %in% ctl]

  # Create data table rows for the DMSO wells
  DMSOwells = t(as.data.frame(mapply(function(well){
                                c( data$currplateid[dim(data)[1]],
                                          well,
                                          "DMSO",
                                          rep("", dim(data)[2] - 3)
                                          )
                              }, DMSO, SIMPLIFY=F, USE.NAMES=F)))
  row.names(DMSOwells) = 1:dim(DMSOwells)[1]
  colnames(DMSOwells) = names(data)

  # Change any mM concentrations to uM
  mili = which(data$CU == "mM")
  data$CONC[mili] = data$CONC[mili] * 1000
  data$CU[mili] = "uM"

  # Bind the DMSO well rows to the main table
  data = rbind(data, DMSOwells)

  # Process the synonyms
  synonyms = unlist(mapply(function(value){
    if(is.na(value)) return(NA)

    first = strsplit(as.character(value), " \\| ")[[1]][1]
    if( nchar(first) > 20 ) first = strtrim(first, 20)

    return(first)

  }, data$synonym, SIMPLIFY=F, USE.NAMES=F))

  # Bind the synonyms to the end of the table
  data$altSynonym = synonyms
  data$comboId = paste0(data$SAMPLE, " (", synonyms, ")")

  return( data )
}

