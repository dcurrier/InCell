library(RColorBrewer)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
add.alpha <- function(col, alpha=1){
	# From
	# http://lamages.blogspot.com/2013/04/how-to-change-alpha-value-of-colours-in.html 
 if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
	}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	
f_plotDistsSimple <- function (zFileNames, plateMap, plateAlias) {
	# Consider only One Stain and One Feature of Interest
	# -hard coded feature of Interest
	# Inputs : zFileNames -  a vector of filenames corresponding to
	#							 zscores from replicated plates
	#		  plateMap - the plate map containing the annotation
	# 
	# (1) Collate Zscores from Replicate Plates and
	#	 Set Distribution Features with Insufficient Objct Number to NA

	for (fileName in zFileNames) {
		load(fileName)
		
	
	# (2) Make Cmpd Well Names Consistent
	cmpdWells  <- wellNames[-idxNegCtrl]
	idxFixName <- which(nchar(cmpdWells) < 3)
	if (length(idxFixName) != 0) {
		splitNames <- sapply(cmpdWells[idxFixName], function(x) strsplit(x,"")[[1]])
		fixedNames <- paste(splitNames[1,], splitNames[2,], sep="0")
		cmpdWells[idxFixName] <- fixedNames
		}
	
	# (3) Loop Over Each Compound
    cmpds       <- unique(as.vector(plateMap$SAMPLE))
	cmpdAliases <- cmpds
	# par(mfcol = c(4, 3))
	
	for (cmpd in cmpds) {
		cmpdNum   <- which(cmpds == cmpd)
		#cmpdMap   <- plateMap[as.vector(plateMap$SAMPLE) %in% cmpd, ]
		cmpdMap   <- plateMap[as.vector(plateMap$SAMPLE) %in% cmpd & 
		                      as.vector(plateMap$well) %in% cmpdWells, ]
		
		idxWells  <- which(cmpdWells %in% as.vector(cmpdMap$well))
		cmpdConc  <- cmpdMap$CONC / 1000 
		# 05.01.2014 
		# Be strict in enforcing the order of compound concentrations
		idxWells <- idxWells[order(cmpdConc, decreasing = TRUE)]
		cmpdConc <- sort(cmpdConc, decreasing = TRUE)
		# 
		cmpdAlias <- as.vector(plateAlias$alias[which(plateAlias$regno %in% cmpd)])
		if (length(cmpdAlias) != 0) cmpdAliases[cmpdNum] <- cmpdAlias
		
		# 05.01.2014
		# Make One PDf for each compound
		plotPrefix <- strsplit(zFileNames, "/")
		plotPrefix <- strsplit(plotPrefix[[1]][length(plotPrefix[[1]])], 
		                       plateResultsSuffix)[[1]][1]
		if (length(idxWells) <= 10) {
			pdf(file = paste(plotDir, plotPrefix, "_",cmpd, "_",cmpdAlias, ".pdf", sep = ""), 
				 pointsize = 8) #, width = 200, height = 200)	
			} else {
			pdf(file = paste(plotDir, plotPrefix, "_",cmpd, "_",cmpdAlias, ".pdf", sep = ""), 
				 pointsize = 8, width = 21, height = 21)	
			}
		# 05.19.2014 [ changed number of columns from 5 
		# par(mfcol = c(length(idxWells), 5), #length(inputHeader) - 1), 
		par(mfcol = c(length(idxWells), length(inputHeader) - 1), 
  	   	    mex =  0.5)# 0.625) # 0.5 is too small | 0.75 is slightly too big	
		# Loop Over Features
		for (j in seq(length(inputHeader) - 1)) {
			for (i in rev(idxWells)){ # order of Increasing Concentration
		    
			idxWell      <- seq(wellIdxStain1[seq(nrow(wellIdxStain1))[-idxNegCtrl][i], 1],
								wellIdxStain1[seq(nrow(wellIdxStain1))[-idxNegCtrl][i], 2])
			idxOtherWells <- f_combIdx(wellIdxStain1, idxNegCtrl)
			
			# Stop Plot if there are < 2 wells
			if (length(idxWell) < 2) {
			   plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
			        main = paste(round(rev(cmpdConc)[which(rev(idxWells == i))], 3), "uM\n",
						             inputHeader[-1][j]))
			   } else{
				cmpdResp <- log10(na.omit(inputData[idxWell, j + 1]))
				dmsoResp <- log10(na.omit(inputData[idxOtherWells, j + 1]))
				# Get Breaks
				myBreaks <- hist(c(cmpdResp, dmsoResp), plot = FALSE, breaks = 100)$breaks
				# Plot DMSO Resp: Hist
				hist(dmsoResp, freq = FALSE, breaks = myBreaks, 
				      
				# 05.19.2014 Changed Title to Simply have Conc \n Feature
				#            and xlabel to have stats
				# lwd = 0.5, xlab = "", ylab = "", 
				#	  main = paste(round(rev(cmpdConc)[which(rev(idxWells == i))], 3), "uM",
				#		             inputHeader[-1][j], "\n",
  	  			#		             "KS =", round(f_ks(cmpdResp, dmsoResp), 2),
				#					 "CVMTS =", round(f_cvmts(cmpdResp, dmsoResp), 2),
				#					 "n =", length(cmpdResp)), cex.main = 1,
					  lwd = 0.5, ylab = "", xlab = "",
				      main = paste(round(rev(cmpdConc)[which(rev(idxWells == i))], 3), "uM\n",
						             inputHeader[-1][j]),
					 cex.main = 1,
					 col = brewer.pal(4,"Paired")[1], border = brewer.pal(4,"Paired")[1],
					 axes = FALSE, ylim = c(0, max(c(pretty(density(dmsoResp)$y), 
					                                 pretty(density(cmpdResp)$y)))))
				mtext(paste("KS =", round(f_ks(cmpdResp, dmsoResp), 2),			 
  	  			 				   "CVMTS =", round(f_cvmts(cmpdResp, dmsoResp), 2),
									"n =", length(cmpdResp)), side = 1, cex = 0.5) 
									 
				# Add Cmpd Resp
				hist(cmpdResp, freq = FALSE, breaks = myBreaks, add = TRUE,
 				     col = add.alpha(brewer.pal(6,"Paired")[5], 0.5), 
					 border = add.alpha(brewer.pal(6,"Paired")[5], 0.5))
				# Add Densities
				lines(density(dmsoResp), lwd = 2, col = brewer.pal(4,"Paired")[2])	 
				lines(density(cmpdResp), lwd = 2, col = brewer.pal(6,"Paired")[6])
				}
			}
			}
		dev.off()	
		}
	}
}


scriptDir		<- "U:/Transfer_TMP/BirgitJon/July7_2014/SourceScripts/"
inputDir		<- "U:/Transfer_TMP/BirgitJon/July7_2014/Input/"
resultsDir		<- "U:/Transfer_TMP/BirgitJon/July7_2014/Results/"
plotDir         <- paste(resultsDir, "Distributions/", sep = "")

plateResultsSuffix   <- "_CS_NoScale_Zscores_featsInterest"
plateResultsPrefix        <- ""

cmpdAliasFile   <- "SJReg_Aliases.txt"
exptlDesignFile <- "ExptlDesign.txt"
plateMapsFile   <- "PlateSummary.txt"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# (1) 	Read Experimental Design and Drop a File
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
exptlDesign  <- read.delim(file = paste(inputDir, exptlDesignFile, sep = ""), stringsAsFactors = FALSE)
names(exptlDesign)[names(exptlDesign) == "Drug.plate"] <- "PlateID"
names(exptlDesign)[names(exptlDesign) == "Control.plate"] <- "CtrlPlateID"
names(exptlDesign)[names(exptlDesign) == "Assay.plate"] <- "AssayPlate"
exptlDesign$FileName <- sub(".CSV","",exptlDesign$FileName)
# Drop CPC DNA H2AX488 20140619.CSV
exptlDesign <- exptlDesign[exptlDesign$FileName != "CPC DNA H2AX488 20140619", ]

uniqPlateIDs <- unique(exptlDesign$PlateID)

plateMaps    <- read.delim(file = paste(inputDir,plateMapsFile, sep=""),
						   header = TRUE, row.names = NULL)
cmpdAliases  <- read.delim(file = paste(inputDir, cmpdAliasFile, sep = ""))

# Loop over each Plate ID
for (plateID in uniqPlateIDs) {
	# Select Biological Replicates of this Plate ID
	assayPlates <- exptlDesign$AssayPlate[exptlDesign$PlateID %in% plateID]
	ctrlPlates  <- exptlDesign$CtrlPlateID[exptlDesign$PlateID %in% plateID]
	zFileNames  <- paste(resultsDir, exptlDesign$FileName[exptlDesign$PlateID %in% plateID],
						 plateResultsSuffix, ".Rdata", sep ="")
	#zFileNames  <- paste(resultsDir, assayPlates, " Scan 1", 
	#					 plateResultsSuffix, ".Rdata", sep ="")					 
	plateMap    <- plateMaps[which(plateMaps$currplateid %in% c(plateID, ctrlPlates)), ]
	plateAlias  <- cmpdAliases[which(cmpdAliases$regno %in% unique(plateMap$SAMPLE)), ]
	
	f_plotDistsSimple(zFileNames, plateMap, plateAlias)
	}
	
