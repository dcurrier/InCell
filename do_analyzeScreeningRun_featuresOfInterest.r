# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#	OBJECTIVE: Given Z-scores from a High-Content Screening Run,
#			  Fit Dose-Response Curves

# NOTE: we didn't anayze uniqPlateIDs[1]

#	Assumes input Files have "Scan 1" in their file name
#			ExptlDesign.txt has two columns: AssayPlate [Col1] and PlateID [Col2]
#			PlateMaps contains the platemaps for all screened plates
#	BSUB:
#	bsub -P ChemBio -q priority -J "HCA2[1-12]" -app R-2.15.0 R --no-save --file=do_plateAnalysis.r
#	ChangeLog
#		[February 19, 2013] 
# 		- originated [modified from Hdo_HCA_CVMTS_forPlateI]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 							DEFINE CONSTANTS	
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 					Directories and File Names
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


scriptDir		<- "U:/Transfer_TMP/BirgitJon/July7_2014/SourceScripts/"
inputDir		<- "U:/Transfer_TMP/BirgitJon/July7_2014/Input/"
resultsDir		<- "U:/Transfer_TMP/BirgitJon/July7_2014/Results/"

plateFile       <- "PlateSummary.txt"

funFile			          <- "screeningRunAnalysis_functions.r"
plateResultsSuffix        <-   "_CS_NoScale_Zscores_featsInterest"
plateResultsPrefix        <- ""

#screeningRunResultsSuffix <- "_CS_NoScale_CurveScores"
screeningRunResultsSuffix <- "_CS_NoScale_CurveScores_featsInterest"
exptlDesignFile <- "ExptlDesign.txt"
plateMapsFile   <- "PlateSummary.txt"
cmpdAliasFile   <- "SJReg_Aliases.txt"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 					Quality-control and Statistical Methods Constants
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

distStatMethod	<- "CVMTS" # "KS" or "CVMTS"
centerRobust	<- TRUE # TRUE = use median/mad; FALSE = meanSD
rSquaredMin     <- 0.8
numObjMin       <- 25

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 					Screen-specific Constants for this Screen
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# stain1		<- "Hoescht"
# stain1Names <- c("Major Axis Length", "Gaussian", "Peak Intensity", "Total Intensity")

# stain2		<- "Alpha-Tubulin"
# stain2Names <- c("Major Axis Length", "Gaussian", "Mean Intensity",  "Total Intensity", "Peak Intensity")
# stain2CntPercent <- TRUE # Take Count as a Percent of Stain 1?

# stain3		<- "H2AX"
# stain3Names <- NULL # Measured these, but Are Disregarded:
#						c("Major Axis Length", "Gaussian", "Total Intensity", "Peak Intensity")
# stain3CntPercent <- TRUE # Take Count as a Percent of Stain 1?

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	
# # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	   
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 							SCRIPT BEGINS HERE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

source(paste(scriptDir, funFile, sep = ""))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# (1) 	Read Experimental Design
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
exptlDesign  <- read.delim(file = paste(inputDir, exptlDesignFile, sep = ""), stringsAsFactors = FALSE)
names(exptlDesign)[names(exptlDesign) == "Drug.plate"] <- "PlateID"
names(exptlDesign)[names(exptlDesign) == "Control.plate"] <- "CtrlPlateID"
names(exptlDesign)[names(exptlDesign) == "Assay.plate"] <- "AssayPlate"
exptlDesign$FileName <- sub(".CSV","",exptlDesign$FileName)

                               



# 0.4.28.2014
# uniqPlateIDs <- unique(exptlDesign[, 2])
uniqPlateIDs <- unique(exptlDesign$PlateID)

plateMaps    <- read.delim(file = paste(inputDir,plateMapsFile, sep=""), stringsAsFactors = FALSE,
						   header = TRUE, row.names = NULL)
# Make Compound Aliases
# write.table(plateMaps[!duplicated(plateMaps$SAMPLE), c("SAMPLE", "synonym")],
#            file = paste(inputDir, "TMP_SJReg_Aliases.txt", sep = ""), row.names = FALSE,
#            quote = FALSE, sep = "\t")
				   
						   
cmpdAliases  <- read.delim(file = paste(inputDir, cmpdAliasFile, sep = ""))

scoredDR       <- NULL
cmpds_DR       <- NULL
cmpdAliases_DR <- NULL


# Loop over each Plate ID
for (plateID in uniqPlateIDs[-1]) {
	# Select Biological Replicates of this Plate ID
	# assayPlates <- exptlDesign[exptlDesign[, 2] %in% plateID, 1] 
	# zFileNames  <- paste(resultsDir, plateResultsPrefix, assayPlates, 
	#					 plateResultsSuffix, ".Rdata", sep ="")
	assayPlates <- exptlDesign$AssayPlate[exptlDesign$PlateID %in% plateID]
    ctrlPlates  <- exptlDesign$CtrlPlateID[exptlDesign$PlateID %in% plateID]
	
	zFileNames  <- paste(resultsDir, exptlDesign$FileName[exptlDesign$PlateID %in% plateID],
						 plateResultsSuffix, ".Rdata", sep ="")
	# [05.1.2014] Include CtrlPlateID 
	# plateMap    <- plateMaps[which(plateMaps$currplateid %in% plateID), ]
	plateMap    <- plateMaps[which(plateMaps$currplateid %in% c(plateID, ctrlPlates)), ]
	
	# [05.1.2014] Go by SJ number (not batch Number
	# plateAlias  <- cmpdAliases[which(cmpdAliases$regno %in% unique(plateMap$SAMPLE)), ]
	plateAlias  <- cmpdAliases[which(cmpdAliases$regno %in% unique(plateMap$SAMPLE)), ]
	
	# Initialize Plot
	pdfFile <- paste(plateID, screeningRunResultsSuffix, ".pdf", sep = "")
	pdf(file = paste(resultsDir, pdfFile, sep = ""), pointsize =8)
	# par(mfrow = c(4,2))
	# par(mfrow = c(5, 3), pty = "s")
	par(pty = "s")
	
	# Fit Curves and Make Plots
	plateIDResults <- f_procPlateMapSimple(zFileNames, plateMap, plateAlias)
	scoredDR       <- rbind(scoredDR, plateIDResults$scoredDR)
	cmpds_DR       <- c(cmpds_DR, plateIDResults$cmpds)
	cmpdAliases_DR <- c(cmpdAliases_DR, plateIDResults$cmpdAliases)
	featNames      <- plateIDResults$featNames
	dev.off()
	}

# Save Screening Run Results
resultsFile <- paste("ScreeningResults", screeningRunResultsSuffix, ".Rdata", sep = "")
save.image(file = paste(resultsDir, resultsFile, sep = ""), compress = TRUE)
quit(save="no") 