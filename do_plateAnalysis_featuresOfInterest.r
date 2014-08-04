# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#	OBJECTIVE: Calculate Z-scores from a High-Content Screening Run
#              Using InCell with Only One Stain
#   Notes:
#         I ran this by simply hard-coding the plateNum == 1
#                  then cleared memory and repeated for plateNum == 2
# ChangeLog:
# [05.19.2014]
# - Restrict Analysis to Features of Interest
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


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
funFile			<- "plateAnalysis_functions.r"
# resultsSuffix   <- "_CS_NoScale_Zscores_FixedDMSO"
resultsSuffix   <- "_CS_NoScale_Zscores_featsInterest"
plateFile       <- "PlateSummary.txt"
exptFile		<- "ExptlDesign.txt"            

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 					Quality-control and Statistical Methods Constants
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

distStatMethod	<- "CVMTS" # "KS" or "CVMTS"
useRobust		<- TRUE # TRUE = use median/mad; FALSE = meanSD
scaleDistBoot   <- FALSE # Scale Dist Stats by bootsrap SD [Perlman] ?
scaleDistFeat   <- FALSE # Should Distribution z-scores be scaled ?

doNSFilter		<- FALSE # Filter non-specific secondary or Tertiary Stains?
makeQCPlot		<- FALSE # Generate a PDF containing the non-specific filter plots
exclEdgeNegCtrl <- FALSE # # Exclude the negative control wells on the edge?

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 					 				 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	   
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 							SCRIPT BEGINS HERE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# (0)							Read Exptl Design
#                             Because Not every Plate may have same NegCtrl plte
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

exptlDesign <- read.delim(paste(inputDir, exptFile, sep = ""), 
                          stringsAsFactors = FALSE)
for (plateNum in seq(nrow(exptlDesign))) {
     inputFile <- exptlDesign$FileName[plateNum]	
     drugPlateID <- exptlDesign$Drug.plate[plateNum]		 
	 ctrlPlateID <- exptlDesign$Control.plate[plateNum]		 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# (1) 	Read Plate Data 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#plateNum	<- 1# + 1

#inputFile	<- dir(inputDir, pattern = ".CSV")[plateNum] # Assumes Scan 1 in Name
tmpHeader	<- read.csv(file = paste(inputDir, inputFile, sep = ""), 
						 header = FALSE, row.names=NULL, stringsAsFactors = FALSE,
			 			 nrows = 2, skip = 18) 
						 
# Read Input Data (2 steps:
#                   (1) All lines to figure out the last cell-level line
#                   (2) Use this specification to re-read data
inputData	<- read.csv(file = paste(inputDir, inputFile, sep = ""),
						 header = FALSE, row.names = NULL, stringsAsFactors = FALSE,
						 skip = 20)
idxLastRow <- which(!grepl(" - ", inputData[,1]))[1] - 1
#idxLastRow <- which(inputData[,1] == "Well")[1] - 3
inputData	<- read.csv(file = paste(inputDir, inputFile, sep = ""),
						 header = FALSE, row.names = NULL, stringsAsFactors = FALSE,
						 skip = 20, nrows = idxLastRow)
						 
# Do Housekeeping on Feature Names
inputHeader	<- rep(NA, ncol(tmpHeader))					 
for (i in 1:ncol(tmpHeader)) {
     inputHeader[i] <- as.vector(tmpHeader[2, i])
     if (!is.na(tmpHeader[1, i]) & as.vector(tmpHeader[1, i]) != ""){
	     inputHeader[i] <- paste(as.vector(tmpHeader[1, i]), as.vector(tmpHeader[2, i]), sep =": ")
		 }
	}
rm(tmpHeader)	 

# Rename Row names to be be consistent with Acumen Format
tmpNames <- as.vector(inputData[, 1])
tmpNames <- sub(" - ", "", tmpNames)
tmpNames <- sub("(fld 1)","", tmpNames,fixed = TRUE)
tmpNames <- sub("(fld 2)","", tmpNames,fixed = TRUE)
inputData[, 1] <- tmpNames
gc()						 					

# Pad Well Names with Zeros to be consistent with PlateMaps
wellNames <- inputData[, 1] 											
idxFix    <- which(nchar(wellNames) == 2)
fixedNames  <- sapply(wellNames[idxFix], 
                      function(x) paste(strsplit(x,"")[[1]],collapse = "0")) 
wellNames[idxFix] <- as.vector(fixedNames)
inputData[, 1] <- wellNames
rm(wellNames, idxFix, fixedNames)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# (2) Perform Quality-Control
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

plateMaps <- read.delim(paste(inputDir, plateFile, sep = ""), stringsAsFactors = FALSE)
# Limit Plate Maps to drugPlateID and ctrlPlateID
plateMaps <- plateMaps[plateMaps$currplateid %in% c(drugPlateID, ctrlPlateID), ]

negCtrlWells <- unique(inputData[! inputData[, 1] %in% plateMaps$well, 1]) 


# (A) Exclude Negative Control Wells on Edges
if (exclEdgeNegCtrl) {
	# edgeWells <- c("A21", "A22", "A23", "A24", "P21", "P22", "P23", "P24",
					# "B24", "C24", "D24", "E24", "F24", "G24", "H24", "I24", 
					# "J24", "K24","L24", "M24", "N24", "O24")
	
	 edgeWells <- c(paste(LETTERS[1:4], 21, sep = ""), # Cycloheximide
	               paste(LETTERS[9:12],21, sep = ""), # Cycloheximide
				   paste(LETTERS[5:8], 22, sep = ""), # Staurosporine
	               paste(LETTERS[13:16],22, sep = ""), # Staurosporin
				   paste(LETTERS[1:16], 23, sep = ""), # Cycloheximide
				   paste(LETTERS[1:16], 24, sep = "")) # Staurosporine
				   
	idxExcl   <- which(inputData[, 1] %in% edgeWells)
	inputData <- inputData[-idxExcl, ]
	}


# (C) Remove Unnecessary Features [any feature NOT in the following vector]
# [04.28.2014]
# inputData <- inputData[, -11] # Last Column contains no infomation
# inputData <- inputData[, -c(2, 9:10)] # Col 2 = Cell Number | Cols 9 and 10 are binary
# A Feature is Unnecesary if there are < 100 unique values across all measurmenents
# idxFeatRemove <- apply(inputData, 2,function(x) length(unique(x))) < 100  
# [05.19.2014]
featInterest <- c("Nuclei: Nuc/Cell Intensity", # : Ratio of intensities sampled in the nuclear and cytoplasm regions
                 "Nuclei: Nuc Area", # Area of identified nucleus
				 "Nuclei: Nuc Elongation", # Mean ratio of the short axis of the nucleus to the long axis of the nucleus. If value is 1 object is center-symmetric
				  "Nuclei: Nuc Intensity", # Average nuclear intensity
				  "Nuclei: Compactness", # Characterizes shape. Calculated by 2*PI/area (gyration radius *gyration radius2). 
				  "Nuclei: Light Flux", #   Normalized amount of light emitted by the whole nuclei. It is equal to nucleus average intensity multiplied by area and normalized by cytoplasm average intensity. Nucleus are is taken in pixels. 				  
				  "Nuclei: IxA (Nuc)")  # Amount of light emited by nucleus. Is equal to nucleus average intensity multiplied by nucleus area (= IxA Nuc)				  
if (any(!featInterest %in% inputHeader)) {
    stop("Not all Features of Interest are Present in Input Data")
	}

inputData <- inputData[, inputHeader %in% c("Well", featInterest)]	
inputHeader <- inputHeader[inputHeader %in% c("Well", featInterest)]	


# (D) Filter out Non-Specific Stains
if (doNSFilter) {
	wellNames <- as.vector(unique(inputData[, 1]))
	idxExclAll <- NULL
	if	(makeQCPlot) {
		pdfFile <- paste(strsplit(inputFile,".csv")[[1]], "_QCGeom.pdf", sep = "")
		pdf(file = paste(resultsDir, pdfFile, sep =""))
		}
	for (i in 1:length(wellNames)) {
		well    <- wellNames[i]
		idxWell <- which(inputData[, 1] == well)
		idxExcl <- f_exclNonSpecific(wellData = inputData[idxWell, 3:14],
									 plotText = c(substr(c(stain1, stain2, stain3), 1, 8), well),
									 makePlot = makeQCPlot)
		idxExclAll <- c(idxExclAll, idxWell[idxExcl])
		}
	if	(makeQCPlot) dev.off()
	if  (length(idxExclAll) != 0) inputData <- inputData[-idxExclAll, ]
	}
	
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# (3) Calculate Count Z-scores
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

source(paste(scriptDir, funFile, sep = ""))

# (A) Determine Indices for Well / Stain and Number of Objects
wellNames <- as.vector(unique(inputData[, 1]))
wellIdxStain1 <- f_wellIdx(cbind(inputData[, 1], "Hoescht"), wellNames, "Hoescht")
# wellIdxStain1 <- f_wellIdx(inputData[, 1:2], wellNames, stain1)
# wellIdxStain2 <- f_wellIdx(inputData[, 1:2], wellNames, stain2)
# wellIdxStain3 <- f_wellIdx(inputData[, 1:2], wellNames, stain3)
wellNumStain1 <- f_wellNumObj(wellIdxStain1)
#wellNumStain2 <- f_wellNumObj(wellIdxStain2)
#wellNumStain3 <- f_wellNumObj(wellIdxStain3)

# (B) Calculate Count / Percent Z-scores
idxNegCtrl <-  which(wellNames %in% negCtrlWells)
zCntStain1  <- f_cnt(wellNumStain1, idxNegCtrl, useRobust)
# if (stain2CntPercent) {
	# zCntStain2  <- f_cnt(wellNumStain2 / wellNumStain1, 
							# idxNegCtrl, useRobust)
	# }else{
	# zCntStain2  <- f_cnt(wellNumStain2, idxNegCtrl, useRobust)
	# }
# if (stain3CntPercent) {
	# zCntStain3  <- f_cnt(wellNumStain3 / wellNumStain1, 
							# idxNegCtrl, useRobust)
	# }else{
	# zCntStain3  <- f_cnt(wellNumStain3, idxNegCtrl, useRobust)
	# }	
	
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# (4) Calculate Distribution Z-scores
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	
if (scaleDistBoot) {
	zDistStain1 <- f_cmpdWellZBootDist(wellIdxStain1, stain1, stain1Names, 
									   idxNegCtrl, inputHeader, inputData, distStatMethod)
	# zDistStain2 <- f_cmpdWellZBootDist(wellIdxStain2, stain2, stain2Names, 
									   # idxNegCtrl, inputHeader, inputData, distStatMethod)
	
	}else{	
	# (A) Compare each Negative Control Well to Population of  all Negative Control Wells w/o that well
	stain1 <- "Hoescht"
	stain1Names <- paste(stain1, inputHeader[3:8], sep = ":")
	negCtrlDistStain1 <- f_wellNegCtrlDistSimple(wellIdxStain1, idxNegCtrl, inputData[, -1], distStatMethod)
	#negCtrlDistStain1 <- f_wellNegCtrlDist(stain1, stain1Names, wellIdxStain1, idxNegCtrl, inputHeader, inputData, distStatMethod)
	#negCtrlDistStain2 <- f_wellNegCtrlDist(stain2, stain2Names, wellIdxStain2, idxNegCtrl, inputHeader, inputData, distStatMethod)	

	# (B) Calculate location and Spread of Distribution Statistics
	if (useRobust) {
		locDistStain1 <- apply(negCtrlDistStain1, 2, median)
		#locDistStain2 <- apply(negCtrlDistStain2, 2, median)
		
		spreadDistStain1 <- apply(negCtrlDistStain1, 2, mad)
		#spreadDistStain2 <- apply(negCtrlDistStain2, 2, mad)
		}else{
		locDistStain1 <- apply(negCtrlDistStain1, 2, mean)
		#locDistStain2 <- apply(negCtrlDistStain2, 2, mean)
		
		spreadDistStain1 <- apply(negCtrlDistStain1, 2, sd)
		#spreadDistStain2 <- apply(negCtrlDistStain2, 2, sd)	
		}
	
	if (!scaleDistFeat) {
		locDistStain1    <- rep(0,  ncol(negCtrlDistStain1))
		#locDistStain2    <- rep(0, length(stain2Names))
		spreadDistStain1 <- rep(1,  ncol(negCtrlDistStain1))
		#spreadDistStain2 <- rep(1, length(stain2Names))
		}
		
	# (C) Calculate Z-scores	
	zDistStain1 <- f_cmpdWellZDistSimple(wellIdxStain1, idxNegCtrl, inputData[, -1], 
	                                     distStatMethod, locDistStain1, spreadDistStain1)
	#zDistStain1 <- f_cmpdWellZDistSimple(stain1, stain1Names, wellIdxStain1, idxNegCtrl, inputHeader, 
	#						   inputData, distStatMethod, locDistStain1, spreadDistStain1)
	#zDistStain2 <- f_cmpdWellZDist(stain2, stain2Names, wellIdxStain2, idxNegCtrl, inputHeader, 
	#						   inputData, distStatMethod, locDistStain2, spreadDistStain2)
	}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# (5) Save Workspace
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	
resultsFile <- paste(strsplit(inputFile,".CSV")[[1]], resultsSuffix, ".Rdata", sep = "")
save.image(file = paste(resultsDir, resultsFile, sep = ""), compress = TRUE)
}
#quit(save="no")

