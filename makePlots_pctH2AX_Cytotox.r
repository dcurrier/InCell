library(RColorBrewer)
featInterest <- "Cells: IxA (Nuc)"


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 					Directories and File Names
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

baseDir		<- "U:/Transfer_TMP/BirgitJon/July31_2014"
scriptDir	<- paste(baseDir, "SourceScripts", sep = "/")
inputDir	<- paste(baseDir, "Input", sep = "/")
resultsDir	<- paste(baseDir, "Results", sep = "/")


inputFiles <- c("HCS000100056_CPC DNA H2AX 10X.CSV", "HCS000101483_CPC DNA H2AX 10X.CSV")

for (inputFile in inputFiles) {
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	
# Load Input Data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 					 				 

tmpHeader	<- read.csv(file = paste(inputDir, inputFile, sep = "/"),
						 header = FALSE, row.names=NULL,  stringsAsFactors = FALSE,
			 			 nrows = 2, skip = 18) 
# Read Input Data (2 steps:
#                   (1) All lines to figure out the last cell-level line
#                   (2) Use this specification to re-read data
inputData	<- read.csv(file = paste(inputDir, inputFile, sep = "/"),
						 header = FALSE, row.names = NULL, stringsAsFactors = FALSE,
						 skip = 20)
idxLastRow <- which(!grepl(" - ", inputData[,1]))[1] - 1
#idxLastRow <- which(inputData[,1] == "Well")[1] - 3
inputData	<- read.csv(file = paste(inputDir, inputFile, sep = "/"),
						 header = FALSE, row.names = NULL,  stringsAsFactors = FALSE,
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
tmpNames <- sub("(fld 3)","", tmpNames,fixed = TRUE)
tmpNames <- sub("(fld 4)","", tmpNames,fixed = TRUE)
inputData[, 1] <- tmpNames
gc()						
						 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	
# Restrict to Feature of Interest: Only Well Name and Feature
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 					 				 						 
						
# Column 1 == well
featDF <- inputData[, c(1, which(inputHeader == featInterest))] 						
names(featDF) <- c("well", "y")

# log10 transform y
featDF$y <- log10(featDF$y)

# Make well names consistent with plate Maps
idxFix <- which(nchar(featDF$well) == 2)
oldWells <- featDF$well[idxFix]
newWells <- sapply(oldWells, function(x) {
                 paste(strsplit(x, "")[[1]][1], "0", strsplit(x, "")[[1]][2], 
				       sep = "") })
featDF$well[idxFix] <- newWells				 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	
# Load Annotations
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 					 				 						 

plateMaps <- read.delim(paste(inputDir, "PlateSummary.txt", sep = "/"), 
                        stringsAsFactors = FALSE)
# Restrict to Plates Screened
plateMaps <- plateMaps[plateMaps$currplateid %in% c(500000206131, 500000220168 ), ]						

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	
# Use Pos and Neg Ctrl Wells to Determine Cutoff
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 					 				 						 						

yCuts <- pretty(featDF$y, 10)

# Look up negative control wells
negCtrlWells <- unique(featDF$well[!featDF$well %in% plateMaps$well])

# Define a Data Frame to Contain well names & Label
ctrlDF <- data.frame(well = negCtrlWells, label = "DMSO")

# Define Positive Controls and Look up well names
posCtrlCmpds <- unique(plateMaps$SAMPLE[plateMaps$currplateid == 500000220168])
for (thisCmpd in posCtrlCmpds) {
     cmpdDat <- plateMaps[plateMaps$currplateid == 500000220168 & plateMaps$SAMPLE == thisCmpd, ]
	 thisLabel <- strsplit(cmpdDat$synonym[1]," | ")[[1]][1]
	 # Select the Wells with Highest Concentration
	 cmpdWells <- cmpdDat$well[cmpdDat$CONC == max(cmpdDat$CONC)]
	 ctrlDF <- rbind(ctrlDF, data.frame(well = cmpdWells, label = thisLabel))
	 }

# Compute Percent H2AX + for all wells	 
pctMat <- matrix(NA, length(yCuts), nrow(ctrlDF))	 
for (i in seq(nrow(pctMat))) {
     thisCut <- yCuts[i] 
     for (j in seq(ncol(pctMat))) {
	      thisWell <- ctrlDF$well[j]
		  thisY    <- featDF$y[featDF$well == thisWell]
		  thisPct  <- mean(thisY > thisCut)
		  pctMat[i, j] <- thisPct
		}
	}
	
pctH2AXcut <- 5.5
	
# Make Plot
#resPrefix   <- strsplit(inputFile," ")[[1]][3]
resPrefix   <- strsplit(inputFile,"_")[[1]][1]
resultsFile <- paste(resPrefix, "Cutoff_Selection.pdf", sep = "_")

pdf(file = paste(resultsDir, resultsFile, sep = "/"))
plot(rep(seq(length(yCuts)), ncol(pctMat)), 100*pctMat, typ = "n", xaxt = "n",
     xlab = paste("log10", featInterest, "cutoff"), ylab = "% H2AX Positive",
 	 main = paste("Cutoff =", pctH2AXcut))
axis(1, seq(length(yCuts)), yCuts)
abline(v = seq(length(yCuts))[yCuts == pctH2AXcut], lwd = 2, lty = 2)
for (i in seq(length(unique(ctrlDF$label)))) {
     thisLabel <- unique(ctrlDF$label)[i]
     labelPctMat <- 100* pctMat[, ctrlDF$label == thisLabel]
	 lines(seq(yCuts), rowMeans(labelPctMat), lwd = 2, 
	       col = c("black", brewer.pal(4, "Paired")[c(2,4)])[i])
	 boxplot(t(labelPctMat), add = TRUE, xaxt = "n", outline = FALSE, col = "white",
	         border = c("black", brewer.pal(4, "Paired")[c(2,4)])[i])
	 points(jitter(rep(seq(length(yCuts)), ncol(labelPctMat))), labelPctMat,
            pch = 19, col = c("gray", brewer.pal(4, "Paired")[c(1,3)])[i]) 	 
	 }
legend("topright", legend = unique(ctrlDF$label), 
       fill = c("black", brewer.pal(4, "Paired")[c(1, 3) + 1]), bty = "n")			 
dev.off()	



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	
# Use Cutoff to Calculate %H2AX Positive
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 					 				 						 						


pctHP <- rep(NA, length(unique(featDF$well)))
for (i in seq(length(pctHP))) {
     thisWell <- unique(featDF$well)[i]
	 thisY    <- featDF$y[featDF$well == thisWell]
	 thisPct  <- mean(thisY > pctH2AXcut)
	  pctHP[i] <- thisPct	
	}
	
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 	
# Plot Concentration-Response for Non - Ctrl Compounds
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 					 				 						 							
resultsFile <- paste(resPrefix, "PctH2AX_and_Cytotox.pdf", sep = "_")
pdf(file = paste(resultsDir, resultsFile, sep = "/"))

yLimCnt <- range(pretty(table(featDF$well)))
yLimH2AX <- c(0, 50)	

# Loop over Each Compound Not on Ctrl Plate
unCmpds <- unique(plateMaps$SAMPLE[plateMaps$currplateid == 500000206131])
for (thisCmpd in unCmpds) {

	cmpdDat <- plateMaps[plateMaps$SAMPLE == thisCmpd & plateMaps$currplateid == 500000206131,]
	conc   <- cmpdDat$CONC /1000
	wells <- cmpdDat$well
	act   <- rep(NA, length(wells))
	n     <- rep(NA, length(wells))
	for (j in seq(length(conc))) {
		act[j] <- pctHP[which(unique(featDF$well) == wells[j])] * 100
		n[j]   <- sum(featDF$well == wells[j])
		}   

	# yLimCnt <- c(0, 500)	
	
	colCnt <- brewer.pal(4, "Set1")[2]
	colH2AX <- brewer.pal(4, "Set1")[1]

	plot(conc, n, pch = 19, type  = "b", 
		 main = paste(thisCmpd,strsplit(cmpdDat$synonym[1], " | ")[[1]][1]),
		 log = "x", bty = "n", xlab = "", ylab = "", axes = FALSE, 
		 ylim = yLimCnt, col = colCnt)
	axis(4, pretty(yLimCnt, 7), col=colCnt, col.ticks=colCnt,col.axis=colCnt)
	mtext("Cell Count", 4, line =2.5,col=colCnt)

	par(new=TRUE)
	plot(conc, act, pch=19, cex = 1.25, bty="n", axes=FALSE, type = "b",
		col=colH2AX, xlab="log10 Conc [uM]", lwd=2, log="x", ylab="", ylim=yLimH2AX)
	axis(2,pretty(yLimH2AX, 7), paste(pretty(yLimH2AX, 7), "%"), col=colH2AX,col.ticks=colH2AX,col.axis=colH2AX)
	axis(1,axTicks(1), log10(axTicks(1)), col='black',col.ticks='black',col.axis='black')				
	}
dev.off()	
}
