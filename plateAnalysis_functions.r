# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##
#	 FUNCTIONS TO PROCESS z-SCORES FROM A HIGH-CONTENT SCREENING RUN
# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##

# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##
f_combIdx <- function(x, idx){
	out <- NULL
	for (i in idx){
		idxI <- seq(x[i, 1], x[i, 2])
		out  <- c(out, idxI)
		}
	return(out)
	}
# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##

# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##
f_cvmts <- function(x, y) {
	# (1) Sort x and y
	xS <- sort(x)
	yS <- sort(y)
	# (2) Get Ranks within x and y
	i <- rank(xS)
	j <- rank(yS)
	xDF <- data.frame(val = xS, ID = "x")
	yDF <- data.frame(val = yS, ID = "y")
	# (3) Join x and y
	DF	<- rbind(xDF, yDF)
	# (4) Get Ranks within Combined Data
	r <- rank(DF$val)[which(DF$ID == "x")]
	s <- rank(DF$val)[which(DF$ID == "y")]
	# (5) Calculate N and M
	N <- length(x)
	M <- length(y)
	# (6) Calculate U
	U <- N*sum((r -i)^2) + M*sum((s - j)^2)
	# (7) Calculate T
	T <- U/N/M/(N+M) - (4*M*N-1)/(6*(M+N))
	# (8) Scale by Maximum Theoretical Value
	Tmax <- (2*M*N + 1)/(6*(M+N))
	Tscaled <- T / Tmax
	# (9) Correct Sign
	Tscaled <- Tscaled * sign(median(y) - median(x))
	return(Tscaled)
	}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #

## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #
f_wellNegCtrlDistSimple <- function (wellIdxStain, idxNegCtrl, inputData, distStatMethod) {

	if(distStatMethod == "KS")	  f_distStat <- f_ks
	if(distStatMethod == "CVMTS") f_distStat <- f_cvmts

	out <- matrix(NA, length(idxNegCtrl), ncol(inputData))
	for (i in 1:length(idxNegCtrl)) {
		idxWell	      <- seq(wellIdxStain[idxNegCtrl[i], 1], wellIdxStain[idxNegCtrl[i], 2])
		idxOtherWells <- f_combIdx(wellIdxStain, idxNegCtrl[-i])
		for (j in 1:ncol(inputData)) {
			distStat <- f_distStat(inputData[idxWell, j], inputData[idxOtherWells, j])
			out[i, j] <- distStat
		}
	}

	return(out)
	}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #

## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #
f_wellNegCtrlDist <- function (stain, stainNames, wellIdxStain, idxNegCtrl, inputHeader, inputData, distStatMethod) {

	if(distStatMethod == "KS")	  f_distStat <- f_ks
	if(distStatMethod == "CVMTS") f_distStat <- f_cvmts

	out <- matrix(NA, length(idxNegCtrl), length(stainNames))
	for (i in 1:length(idxNegCtrl)) {
		idxWell	      <- seq(wellIdxStain[idxNegCtrl[i], 1], wellIdxStain[idxNegCtrl[i], 2])
		idxOtherWells <- f_combIdx(wellIdxStain, idxNegCtrl[-i])
		idxFeatStain  <- f_getIdx(stain, stainNames, inputHeader)
		for (j in 1:length(idxFeatStain)) {
			idxFeat	 <- idxFeatStain[j]
			distStat <- f_distStat(inputData[idxWell, idxFeat], inputData[idxOtherWells, idxFeat])
			out[i, j] <- distStat
		}
	}

	return(out)
	}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #

## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #
f_cmpdWellZDist <- function (stain, stainNames, wellIdxStain, idxNegCtrl, inputHeader, inputData, distStatMethod,
							locDist, spreadDist) {

	if(distStatMethod == "KS")  f_distStat <- f_ks
	if(distStatMethod == "CVMTS") f_distStat <- f_cvmts

	idxCmpdWells <- seq(nrow(wellIdxStain))[-idxNegCtrl]
	out <- matrix(NA, length(idxCmpdWells), length(stainNames))
	for (i in 1:(length(idxCmpdWells))) {
		idxWell      <- seq(wellIdxStain[idxCmpdWells[i], 1], wellIdxStain[idxCmpdWells[i], 2])
		idxOtherWells <- f_combIdx(wellIdxStain, idxNegCtrl)
		idxFeatStain  <- f_getIdx(stain, stainNames, inputHeader)
		for (j in 1:length(idxFeatStain)) {
			idxFeat <- idxFeatStain[j]
			distStat <- f_distStat(inputData[idxWell, idxFeat], inputData[idxOtherWells, idxFeat])
			out[i, j] <- (distStat - locDist[j]) / spreadDist[j]
			}
		}

	return(out)
	}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #


## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #
f_cmpdWellZDistSimple <- function (wellIdxStain, idxNegCtrl, inputData, distStatMethod,
							locDist, spreadDist) {

	if(distStatMethod == "KS")  f_distStat <- f_ks
	if(distStatMethod == "CVMTS") f_distStat <- f_cvmts

	idxCmpdWells <- seq(nrow(wellIdxStain))[-idxNegCtrl]
	out <- matrix(NA, length(idxCmpdWells), ncol(inputData))
	for (i in 1:(length(idxCmpdWells))) {
		idxWell      <- seq(wellIdxStain[idxCmpdWells[i], 1], wellIdxStain[idxCmpdWells[i], 2])
		idxOtherWells <- f_combIdx(wellIdxStain, idxNegCtrl)
		for (j in 1:ncol(inputData)) {
			distStat <- f_distStat(inputData[idxWell, j], inputData[idxOtherWells, j])
			out[i, j] <- (distStat - locDist[j]) / spreadDist[j]
			}
		}

	return(out)
	}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #


## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #
f_ks <- function(x, y) {
    ksTwo <- as.vector(ks.test(x, y)$statistic)
    ksGreater <- as.vector(ks.test(x, y, alternative = "greater")$statistic)
    if (ksTwo == ksGreater) {
        ks <- ksTwo
    }else{
        ks <- -ksTwo
    }
    return(ks)
}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #

## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #
f_bootDistStat <- function(x, n, f_distStat){
	# Input: Data to bootstrap [x] a vector
	#		Sample Size [n] a scalar

	idx <- sample(length(x), n, replace = TRUE);
	out <- f_distStat(x[idx], x[-idx])
	return(out)
	}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #

## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #
f_cmpdWellZBootDist <- function (wellIdxStain, stain, stainNames, idxNegCtrl, inputHeader, inputData, distStatMethod) {
	#  Use Bootstrap resampling to estimate standard deviation of distribution statistic as f(sample size)
	rowIdx  <- f_combIdx(wellIdxStain, idxNegCtrl)
	colIdx 	<- f_getIdx(stain, stainNames, inputHeader)
	negData <- inputData[rowIdx, colIdx]
	if(distStatMethod == "KS")	  f_distStat <- f_ks
	if(distStatMethod == "CVMTS") f_distStat <- f_cvmts
	nVec    <- seq(10, 2000, length.out = 20)
	sdMat	<- matrix(NA, length(nVec), length(colIdx))
	for (i in 1:length(nVec)) {
		bootDat <- matrix(NA, 200, length(colIdx))
		for (j in 1:200) bootDat[j, ] <-  apply(negData, 2, f_bootDistStat, nVec[i], f_distStat)
		sdMat[i, ]  <- apply(bootDat, 2, sd)
		}

	# Calculate Z-sores by scaling Cmpd wells to SD
	idxCmpdWells <- seq(nrow(wellIdxStain))[-idxNegCtrl]
	out <- matrix(NA, length(idxCmpdWells), length(colIdx))
	for (i in 1:(length(idxCmpdWells))) {
		idxWell   <- seq(wellIdxStain[idxCmpdWells[i], 1], wellIdxStain[idxCmpdWells[i], 2])
		for (j in 1:length(colIdx)) {
			idxFeat <- colIdx[j]
			distStat <- f_distStat(inputData[idxWell, idxFeat], inputData[rowIdx, idxFeat])
			out[i, j] <- distStat/approx(nVec, sdMat[, j], length(idxWell))$y
			}
		}
	return(out)
	}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #


# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##
# Utility function to Identify the columns of interest for a given raw data file
f_getIdx <- function (Stain, StainFeatures, data.colnames){
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	# Input:
	#   Stain - Name of Stain (Hoescht, H2AX, etc)
	#	StainFeatures - Names of Features of Interest for this stain [Width / Depth/ Major Axis Length]
	#	data.colnames - Names of the Columns in the Raw data Files
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

	featureIdx <- rep(NA, length(StainFeatures))
	for(feature.num in 1:length(StainFeatures)){
		this.feature <- StainFeatures[feature.num]
		# Return the Index of the Feature that has Stain and this.feature by Using Grep
		this.idx <- which(sapply(data.colnames, function(x) grepl(Stain, as.vector(x))) &
						  sapply(data.colnames, function(x) grepl(this.feature, as.vector(x))))
		featureIdx[feature.num] <- as.vector(this.idx)
		}
	return(featureIdx)
	}
# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##

# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##
# Function To Link a Given Hoescht Object to One [2nd stain] and One [third stain] Object or NULL Using the Objects' Locations
f_overlap <- function(i, geomAnchor, geomSecondary, geomTertiary){

		thisAnchor <- as.numeric(geomAnchor[i, ])
		thisAnchor[3:4] <- thisAnchor[3:4]/1000/2
		geomSecondary[, 3:4] <- geomSecondary[, 3:4]/1000/2
		geomTertiary[, 3:4] <- geomTertiary[, 3:4]/1000/2

		# Objects intersect if any of the following conditions are NOT TRUE:
		#Note that we're correcting for the factor of 2 Above:
		# (1) Xcaspase - WidCaspase/2 > Xhoescht + WidHoescht/2 (Caspase object is to the Right of Hoescht Cell)
		# (2) Xcaspase + WidCaspase/2 < Xhoescht - WidHoescht/2 (Caspase object is to the Left of Hoescht Cell)
		# (3) Ycaspase - DepthCaspase/2 > Yhoescht + DepthHoescht/2 (Caspase object is Above the Hoescht Cell)
		# (4) Ycaspase + DepthCaspase/2 < Yhoescht - DepthHoescht/2 (Caspase object is Below the Hoescht Cell)

		Stain2Overlap<- which(!(geomSecondary[,1] - geomSecondary[,3] > thisAnchor[1] + thisAnchor[3] |
								geomSecondary[,1] + geomSecondary[,3] < thisAnchor[1] - thisAnchor[3] |
								geomSecondary[,2] - geomSecondary[,4] > thisAnchor[2] + thisAnchor[4] |
								geomSecondary[,2] + geomSecondary[,4] < thisAnchor[2] - thisAnchor[4]))

		Stain3Overlap <- which(!(geomTertiary[,1] - geomTertiary[,3] > thisAnchor[1] + thisAnchor[3] |
								geomTertiary[,1] + geomTertiary[,3] < thisAnchor[1] - thisAnchor[3] |
								geomTertiary[,2] - geomTertiary[,4] > thisAnchor[2] + thisAnchor[4] |
								geomTertiary[,2] + geomTertiary[,4] < thisAnchor[2] - thisAnchor[4]))

		# If empty, replace with NA
		if(length(Stain2Overlap) == 0){ Stain2Overlap <- NA };
		if(length(Stain3Overlap) == 0){ Stain3Overlap <- NA };

		if(length(Stain2Overlap) > length(Stain3Overlap)){
			Stain3Overlap <- c(Stain3Overlap,rep(NA,length(Stain2Overlap) - length(Stain3Overlap)))};

		if(length(Stain3Overlap) > length(Stain2Overlap)){
					Stain2Overlap <- c(Stain2Overlap,rep(NA,length(Stain3Overlap) - length(Stain2Overlap)))};

		# 1st Column is the Hoescht Index (i); 2nd is the matching Stain2 Overlap; 3rd is the matching H3
		return(cbind(i, Stain2Overlap, Stain3Overlap));
		};
# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##

# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##
f_exclStain <- function (stainNum, overlap) {
		# (1) Select Hoescht [col 1] and Stain of Interest [Col 2]
		overlap <- overlap[, c(1, stainNum)]
		# (2) Remove Hoescht objects [rows] that Overlap No Objects of Interest
		overlap <- as.matrix(na.omit(overlap))
		# (2A) Consider Special Cases: No Overlap or One Overlap [Therefore specific]
		if (nrow(overlap) == 1) idxOverlap <- overlap[2]
		if (nrow(overlap) == 0) idxOverlap <- NULL
		# (3) Identify Hoescht Objects intersecting > 1 ObjOfInterest [idxPromAnchor]
		#		and ObjsofInterst intersecting > 1 Hoescht Object [idxPromStain]
		if (nrow(overlap) > 1) {
			idxPromAnch  <- unique(overlap[duplicated(overlap[, 1]), 1])
			idxPromStain <- unique(overlap[duplicated(overlap[, 2]), 2])
			overlap <- overlap[!overlap[, 1] %in% idxPromAnch, ]
			if (!is.null(nrow(overlap))) {
				overlap <- overlap[!overlap[, 2] %in% idxPromStain, ]
				}else{
				overlap <- overlap[!overlap[2] %in% idxPromStain]
				}
			if (!is.null(nrow(overlap))) {
				idxOverlap <- overlap[, 2]
				}else{
				idxOverlap <- overlap[2]
				}
			}
		return(idxOverlap)
	}
# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##

# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##
f_plotBox <- function(geom, boxcol, addCrossline) {
	# Top Line
	lines(x = geom[1] + geom[3]/1000*c(-1/2, 1/2), y = geom[2] + geom[4]/1000*c(1/2, 1/2), col=boxcol)
	# Bottom Line
	lines(x = geom[1] + geom[3]/1000*c(-1/2, 1/2), y = geom[2] - geom[4]/1000*c(1/2, 1/2), col=boxcol)
	# Left Line
	lines(x = geom[1] - geom[3]/1000*c(1/2, 1/2), y = geom[2] + geom[4]/1000*c(-1/2, 1/2), col=boxcol)
	# Right Line
	lines(x = geom[1] + geom[3]/1000*c(1/2, 1/2), y = geom[2] + geom[4]/1000*c(-1/2, 1/2), col=boxcol)
	# CrossLines
	if (addCrossline) {
		lines(x= geom[1] + geom[3]/1000*c(-1/2, 1/2), y = geom[2] + geom[4]/1000*c(1/2, -1/2), col=boxcol)
		lines(x= geom[1] + geom[3]/1000*c(-1/2, 1/2), y = geom[2] + geom[4]/1000*c(-1/2, 1/2), col=boxcol)
		}
	}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##
f_plotwellGeom <- function(geomStain1, geomStain2, geomStain3, idxSpecStain2, idxSpecStain3, plotText) {

	# Define Plot Area: Black Background , with white axes and Main
	par(bg="black",fg="white",col.axis="white",col.lab="white",col.main = "white")
	plot(geomStain1[,1], geomStain1[,2], type="n",
	     pch=19,cex=0.125,col='blue',main=plotText[4],
	     xlab = "millimeters",ylab = "millimeters"); #xlim=c(1.5,2),ylim=c(0.5,0.8)

	# Annotate the Stains
	mtext(plotText[1], side=3, line=0.5, col='blue') # Stain 1 Name
	mtext(plotText[3], side=3, line=0.5, col='green', at=1.5) # Stain 3 Name
	mtext("Unspecific ", side=3, line=0.5,col='yellow', at=1)
	mtext(plotText[2], side=3, line=0.5, col='red',at=2.4) # Stain 2 name
	mtext("Unspecific", side=3, line=0.5, col='pink',at=2.9)

	# Plot Hoescht
	apply(geomStain1, 1, f_plotBox, boxcol = "blue", addCrossline = FALSE)

	# Plot H2AX
	if (is.null(idxSpecStain3)) {
		# No Specific Stains; All Non-Specific
		apply(geomStain3, 1, f_plotBox, boxcol = "yellow", addCrossline = TRUE)
		}else{
		if (length(idxSpecStain3) == nrow(geomStain3)) {
			# All Speciic Staining
			apply(geomStain3, 1, f_plotBox, boxcol = "green", addCrossline = TRUE)
			}else{
			# Some Specific, some non-specific
			apply(geomStain3[idxSpecStain3, ], 1, f_plotBox, boxcol = "green", addCrossline = TRUE)
			apply(geomStain3[-idxSpecStain3, ], 1, f_plotBox, boxcol = "yellow", addCrossline = TRUE)
			}
		}
	# Plot Tubulin
	if (is.null(idxSpecStain2)) {
		# No Specific Stains; All Non-Specific
		apply(geomStain2, 1, f_plotBox, boxcol = "pink", addCrossline = TRUE)
		}else{
		if (length(idxSpecStain2) == nrow(geomStain2)) {
			# All Speciic Staining
			apply(geomStain2, 1, f_plotBox, boxcol = "red", addCrossline = TRUE)
			}else{
			# Some Specific, some non-specific
			apply(geomStain2[idxSpecStain2, ], 1, f_plotBox, boxcol = "red", addCrossline = TRUE)
			apply(geomStain2[-idxSpecStain2, ], 1, f_plotBox, boxcol = "pink", addCrossline = TRUE)
			}
		}
	}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
f_exclNonSpecific <- function (wellData, plotText, makePlot) {
	idxStain1 <- which(!is.na(wellData[, 1]))
	idxStain2 <- which(!is.na(wellData[, 5]))
	idxStain3 <- which(!is.na(wellData[, 9]))

	geomStain1 <- wellData[idxStain1, 1:4]
	geomStain2 <- wellData[idxStain2, 5:8]
	geomStain3 <- wellData[idxStain3, 9:12]

	# Loop Over Each Anchor Object (Stain 1) to Find overlap
	overlap <- NULL
	for (i in 1:nrow(geomStain1)) {
		thisOverlap <- f_overlap(i, geomStain1, geomStain2, geomStain3)
		overlap <- rbind(overlap, thisOverlap)
		}

	# Identify Secondary and Tertiary Objects that have 1-1 overlap with Anchor
	idxSpecStain2 <- f_exclStain(2, overlap)
	idxSpecStain3 <- f_exclStain(3, overlap)

	# Get the Indices Identifying Non-specific Stains
	if	(!is.null(idxSpecStain2)) {
		idxNonSpecStain2 <- idxStain2[-idxSpecStain2]
		}else{
		idxNonSpecStain2 <- idxStain2
		}
	if	(!is.null(idxSpecStain3)) {
		idxNonSpecStain3 <- idxStain3[-idxSpecStain3]
		}else{
		idxNonSpecStain3 <- idxStain3
		}

	# Make Geometry Plot
	if (makePlot) {
		f_plotwellGeom(geomStain1, geomStain2, geomStain3,
					   idxSpecStain2, idxSpecStain3, plotText)
		}

	# Return the Indices Identifying Non-specific Stains
	return(c(idxNonSpecStain2, idxNonSpecStain3))
	}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
f_wellIdx <- function(wellScan, wellNames, stain) {
	out <- matrix(NA, length(wellNames), 2)
	for (i in 1:length(wellNames)) {
		well <- wellNames[i]
		idx	 <- which(wellScan[, 1] == well & wellScan[, 2] == stain)
		if	(length(idx) != 0)	out[i, ] <- range(idx)
		}
	return(out)
	}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
f_wellNumObj <- function(x) {
	numObj <- x[, 2] - x[, 1] + 1
	numObj[which(is.na(numObj))] <- 0
	return(numObj)
	}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
f_cnt <- function(wellNum, idxNegCtrl, useRobust) {
	negCtrl <- wellNum[idxNegCtrl]
	if (useRobust){
		loc     <- median(negCtrl)
		spread  <- mad(negCtrl)
		}else{
		loc     <- mean(negCtrl)
		spread  <- sd(negCtrl)
		}
	cnt <- (wellNum[-idxNegCtrl] - loc)/spread
	return(cnt)
	}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #