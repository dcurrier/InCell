library(drc)

if (centerRobust){
	f_center <- median
	}else{
	f_center <- mean
	}
	

f_aic <- function(model,n,k) {
  2*k - 2*(as.vector(logLik(model))) + 2*k*(k + 1)/(n - k - 1)
	}

f_fitDR <- function(concResp) {							
	# ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # 
	# Changelog:
	# 	- Imposed Constraint that the EC50 must lie within the concentration range tested
	#   - If a fit doesn't converge, Select the best of the surviving fits
	#	- Added a Linear Fit to Compare
	#	- Removed the Cases where Either 0 or 1 fit Converged
	#   - Added a Check at the end to Compare AIC(fit) to AIC(fitLinear)
	#	 -Use the corrected AIC function Anang had written
	# ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # 
	
	# ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # 
	# Fit 1 is the fit to a 4-parameter log-logisitic function with constraints on the EC50 [in log units]
	fit1 <- drm(y ~ x, data = concResp, fct = LL2.4(fixed=c(NA,NA,NA,NA)),control = drmc(errorm = FALSE,otrace=TRUE),
				lowerl = c(-Inf,-Inf,-Inf,log(min(concResp$x))), upperl = c(Inf,Inf,Inf,log(max(concResp$x))));
	
	# Fit 2 adds fixed values for our best guesses for Max and Min Response [@the lowest conc, response is 0; at the hightes conc; and use the response @ highest conc];	             
	fit2 <- drm(y ~ x, data = concResp, fct = LL2.4(fixed=c(NA,0,concResp$y[which.max(concResp$x)],NA)),control = drmc(errorm = FALSE,otrace=TRUE),
				lowerl = c(-Inf,log(min(concResp$x))), upperl = c(Inf,log(max(concResp$x))));
	
	# Fit 3 is the fit to a 4-parameter log-logistic model with constraints on the EC50 [in log units]
	fit3 <- drm(y ~ x, data = concResp, fct = LL2.3(fixed=c(NA,NA,NA)),control = drmc(errorm = FALSE,otrace=TRUE),
				lowerl = c(-Inf,-Inf,log(min(concResp$x))), upperl = c(Inf,Inf,log(max(concResp$x))));
	
	# Fit Linear is the Fit to a Horizontal, Concentration-Independent Response
	fitLinear <- lm(y ~ 1, data = concResp)
	# ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # 
	
	# ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # 
	# Count How Many Fits Converged
	numCovergedFits <- length(which(c(is.null(fit1$convergence), is.null(fit2$convergence), is.null(fit3$convergence))))
	# If all Fits Converged, Select by AIC
	if(numCovergedFits == 3){
		AICs <- c(f_aic(fit1, nrow(concResp), 4), f_aic(fit2, nrow(concResp), 2), f_aic(fit3, nrow(concResp), 3))
		bestFit <- which.min(AICs);
		if(bestFit == 1){
			fit <- fit1;
			}else{
				if(bestFit == 2){
					fit <- fit2;
				}else{
				fit <- fit3;
				};
			};
		}
	
	# If Two Fits Converged, Select the better of the Two by AIC
	if(numCovergedFits == 2){
		# Fit 1 Failed to Converge
		if(!is.null(fit1$convergence)){  
			AICs <- c(f_aic(fit2, nrow(concResp), 2), f_aic(fit3, nrow(concResp), 3))
			bestFit <- which.min(AICs);
			if(bestFit == 1) fit <- fit2
			if(bestFit == 2) fit <- fit3
			}
		# Fit 2 Failed to Converge
		if(!is.null(fit2$convergence)){  
			AICs <- c(f_aic(fit1, nrow(concResp), 4), f_aic(fit3, nrow(concResp), 3))
			bestFit <- which.min(AICs);
			if(bestFit == 1) fit <- fit1
			if(bestFit == 2) fit <- fit3
			}	
		# Fit 3 Failed to Converge
		if(!is.null(fit3$convergence)){  
			AICs <- c(f_aic(fit1, nrow(concResp), 4), f_aic(fit2, nrow(concResp), 2))
			bestFit <- which.min(AICs);
			if(bestFit == 1) fit <- fit1
			if(bestFit == 2) fit <- fit2
			}
		}
	# Compare to a Linear fit
	if(numCovergedFits > 0){
		if(f_aic(fitLinear, nrow(concResp), 1) < min(AICs)){
			fit <- fitLinear
			}
		}else{ # No fit converged, take the linear fit
		fit <- fitLinear # 
		}
	# ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # # ### # 
	 return(fit)
	 }
	 

f_evalFit <- function(fit, rSquaredMin, concResp){
	# ChangeLog
	# [August 10, 2012]
	#	- Removed Excessive Comments [given in previous version of functions_EN_OverParameterSpace.r]
	#	- After Convergence, Check, Fit a Linear Model
	#   - Removed nPassed.min input and DMSOstats.in
	
	
	# Did the fit Converge?
	fitConverged <- is.null(fit$convergence);
	
	# Is the Returned Fit Linear?
	fit.linear <- length(coef(fit)) == 1
	
	if(fitConverged & !fit.linear){
		
		# Calculate Rsquared
		fitPredict <- fitted(fit);
		ssTot <- sum((fit$data[,2] - mean(fit$data[,2]))^2);
		ssErr <- sum((fit$data[,2] - fitPredict)^2);
		rSquared <- (1 - (ssErr/ssTot));
		
		# Is the compound response within the linear range of the Dose-response curve?
		# could test this by setting a cutoff for the slope of the first/last points or 
		# -or- by determining the span of EC values 
		#x.interp.range <- seq(2.5,97.5,2.5); # get ec 2.5 to ec 97.5 by 2.5 
		#ECs <- as.vector(exp(ED(fit,x.interp.range,display=FALSE)[,1]));
		#EC.high <- approx(ECs, x.interp.range, max(fit$data$x)); # The EC corresponding to the max cmpd concentration
		#EC.low <- approx(ECs, x.interp.range, min(fit$data$x)); # The EC corresponding to the min cmpd concentration
		EC20 <- exp(ED(fit,20,display=FALSE)[1]);
		EC50 <- exp(ED(fit,50,display=FALSE)[1]);
		EC80 <- exp(ED(fit,80,display=FALSE)[1]);
		
		if(rSquared < rSquaredMin | min(fit$data$x) > EC20 | max(fit$data$x) < EC80 | 
			EC50 == min(fit$data$x) | EC50 == max(fit$data$x)){ # Active, not curve-like; cmpd response may lie in linear region
			efficacy <- fit$data$y[which.max(fit$data$x)] - fit$data$y[which.min(fit$data$x)]; # Response@HighestConc - Response@LowestConc
			potency <- -log10(max(fit$data[,1])/10^6); # Set potency to be the highest concentration tested [convert to -log10 and from uM to Molar
			}else{  # Active and Acceptable fit
				efficacy <- PR(fit,1e20)[[1]] - PR(fit,1e-20)[[1]]; # Efficacy is response at a High Conc - response at a Low Conc
				potency <- -log10(EC50/10^6); # -log10 of the ec50 value, [convert from uM to M]
		
			};
		};
	
	# Fit Did Not Converge, but Is Not Linear
	if(!fitConverged & !fit.linear){ 
		efficacy <- concResp$y[which.max(concResp$x)] - concResp$y[which.min(concResp$x)] # use concResp; efficacy is Highest - Lowest Conc response
		potency <- -log10(max(concResp$x)/10^6); # Set potency to be the highest concentration tested [convert to -log10 and from uM to Molar
		}
	
	# Fit was Linear
	if(fit.linear){
		potency <- -log10(max(concResp$x)/10^6); # Set potency to be the highest concentration tested [convert to -log10 and from uM to Molar
		efficacy <- 0
		}
	
	output <- c(efficacy,potency);
	return(output)
};


f_procPlateMap <- function (zFileNames, plateMap, plateAlias) {
	# Inputs : zFileNames -  a vector of filenames corresponding to
	#							 zscores from replicated plates
	#		  plateMap - the plate map containing the annotation
	
	# (1) Collate Zscores from Replicate Plates and
	#	 Set Distribution Features with Insufficient Objct Number to NA
	repCntStain1 <- NULL
	repCntStain2 <- NULL
	repCntStain3 <- NULL
	repDistStain1  <- NULL
	repDistStain2  <- NULL
	repDistStain3  <- NULL
	for (fileName in zFileNames) {
		load(fileName)
		#  Set Distribution Features with Insufficient Objct Number to NA
		zDistStain1[which(wellNumStain1[-idxNegCtrl] < numObjMin), ] <- NA
		zDistStain2[which(wellNumStain2[-idxNegCtrl] < numObjMin), ] <- NA
		# Collate Z-scores
		repCntStain1 <- cbind(repCntStain1, zCntStain1)
		repCntStain2 <- cbind(repCntStain2, zCntStain2)
		repCntStain3 <- cbind(repCntStain3, zCntStain3)
		repDistStain1 <- cbind(repDistStain1, zDistStain1)	
		repDistStain2 <- cbind(repDistStain2, zDistStain2)	
		}
	numRep <- ncol(repCntStain1)
	
	# (2) Make Cmpd Well Names Consistent
	cmpdWells  <- wellNames[-idxNegCtrl]
	idxFixName <- which(nchar(cmpdWells) < 3)
	splitNames <- sapply(cmpdWells[idxFixName], function(x) strsplit(x,"")[[1]])
	fixedNames <- paste(splitNames[1,], splitNames[2,], sep="0")
	cmpdWells[idxFixName] <- fixedNames
	
	# (3) Loop Over Each Compound
    cmpds       <- unique(as.vector(plateMap$SAMPLE))
	cmpdAliases <- cmpds
    numScores   <- 2*(3 + ncol(zDistStain1) + ncol(zDistStain2))
	scoredDR    <- matrix(NA, nrow = length(cmpds), ncol = numScores)
	featNames   <- rep(NA, length(numScores))
	for (cmpd in cmpds) {
		cmpdNum   <- which(cmpds == cmpd)
		cmpdMap   <- plateMap[as.vector(plateMap$SAMPLE) %in% cmpd, ]
		idxWells  <- which(cmpdWells %in% as.vector(cmpdMap$well))
		cmpdConc  <- cmpdMap$CONC / 1000 
		cmpdAlias <- as.vector(plateAlias$alias[which(plateAlias$regno %in% cmpd)])
		if (length(cmpdAlias) != 0) cmpdAliases[cmpdNum] <- cmpdAlias
		
		# Plot 1: Annotation
		plot(5,axes=FALSE, type="n", main = paste(cmpd, cmpdAlias, sep = "\n"),
		     xlab="", ylab="", cex.main=1.5, bty="n")	
		
		# Fit Count Stain 1
		concResp               <- data.frame(y = apply(repCntStain1, 1, f_center)[idxWells], x = cmpdConc)
		cntFit1                <- f_fitDR(concResp)
		paramCnt1              <- f_evalFit(cntFit1, rSquaredMin, concResp)
		featureName            <- paste("No. [", stain1, "] Objects", sep = "")
		featNames[1:2]         <- paste(featureName, c("Efficacy","Potency"),sep = " - ")
		scoredDR[cmpdNum, 1:2] <- paramCnt1
		# Loop Over the Additional Features
		for (j in 1:(2 + ncol(zDistStain1) + ncol(zDistStain2))) {
			
			if (j == 1) { # Cnts Stain 2	
				featureName <- ifelse(stain2CntPercent, paste("%[", stain2, "]+", sep = ""), 
														paste("No. [", stain2, "]+", sep = ""))
				resp <- repCntStain2[idxWells, ]
				concResp <- data.frame(x = cmpdConc, y = apply(resp, 1, f_center))
				fit  <- f_fitDR(concResp)
				yLim <- range(repCntStain2)
				scoredParam <- f_evalFit(fit, rSquaredMin, concResp)
				}
				
			if (j == 2) { # Cnts Stain 3
				featureName <- ifelse(stain3CntPercent, paste("%[", stain3, "]+", sep = ""), 
														paste("No. [", stain3, "]+", sep = ""))
				resp <- repCntStain3[idxWells, ]
				concResp <- data.frame(x = cmpdConc, y = apply(resp, 1, f_center))
				fit  <- f_fitDR(concResp)
				yLim <- range(repCntStain3)
				scoredParam <- f_evalFit(fit, rSquaredMin, concResp)
				}
				
			if (j > 2 & j < (2 + length(stain1Names) + 1)) {
				# Stain 1 Dist Features	
				stainIdx <- length(stain1Names)*c(0, 1, 2) + (j - 2)
				featureName <- paste(stain1, stain1Names[j - 2], sep = ":")
				resp <- repDistStain1[idxWells, stainIdx]
				
				# If Only One Replicate has data, set to NA
				resp[which(apply(!is.na(resp), 1, sum) == 1), ] <- NA 
							
				concResp <- data.frame(x = cmpdConc, y = apply(resp, 1, f_center, na.rm = TRUE))
				concResp <- na.omit(concResp)
				yLim <- range(repDistStain1[, stainIdx], na.rm = TRUE)
		
				# If <2 concentration has a non-NA response, Efficacy = 0, potency = highest conc
				if (length(which(!is.nan(rowMeans(resp,na.rm=TRUE)))) > 1) {	
					fit  <- f_fitDR(concResp)
					scoredParam <- f_evalFit(fit, rSquaredMin, concResp)
					}else{ # Fit Response
					fit  <- NULL
					scoredParam <- c(0, -log10(max(cmpdConc)/10^6))
					}
				}
				
			if (j > (2 + length(stain1Names))) { 
				# Stain 2 Dist Features	
				stainIdx <- length(stain2Names)*c(0, 1, 2) + (j - 2 - length(stain1Names))
				featureName <- paste(stain2, stain2Names[j - 2 - length(stain1Names)], sep = ":")
				resp <- repDistStain2[idxWells, stainIdx]
				
				# If Only One Replicate has data, set to NA
				resp[which(apply(!is.na(resp), 1, sum) == 1), ] <- NA 
				
				concResp <- data.frame(x = cmpdConc, y = apply(resp, 1, f_center, na.rm = TRUE))
				concResp <- na.omit(concResp)
				yLim <- range(repDistStain2[, stainIdx], na.rm = TRUE)
				
				# If <2 concentration has a non-NA response, Efficacy = 0, potency = highest conc
				if (length(which(!is.nan(rowMeans(resp,na.rm=TRUE)))) > 1) {	
					fit  <- f_fitDR(concResp)
					scoredParam <- f_evalFit(fit, rSquaredMin, concResp)
					}else{ # Fit Response
					fit  <- NULL
					scoredParam <- c(0, -log10(max(cmpdConc)/10^6))
					}
				}	
				
				# Make Plot
				f_plotRightAxis(cmpdConc, repCntStain1[idxWells, ], cntFit1, range(repCntStain1), "Cell Count", "gray")
				f_plotLeftAxis(cmpdConc, resp, fit, yLim, scoredParam, featureName, "black")
			
				# Assign Outputs
				featNames[1:2 + j*2] <- paste(featureName, c("Efficacy","Potency"),sep = " - ")
				scoredDR[cmpdNum, 1:2 + j*2] <- scoredParam
				}
			}
		return(list(featNames = featNames, scoredDR = scoredDR, cmpds = cmpds, cmpdAliases = cmpdAliases))
		}

f_procPlateMapSimple <- function (zFileNames, plateMap, plateAlias) {
	# Consider only One Stain
	# Inputs : zFileNames -  a vector of filenames corresponding to
	#							 zscores from replicated plates
	#		  plateMap - the plate map containing the annotation
	# 
	# (1) Collate Zscores from Replicate Plates and
	#	 Set Distribution Features with Insufficient Objct Number to NA
	repCntStain1 <- NULL
	repDistStain1  <- NULL
	for (fileName in zFileNames) {
		load(fileName)
		#  Set Distribution Features with Insufficient Objct Number to NA
		zDistStain1[which(wellNumStain1[-idxNegCtrl] < numObjMin), ] <- NA
		# Collate Z-scores
		repCntStain1 <- cbind(repCntStain1, zCntStain1)
		repDistStain1 <- cbind(repDistStain1, zDistStain1)	
		}
	numRep <- ncol(repCntStain1)
	
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
    numScores   <- 2*(1 + ncol(zDistStain1))
	scoredDR    <- matrix(NA, nrow = length(cmpds), ncol = numScores)
	featNames   <- rep(NA, numScores)
	for (cmpd in cmpds) {
		cmpdNum   <- which(cmpds == cmpd)
		# 05.01.2014
		# cmpdMap   <- plateMap[as.vector(plateMap$SAMPLE) %in% cmpd, ]
		cmpdMap   <- plateMap[as.vector(plateMap$SAMPLE) %in% cmpd & 
		                      as.vector(plateMap$well) %in% cmpdWells, ]
		idxWells  <- which(cmpdWells %in% as.vector(cmpdMap$well))
		cmpdConc  <- cmpdMap$CONC / 1000 
		cmpdAlias <- as.vector(plateAlias$alias[which(plateAlias$regno %in% cmpd)])
		if (length(cmpdAlias) != 0) cmpdAliases[cmpdNum] <- cmpdAlias
		
		# Plot 1: Annotation
		plot(5,axes=FALSE, type="n", main = paste(cmpd, cmpdAlias, sep = "\n"),
		     xlab="", ylab="", cex.main=1.5, bty="n")	
		
		# Fit Count Stain 1
		concResp               <- data.frame(y = apply(repCntStain1, 1, f_center)[idxWells], x = cmpdConc)
		cntFit1                <- f_fitDR(concResp)
		paramCnt1              <- f_evalFit(cntFit1, rSquaredMin, concResp)
		featureName            <- paste("No. [", stain1, "] Objects", sep = "")
		featNames[1:2]         <- paste(featureName, c("Efficacy","Potency"),sep = " - ")
		scoredDR[cmpdNum, 1:2] <- paramCnt1
		# Loop Over the Additional Features
		for (j in 1:ncol(zDistStain1)) {
							
			# Stain 1 Dist Features	
			stainIdx <- j
			# 05.01.2014
			# featureName <- stain1Names[j]
			featureName <- inputHeader[-1][j]
			resp <- repDistStain1[idxWells, stainIdx]
			concResp <- data.frame(x = cmpdConc, y = resp)
			concResp <- na.omit(concResp)
			yLim <- range(repDistStain1[, stainIdx], na.rm = TRUE)
		
			# If <2 concentration has a non-NA response, Efficacy = 0, potency = highest conc
			if (length(which(!is.na(resp))) > 1) {	
				fit  <- f_fitDR(concResp)
				scoredParam <- f_evalFit(fit, rSquaredMin, concResp)
				}else{ # Fit Response
				fit  <- NULL
				scoredParam <- c(0, -log10(max(cmpdConc)/10^6))
				}
				
				
				# Make Plot
				f_plotRightAxisSimple(cmpdConc, repCntStain1[idxWells, ], cntFit1, range(repCntStain1), "Cell Count", "gray")
				f_plotLeftAxisSimple(cmpdConc, resp, fit, yLim, scoredParam, featureName, "black")
			
				# Assign Outputs
				featNames[1:2 + j*2] <- paste(featureName, c("Efficacy","Potency"),sep = " - ")
				scoredDR[cmpdNum, 1:2 + j*2] <- scoredParam
				}
		# Final Plot is Blank
		plot(5,axes=FALSE, type="n", main = "", xlab="", ylab="", cex.main=1.5, bty="n")	
			}
			
		return(list(featNames = featNames, scoredDR = scoredDR, cmpds = cmpds, cmpdAliases = cmpdAliases))
		}

# f_plotRightAxis <- function() {
			# plot(rep(cmpdConc, numRep), as.vector(repCntStain1[idxWells, ]), col = "gray", pch = 18,
				 # log = "x", bty = "n", xlab = "Conc [uM]", ylab = "", axes = FALSE, 
				 # ylim = range(repCntStain1))
			# if (length(coef(cntFit1)) == 1) {
				# abline(h = coef(cntFit1), lwd = 2, col = "gray")
				# }else{
				# plot(cntFit1, add = TRUE, col = "gray", lwd = 2, type = "n")
				# }
			# axis(4, pretty(range(repCntStain1), 7), col='gray', col.ticks='gray',col.axis="gray")
			# mtext("Cell Count",4, line =2.5,col="gray")
			# }

f_plotRightAxis <- function(cmpdConc, resp, fit, yLim, featureName, ptsCol) {
			plot(rep(cmpdConc, ncol(resp)), as.vector(resp), col = ptsCol, pch = 18,
				 log = "x", bty = "n", xlab = "Conc [uM]", ylab = "", axes = FALSE, 
				 ylim = yLim)
			if (length(coef(fit)) == 1) {
				abline(h = coef(fit), lwd = 2, col = "gray")
				}else{
				plot(fit, add = TRUE, col = "gray", lwd = 2, type = "n")
				}
			axis(4, pretty(yLim, 7), col='gray', col.ticks='gray',col.axis="gray")
			mtext(featureName, 4, line =2.5,col="gray")
			}
f_plotRightAxisSimple <- function(cmpdConc, resp, fit, yLim, featureName, ptsCol) {
			plot(cmpdConc, resp, col = ptsCol, pch = 18,
				 log = "x", bty = "n", xlab = "Conc [uM]", ylab = "", axes = FALSE, 
				 ylim = yLim)
			if (length(coef(fit)) == 1) {
				abline(h = coef(fit), lwd = 2, col = "gray")
				}else{
				plot(fit, add = TRUE, col = "gray", lwd = 2, type = "n")
				}
			axis(4, pretty(yLim, 7), col='gray', col.ticks='gray',col.axis="gray")
			mtext(featureName, 4, line =2.5,col="gray")
			}		
# f_plotLeftAxis <- function(resp, fit, yLim, scoredParam, featureName, ptsCol) {
			# par(new=TRUE)
			# plot(rep(cmpdConc, numRep), as.vector(resp), pch=19, cex = 1.25, bty="n", axes=FALSE,
					 # col="black", xlab="Conc [uM]", lwd=2, log="x", ylab="", ylim=yLim,
					 # col.main = ifelse(scoredParam[1] == 0, "black", "black"),
					 # main = paste(featureName,"\n Efficacy = ", round(scoredParam[1],1),
					              # "; Potency = ", round(scoredParam[2],1), sep=""))
			# points(rep(cmpdConc, numRep), as.vector(resp), pch=19, cex=1.25, col=ptsCol)
			# if(is.null(fit$convergence) == TRUE){
				 # if(length(coef(fit)) > 1){
					 # plot(fit, add=TRUE, col="red", lwd=2, type = "n")
					 # }else{
					 # abline(h = coef(fit), lwd = 2, col = "red")
					 # }
					# }
			# axis(2,pretty(yLim, 7), col="black",col.ticks="black",col.axis="black")
			# axis(1,axTicks(1), col='black',col.ticks='black',col.axis='black')
			# }
		
f_plotLeftAxis <- function(cmpdConc, resp, fit, yLim, scoredParam, featureName, ptsCol) {
			par(new=TRUE)
			plot(rep(cmpdConc, ncol(resp)), as.vector(resp), pch=19, cex = 1.25, bty="n", axes=FALSE,
					 col="black", xlab="Conc [uM]", lwd=2, log="x", ylab="", ylim=yLim,
					 col.main = ifelse(scoredParam[1] == 0, "black", "black"),
					 main = paste(featureName,"\n Efficacy = ", round(scoredParam[1],1),
					              "; Potency = ", round(scoredParam[2],1), sep=""))
			points(rep(cmpdConc, ncol(resp)), as.vector(resp), pch=19, cex=1.25, col=ptsCol)
			if(is.null(fit$convergence) == TRUE){
				 if(length(coef(fit)) > 1){
					 plot(fit, add=TRUE, col="red", lwd=2, type = "n")
					 }else{
					 abline(h = coef(fit), lwd = 2, col = "red")
					 }
					}
			axis(2,pretty(yLim, 7), col="black",col.ticks="black",col.axis="black")
			axis(1,axTicks(1), col='black',col.ticks='black',col.axis='black')
			}
					
f_plotLeftAxisSimple <- function(cmpdConc, resp, fit, yLim, scoredParam, featureName, ptsCol) {
			par(new=TRUE)
			plot(cmpdConc, resp, pch=19, cex = 1.25, bty="n", axes=FALSE,
					 col="black", xlab="Conc [uM]", lwd=2, log="x", ylab="", ylim=yLim,
					 col.main = ifelse(scoredParam[1] == 0, "black", "black"),
					 main = paste(featureName,"\n Efficacy = ", round(scoredParam[1],1),
					              "; Potency = ", round(scoredParam[2],1), sep=""))
			points(cmpdConc, resp, pch=19, cex=1.25, col=ptsCol)
			if(is.null(fit$convergence) == TRUE){
				 if(length(coef(fit)) > 1){
					 plot(fit, add=TRUE, col="red", lwd=2, type = "n")
					 }else{
					 abline(h = coef(fit), lwd = 2, col = "red")
					 }
					}
			axis(2,pretty(yLim, 7), col="black",col.ticks="black",col.axis="black")
			axis(1,axTicks(1), col='black',col.ticks='black',col.axis='black')
			}