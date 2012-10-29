#############################################################################
#
#  This file is a part of the R package "frbs".
#
#  Author: Lala Septem Riza
#  Co-author: Christoph Bergmeir
#  Supervisors: Francisco Herrera Triguero and Jose Manuel Benitez
#  Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' Defuzzification is a transformation that extracts the crisp values from the fuzzy terms. 
#' 
#' In this function, there exist two kinds of models which are based on the Mamdani model and the Takagi Sugeno Kang model on the consequent parts. 
#' For the Mamdani model there are five methods for defuzzifying a fuzzy term A of a universe of discourse Z. 
#' They are as follows:
#' \enumerate{
#' \item weighted average method
#' \item first of maxima
#' \item last of maxima
#' \item mean of maxima
#' \item modified COG
#' }
#' 
#' @title Defuzzifier to transform from fuzzy terms to crisp values
#'
#' @param data a matrix(m x n) of data, where m is the number of instances and 
#' n is the number of variables.
#' @param rule a list or matrix of fuzzy IF-THEN rules, as discussed in \code{\link{rulebase}}.
#' @param range.output a matrix(2 x n) containing the range of the output data. 
#' @param names.varoutput a list for giving names to the fuzzy terms. See \code{\link{rulebase}}.
#' @param varout.mf a matrix for constructing the membership function of the output variable. See \code{\link{fuzzifier}}.
#' @param miu.rule the results of the inference module. See \code{\link{inference}}.
#' @param type.defuz the type of defuzzification to be used, where 1 means weighted average method, and 2, 3, 4 and 5 mean first, last, mean maxima and modified COG, respectively.
#' @param type.model the type of the model that will be used in the simulation. Here, 1 or 2 means we use Mamdani or Takagi Sugeno Kang (which includes the possibility for a constant value), respectively.
#' @param func.tsk a matrix used to build the linear equation for the consequent part if we are using Takagi Sugeno Kang. See also \code{\link{rulebase}}.
#' @seealso \code{\link{fuzzifier}}, \code{\link{rulebase}}, and \code{\link{inference}}
#' @return A matrix of crisp values
# @export
defuzzifier <- function(data, rule, range.output, names.varoutput = NULL, varout.mf = NULL, miu.rule, type.defuz = 1, type.model = 1, func.tsk = NULL){

## copy rule
rule.temp <- matrix(unlist(rule), nrow = length(rule), byrow= TRUE)

## Inisialitation
def <- matrix(0, nrow=nrow(data), ncol = 1)
def.temp <- matrix(0, nrow=nrow(data), ncol = 1)

## Mamdani type
if (type.model == 1){

	## check names.varoutput
	if (is.null(names.varoutput)){
		stop("please define the names of output variable")
	}

	## check parameters of membership functions on output variables
	if (is.null(varout.mf)){
		stop("please define the parameters of membership functions on the output variable")
	}

	## Inisialitation
	cum <- matrix(0, nrow=nrow(data), ncol = ncol(varout.mf)) 
	div <- matrix(0, nrow=nrow(data), ncol = ncol(varout.mf))

	for (k in 1 : nrow(data)){
		##Check not zero on miu.rule
		chck <- which(miu.rule[k, ] != 0)
		l.chck <- length(chck)
		cum <- matrix(0, nrow=nrow(data), ncol = l.chck) 
		div <- matrix(0, nrow=nrow(data), ncol = l.chck)
		
		## initialize
		temp <- 0
		temp.d <- 0
		temp1 <- 0
		temp2 <- 0
		
		#if there are some no zero element on miu rule
		if (l.chck != 0) {
			indx <- matrix(nrow = l.chck, ncol = 2)
			temp.indx <- matrix(nrow = 1, ncol = 2)
			## along number of not zero on miu.rule, check and save string related on names.varoutput and its value
			for (ii in 1 : l.chck){
				#get the names of output value on each rule
				strg <- c(rule.temp[chck[ii], ncol(rule.temp)])
				aaa <- which(strg == names.varoutput)
				if (length(aaa) != 0) {
				indx[ii, 1] <- aaa
				indx[ii, 2] <- miu.rule[k, chck[ii]]
				}
			}
			
			## check duplication on indx, choose with has a max of degree
			if (nrow(indx) >= 2){
				for (jj in 1 : (nrow(indx) - 1)){
					temp.indx[1, ] <- indx[jj, ]
					for (kk in (jj + 1) : nrow(indx)){
						ck <- length(which(temp.indx[1, 1] == indx[kk, 1]))
						if (ck != 0){
							if (temp.indx[1, 2] >= indx[kk, 2])
								indx[kk, ] <- NA
							else 
								indx[jj, ] <- NA
						}
					}
				}				
			}
			
			# erase NA value on indx
			indx <- na.omit(indx)

			#defuzzification procedure for  Centroid
			if (type.defuz == 1) {								
				for (i in 1 : nrow(indx)){										
					#calculate modified centroid
					# update for gaussian
					if (varout.mf[1, indx[i, 1]] == 5 || varout.mf[1, indx[i, 1]] == 6 || varout.mf[1, indx[i, 1]] == 7) {
						## center point
						av.point <- varout.mf[2, indx[i, 1]]
						
						## indx is fired miu.rule/rule
						temp1 <- indx[i, 2] * av.point 
						temp2 <- indx[i, 2] 
												
						temp <- temp + temp1
						temp.d <- temp.d + temp2
					}
					
					else {						
						av.point <- varout.mf[3, indx[i, 1]] 	
						
						# check first which one greater between indx[i, 2] with the value of MF (for centroid)
						temp1 <- indx[i, 2] * av.point 
						
						temp2 <- indx[i, 2]
						temp <- temp + temp1
						temp.d <- temp.d + temp2
					}
					cum[k, i] <- temp
					div[k, i] <- temp.d   					
				}
							
				if (sum(div[k, ]) == 0){
					def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
				}
				else{
					def.temp[k, 1] <- sum(cum[k, ]) / sum(div[k, ])
					if (def.temp[k, 1] <= min(range.output, na.rm=TRUE)){
						def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
					}
					else if (def.temp[k, 1] >= max(range.output, na.rm=TRUE)){
						def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
					}
					else {
						def[k, 1] <- def.temp[k, 1]
					}
				}
			}
			## procedure for type.defuz == 2 (fisrt of Maxima)
			else if (type.defuz == 2){
				max.temp <- max(indx[, 2], na.rm = TRUE)
				max.indx <- which(indx[, 2] == max.temp)			
				
				aa <- varout.mf[2, indx[max.indx[1], 1]]
				bb <- varout.mf[3, indx[max.indx[1], 1]]
				cc <- varout.mf[4, indx[max.indx[1], 1]]
				dd <- varout.mf[5, indx[max.indx[1], 1]]
									
				# check shape of membership function
				if (varout.mf[1, indx[max.indx[1], 1]] == 1){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
											
					for (j in 1:length(seqq)){
							if (seqq[j] < aa){
								temp.miu <- 1
							}
							else if (seqq[j] >= aa & seqq[j] < bb){
								temp.miu <- (seqq[j] - aa) / (bb - aa)
							}
							else if (seqq[j] >= bb & seqq[j] < cc) {
								temp.miu <- (seqq[j] - cc) / (bb - cc)
							}
							else 
								temp.miu <- 1
							
							if (temp.miu >= indx[max.indx[1], 2]) {
								def[k, 1] <- seqq[j]
								break
							}
					}					
			
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 2){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					
					for (j in 1:length(seqq)){
						if (seqq[j] < bb){
							temp.miu <- 1
						}
						else if (seqq[j] >= bb & seqq[j] < cc){
							temp.miu <- (seqq[j] - cc) / (bb - cc)
						}
						else if (seqq[j] > cc) {
							temp.miu <- 1
						}
														
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
					}										
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 3){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					for (j in 1:length(seqq)){
						if (seqq[j] < aa){
							temp.miu <- 0
						}
						else if (seqq[j] >= aa & seqq[j] < bb){
							temp.miu <- (seqq[j] - aa) / (bb - aa)
						}
						else if (seqq[j] > cc) {
							temp.miu <- 1
						}
														
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
					}
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 4) {
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					
					for (j in 1:length(seqq)){
						if (seqq[j] < aa){
							temp.miu <- 0
						}
						else if (seqq[j] >= aa & seqq[j] < bb){
							temp.miu <- (seqq[j] - aa) / (bb - aa)
						}
						else if (seqq[j] >= bb & seqq[j] < cc) {
							temp.miu <- 1
						}
						else if (seqq[j] >= cc & seqq[j] < dd) {
							temp.miu <- (seqq[j] - dd) / (cc - dd)
						}
						else {
							temp.miu <- 0
						}
						
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
					}
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 5) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 100)
					
					for (j in 1:length(seqq)){
						
						temp.miu <- exp(- 0.5 * (seqq[j] - aa) ^ 2 / bb ^ 2)
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 6) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					for (j in 1:length(seqq)){
						
						temp.miu <- 1 / (1 + exp(- aa * (seqq[j] - bb)))
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 7) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					for (j in 1:length(seqq)){
						
						temp.miu <- 1 / (1 + abs((seqq[j] - cc)/aa) ^ (2 * bb))
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
			}	
			
			## procedure for type.defuz == 3 (last of maxima)
			else if (type.defuz == 3) {
				max.temp <- max(indx[, 2], na.rm = TRUE)
				max.indx <- which(indx[, 2] == max.temp)			
				
				aa <- varout.mf[2, indx[max.indx[1], 1]]
				bb <- varout.mf[3, indx[max.indx[1], 1]]
				cc <- varout.mf[4, indx[max.indx[1], 1]]
				dd <- varout.mf[5, indx[max.indx[1], 1]]
									
				# check shape of membership function
				if (varout.mf[1, indx[max.indx[1], 1]] == 1){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
											
					for (j in 1:length(seqq)){
							if (seqq[j] < aa){
								temp.miu <- 1
							}
							else if (seqq[j] >= aa & seqq[j] < bb){
								temp.miu <- (seqq[j] - aa) / (bb - aa)
							}
							else if (seqq[j] >= bb & seqq[j] < cc) {
								temp.miu <- (seqq[j] - cc) / (bb - cc)
							}
							else 
								temp.miu <- 1
							
							if (temp.miu >= indx[max.indx[1], 2]) {
								def[k, 1] <- seqq[j]
							}
					}					
			
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 2){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					
					for (j in 1:length(seqq)){
						if (seqq[j] < bb){
							temp.miu <- 1
						}
						else if (seqq[j] >= bb & seqq[j] < cc){
							temp.miu <- (seqq[j] - cc) / (bb - cc)
						}
						else if (seqq[j] > cc) {
							temp.miu <- 1
						}
														
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
					}										
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 3){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					for (j in 1:length(seqq)){
						if (seqq[j] < aa){
							temp.miu <- 0
						}
						else if (seqq[j] >= aa & seqq[j] < bb){
							temp.miu <- (seqq[j] - aa) / (bb - aa)
						}
						else if (seqq[j] > cc) {
							temp.miu <- 1
						}
														
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
					}
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 4) {
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					
					for (j in 1:length(seqq)){
						if (seqq[j] < aa){
							temp.miu <- 0
						}
						else if (seqq[j] >= aa & seqq[j] < bb){
							temp.miu <- (seqq[j] - aa) / (bb - aa)
						}
						else if (seqq[j] >= bb & seqq[j] < cc) {
							temp.miu <- 1
						}
						else if (seqq[j] >= cc & seqq[j] < dd) {
							temp.miu <- (seqq[j] - dd) / (cc - dd)
						}
						else {
							temp.miu <- 0
						}
						
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
					}
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 5) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					
					for (j in 1:length(seqq)){
						
						temp.miu <- exp(- 0.5 * (seqq[j] - aa) ^ 2 / bb ^ 2)
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 6) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					for (j in 1:length(seqq)){
						
						temp.miu <- 1 / (1 + exp(- aa * (seqq[j] - bb)))
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 7) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					for (j in 1:length(seqq)){
						
						temp.miu <- 1 / (1 + abs((seqq[j] - cc)/aa) ^ (2 * bb))
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}	
			}
			
			## procedure for type.defuz == 4 (mean of maxima)
			else if (type.defuz == 4) {
				max.temp <- max(indx[, 2], na.rm = TRUE)
				max.indx <- which(indx[, 2] == max.temp)			
							
				def[k, 1] <- 0.5 * (max(varout.mf[2:5, indx[max.indx[1], 1]], na.rm = TRUE) + (varout.mf[2, indx[max.indx[1], 1]]))
			}
			
			else if (type.defuz == 5) {
				for (i in 1 : nrow(indx)){										
					#calculate modified centroid
					# update for gaussian
					if (varout.mf[1, indx[i, 1]] == 5 || varout.mf[1, indx[i, 1]] == 6 || varout.mf[1, indx[i, 1]] == 7) {
						## center point
						
						av.point <- varout.mf[2, indx[i, 1]]
						
						## indx is fired miu.rule/rule
						temp1 <- indx[i, 2] * av.point * varout.mf[3, indx[i, 1]]
						temp2 <- indx[i, 2] * varout.mf[3, indx[i, 1]]

						temp <- temp + temp1
						temp.d <- temp.d + temp2
					}
					
					else {				
						average <- (varout.mf[2, indx[i, 1]] + max(varout.mf[2 : 5, indx[i, 1]], na.rm = TRUE)) / 2
						av.point <- varout.mf[3, indx[i, 1]]
						# check first which one greater between indx[i, 2] with the value of MF (for centroid)
						temp1 <- indx[i, 2] * av.point #* abs((average - max(varout.mf[2 : 5, indx[i, 1]], na.rm = TRUE))) / 2
						temp2 <- indx[i, 2] # * (abs((average - max(varout.mf[2 : 5, indx[i, 1]], na.rm = TRUE))) / 2)
						temp <- temp + temp1
						temp.d <- temp.d + temp2
					}
					cum[k, i] <- temp
					div[k, i] <- temp.d   					
				}
				
				if (sum(div[k, ]) == 0){
					def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
				}
				else{
					def.temp[k, 1] <- sum(cum[k, ]) / sum(div[k, ])
					if (def.temp[k, 1] <= min(range.output, na.rm=TRUE)){
						def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
					}
					else if (def.temp[k, 1] >= max(range.output, na.rm=TRUE)){
						def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
					}
					else {
						def[k, 1] <- def.temp[k, 1]
					}
				}
			}
			
		}
		else {
			def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
		}
	}	

}

#TSK type
else if (type.model == 2){
	
	## check the linear equation on consequent part
	if (is.null(func.tsk)){
		stop("please define the linear equation as the consequent part on fuzzy IF-THEN rules")
	}
		
	for (k in 1 : nrow(data)){
		data.k <-(data[k, ])
		data.m <- as.matrix(data.k)
		
		if (ncol(func.tsk) > 1){
			func.tsk.var <- func.tsk[, 1: ncol(func.tsk) - 1]
			func.tsk.cont <- t(t(func.tsk[, ncol(func.tsk)]))
		
			ff <- func.tsk.var %*% data.m + func.tsk.cont
		}
		else if (ncol(func.tsk) == 1){
			ff <- func.tsk 
		}
		
		miu.rule.t <- as.matrix(miu.rule[k, ])
		cum <- sum(miu.rule.t * ff)	
		div <- sum(miu.rule.t)
		
		def[k, 1] <- cum / div
		if (div == 0)
			def[k, 1] <- 0
		else
		{
			def[k, 1] <- cum / div
			if (def[k, 1] > max(range.output, na.rm=TRUE))
				def[k, 1] <- max(range.output, na.rm=TRUE)
			else if (def[k, 1] < min(range.output, na.rm=TRUE))
				def[k, 1] <- min(range.output, na.rm=TRUE)
		}
	}
}

return(def)
}

