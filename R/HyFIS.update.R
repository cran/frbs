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
#' This function is called by \code{\link{HyFIS}} to update the parameters within the HyFIS method.
#'
#' @title HyFIS updating function
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param def matrix of defuzzification results. See \code{\link{defuzzifier}}.
#' @param rule fuzzy IF-THEN rules. See \code{\link{rulebase}}.
#' @param names.varoutput a list of names of the output variable. 
#' @param var.mf a matrix of parameters of the membership functions. Please see \code{\link{fuzzifier}}.
#' @param miu.rule a matrix of degree of rules which is a result of the \code{\link{inference}}.
#' @param num.labels a matrix(1 x n) whose elements represent the number of labels (or fuzzy terms), where n is the number of variables.
#' @param MF a matrix of parameters of the membership functions which is a result of the \code{\link{fuzzifier}}.
#' @param step.size a real number, the step size of the gradient descent. 
#' @param degree.rule a matrix of degrees of rules. See \code{\link{frbs-object}}.
#' @seealso \code{\link{HyFIS}}
# @return a matrix of parameters of membership of function
# @export
HyFIS.update <- function(data.train, def, rule, names.varoutput, var.mf, miu.rule, num.labels, MF, step.size= 0.001, degree.rule){


rule.temp <- rule

indx <- length(num.labels)
l.output <- (ncol(var.mf) - num.labels[1, indx])

data <- data.train

## get varinp.mf
varinp.mf <- matrix(var.mf[, 1 : l.output], nrow = 5)

## get varout.mf
varout.mf <- matrix(var.mf[, (l.output + 1) : ncol(var.mf)], nrow = 5)

##Check not zero on miu.rule
chck <- which(miu.rule[1, ] > 0.001)


l.chck <- length(chck)
			
## initialize
temp <- 0
temp.d <- 0
temp1 <- 0
temp2 <- 0
			
#if there are some no zero element on miu rule
if (l.chck != 0) {
	indx <- matrix(nrow = l.chck, ncol = 3)
	temp.indx <- matrix(nrow = 1, ncol = 3)

	## along number of not zero on miu.rule, check and save string related on names.varoutput and its value
	for (ii in 1 : l.chck){
		#get the names of output value on each rule
		strg <- c(rule.temp[chck[ii], ncol(rule.temp)])
		aaa <- which(strg == names.varoutput)
		if (length(aaa) != 0) {
			indx[ii, 1] <- aaa
			indx[ii, 2] <- miu.rule[1, chck[ii]]
			indx[ii, 3] <- degree.rule[chck[ii]]
		}
	}			
	
## check duplication on indx, choose with has a max of degree
if (nrow(indx) >= 2){
	for (jj in 1 : (nrow(indx) - 1)){
		temp.indx[1, ] <- indx[jj, ]
		for (kk in (jj + 1) : nrow(indx)){
			ck <- length(which(temp.indx[1, 1] == indx[kk, 1]))
			if (ck != 0){
				if ((temp.indx[1, 2] * temp.indx[1, 3]) >= (indx[kk, 2] * indx[kk, 3]))
					indx[kk, ] <- NA
				else 
					indx[jj, ] <- NA
				}
			}
		}				
	}

# erase NA value on indx
indx <- na.omit(indx)

	#####################
	## Update center(mean) and width(variance) on output variable (layer 5)
	####################

	eta.o <- step.size
	eta.i <- step.size
	delta <- (data.train[1, ncol(data.train)] - def)

	####jacobian 
	ss <- seq(from = 1, to = num.labels[1, ncol(num.labels)])

	indx.comp <- cbind(ss, 0)
	for (i in 1:nrow(indx)){
		indx.comp[indx[i, 1], 2] <- (indx[i, 2] * indx[i, 3])
	}

	mean.t <- t(t(varout.mf[2, ]))
	var.t <- t(t(varout.mf[3, ]))

	###update output variable
	delta.k <- 0
	for (i in 1:ncol(varout.mf)){
		
		cc <- ceiling(i / num.labels[1, 1] )
		
		if (sum(var.t * indx.comp[, 2]) != Inf) {
			delta.k4 <-  delta * varout.mf[3, i] * (varout.mf[2, i] * sum(indx.comp[, 2] * var.t) - sum(indx.comp[, 2] * mean.t * var.t))* (1 / (sum(indx.comp[, 2] * var.t))^2)
			temp <- delta * indx.comp[i, 2] * (varout.mf[2, i] * sum(indx.comp[, 2] * var.t) - sum(indx.comp[, 2] * mean.t * var.t))* (1 / (sum(indx.comp[, 2] * var.t))^2)
			varout.mf[3, i] <- varout.mf[3, i] + eta.o * temp
			varout.mf[2, i] <- varout.mf[2, i] + eta.o * delta * (varout.mf[2, i] * indx.comp[i, 2]) / (sum(var.t * indx.comp[, 2]))
		}
		else {
			delta.k4 <-  0
		}
			
		delta.k <- delta.k + delta.k4
		if (varout.mf[2, i] < 0)
		varout.mf[2, i] <- 0
		if (varout.mf[3, i] < 0)
		varout.mf[3, i] <- 0.001
	
		if (varout.mf[2, i] > 1)
		varout.mf[2, i] <- 1
		if (varout.mf[3, i] > 1)
		varout.mf[3, i] <- 1

	}
	######################
	## Update center(mean) and width(variance) on input variable (layer 3)
	######################

	num.varinput <- ncol(data.train) - 1
	for (j in 0 : (num.varinput - 1)){
		

		start.v <- j * num.labels[1, j + 1] + 1
		end.v <- start.v + num.labels[1, j + 1] - 1
		term.min <- which.min(MF[1, start.v : end.v]) 
		
		## get data correspondent with input variable
		indxx <- j * num.labels[1, j + 1] + term.min
		
		################
			a <- MF[1, indxx] * 2 *(data.train[1, j  + 1] -  varinp.mf[2, indxx])/(varinp.mf[3, indxx])^2		
			varinp.mf[2, indxx] <- varinp.mf[2, indxx] + eta.i * delta.k * a * delta	
			b <- MF[1, indxx] * 2 *(data.train[1, j + 1] -  varinp.mf[2, indxx])^2/(varinp.mf[3, indxx])^3	
			varinp.mf[3, indxx] <- varinp.mf[3, indxx] + eta.i * delta.k * b * delta

			if (varinp.mf[2, indxx] < 0)
			varinp.mf[2, indxx] <- 0
			if (varinp.mf[3, indxx] < 0)
			varinp.mf[3, indxx] <- 0.001
		
			if (varinp.mf[2, indxx] > 1)
			varinp.mf[2, indxx] <- 1
			if (varinp.mf[3, indxx] > 1)
			varinp.mf[3, indxx] <- 1
	}	
}
var.mf <- list(varinp.mf = varinp.mf, varout.mf = varout.mf)
return(var.mf)

}