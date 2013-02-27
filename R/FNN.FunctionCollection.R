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
#' The role of this function is to update parameters in the ANFIS method. 
#' This function is called by the main function of the ANFIS method, \code{\link{ANFIS}}.
#'
#' @title ANFIS updating function
#'
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param range.input the range of the input variables, as a matrix(2 x n).
#' @param range.output the range of the output variable, as a matrix(2 x n).
#' @param rule.data.num a matrix containing the rule base in integer form. 
#' @param miu.rule a matrix with the degrees of rules. See \code{\link{inference}}.
#' @param func.tsk a matrix of parameters of the function on the consequent part using the Takagi Sugeno Kang model.
#' @param varinp.mf a matrix of parameters of membership functions of the input variables.
#' @param step.size a real number between 0 and 1 representing the step size of 
#' the gradient descent. 
#' @seealso \code{\link{ANFIS}}, \code{\link{rulebase}}, \code{\link{inference}}, \code{\link{fuzzifier}} and \code{\link{defuzzifier}}
# @return A list containing the following components:
# \item{func.tsk.new}{a matrix of parameters of the consequence part, if using Takagi Sugeno Kang model}
# \item{varinp.mf}{a matrix to generate the shapes of the membership functions on the input variables}
# @export
ANFIS.update <- function(data.train, range.input, range.output, rule.data.num, miu.rule, func.tsk, varinp.mf, step.size = 0.01){

## for now, we use single output
num.varoutput <- 1 
data <- data.train
alpha <- step.size

## set data into matrix with number of column 1
data.m <-matrix(data[1, 1 : (ncol(data.train) - num.varoutput)], ncol = 1)

## split linear equation on consequent part	
func.tsk.var <- func.tsk[, 1: ncol(func.tsk) - 1]
func.tsk.cont <- func.tsk[, ncol(func.tsk), drop = FALSE]

## calculate linear eq.	
ff <- func.tsk.var %*% data.m + func.tsk.cont

## set into matrix
miu.rule.t <- as.matrix(miu.rule)

## get nominator on consequent part
cum <- sum(ff * t(miu.rule.t))

## get denominator on consequent part	
div <- sum(miu.rule.t)

## set if divide by zero and calculate predicted value "def"
if (div == 0)
	def <- 0
else {
	def <- cum / div
	if (def > max(range.output, na.rm=TRUE))
		def <- max(range.output, na.rm=TRUE)
	else if (def < min(range.input, na.rm=TRUE))
		def <- min(range.input, na.rm=TRUE)
	}

## find max value on miu.rule
indxx <- which.max(miu.rule.t)

## calculate galat/error
gal <- (def - data.train[1, ncol(data.train)])

## initialize 
delta.out <- 0

## update on consequent part, TSK model
	for (m in 1 : nrow(func.tsk.var)){
		for (n in 1 : ncol(func.tsk.var)){
			sum.miu <- sum(miu.rule.t, na.rm = TRUE)
			if (sum.miu != 0)
				## update consequent part
				func.tsk.var[m, n] <- func.tsk.var[m, n] - alpha * gal * (miu.rule[m]/sum.miu) *  data.train[1, n]	
			else
				func.tsk.var[m, n] <- func.tsk.var[m, n] - alpha * gal * (miu.rule[m]) *  data.train[1, n]	
		}	
	}
	
	for (mm in 1 : nrow(func.tsk.cont)){
		sum.miu <- sum(miu.rule.t, na.rm = TRUE)
		if (sum.miu != 0)
			## update constant on linear eq of consequent part
			func.tsk.cont[mm, 1] <- func.tsk.cont[mm, 1] - alpha * gal * (miu.rule[mm]/sum.miu)
		else
			func.tsk.cont[mm, 1] <- func.tsk.cont[mm, 1] - alpha * gal * miu.rule[mm]
 	}

## collect new linear eq. 
func.tsk.new <- cbind(func.tsk.var, func.tsk.cont)

####################################################
#=== update input variable using gradient descent
###################################################

## get number of variable
num.varinp <- ncol(data.train) - num.varoutput

## update input variables
for (ii in 1 : ncol(miu.rule.t)) {
		
	## split linear eq.
	func.tsk.var.indxx <- func.tsk[ii, 1: ncol(func.tsk) - 1]
	func.tsk.cont.indxx <- func.tsk[ii, ncol(func.tsk)]
	
	## calculate predicted value
	ff.indxx <- func.tsk.var.indxx %*% data.m + func.tsk.cont.indxx
	tau.miu <- ff.indxx * miu.rule.t[1, ii]
	if (div == 0 )
		div <- 0.0001
	func.val <- tau.miu / div
	
	
	num.varinput <- ncol(data.train) - num.varoutput
	
	## update parameters in antecedent part
	for (j in 1 : num.varinput){
		
		varinp.mf[2, rule.data.num[ii, j]] <- varinp.mf[2, rule.data.num[ii, j]] - step.size * gal * (func.val - data.train[1, ncol(data.train)] ) * miu.rule.t[1, ii] / div * 
						(data.train[1, j] -  varinp.mf[2, rule.data.num[ii, j]])/(varinp.mf[3, rule.data.num[ii, j]])^2
		varinp.mf[3, rule.data.num[ii, j]] <- varinp.mf[3, rule.data.num[ii, j]] - step.size * gal * (func.val - data.train[1, ncol(data.train)] ) * miu.rule.t[1, ii] / div * 
						(data.train[1, j] -  varinp.mf[2, rule.data.num[ii, j]]) ^ 2/(varinp.mf[3, rule.data.num[ii, j]])^3
		
		if (varinp.mf[2, rule.data.num[ii, j]] < 0)
			varinp.mf[2, rule.data.num[ii, j]] <- 0
		if (varinp.mf[3, rule.data.num[ii, j]] < 0)
			varinp.mf[3, rule.data.num[ii, j]] <- 0.001
		
		if (varinp.mf[2, rule.data.num[ii, j]] > 1)
			varinp.mf[2, rule.data.num[ii, j]] <- 1
		if (varinp.mf[3, rule.data.num[ii, j]] > 1)
			varinp.mf[3, rule.data.num[ii, j]] <- 1
	}	
	
}	

## collect new parameters
param.new <- list(func.tsk.new = func.tsk.new, varinp.mf = varinp.mf)

return(param.new)
}


#' This function is called by \code{\link{HyFIS}} to update the parameters within 
#' the HyFIS method.
#'
#' @title HyFIS updating function
#' @param data.train a matrix(m x n) of data for the training process, 
#' where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param def matrix of defuzzification results. See \code{\link{defuzzifier}}.
#' @param rule fuzzy IF-THEN rules. See \code{\link{rulebase}}.
#' @param names.varoutput a list of names of the output variable. 
#' @param var.mf a matrix of parameters of the membership functions. 
#' Please see \code{\link{fuzzifier}}.
#' @param miu.rule a matrix of degree of rules which is a result of the \code{\link{inference}}.
#' @param num.labels a matrix(1 x n) whose elements represent 
#' the number of labels (or fuzzy terms), where n is the number of variables.
#' @param MF a matrix of parameters of the membership functions 
#' which is a result of the \code{\link{fuzzifier}}.
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
varinp.mf <- var.mf[, 1 : l.output, drop = FALSE]

## get varout.mf
varout.mf <- var.mf[, (l.output + 1) : ncol(var.mf), drop = FALSE]

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



