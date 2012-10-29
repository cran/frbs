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
#' This function is the internal function of the frcs method to compute the predicted values.  
#'
#' @title frcs: prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix(m x n) of data for the prediction process, where m is the number of instances and 
#' n is the number of input variables.
#' @return A matrix of predicted values.
# @export
frcs.eng <- function(object, newdata){
	mod <- object
	data.test <- newdata
	num.varinput <- mod$num.varinput
	num.fvalinput <- mod$num.fvalinput
	names.varinput <- mod$names.varinput
	num.fvaloutput <- mod$num.fvaloutput
	rule <- mod$rule
	varinp.mf <- mod$varinp.mf
	type.model <- 1
	func.tsk <- mod$func.tsk
	type.tnorm <- 1
	type.snorm <- 1
	grade.certainty <- mod$grade.cert

	
	#if (type.model == 2){
	#	rule <- matrix(mod$rule, nrow=length(mod$rule), byrow=TRUE)
	#}
	
	rule.mat <- rule
	
	data <- data.test
	##################
	### I. Rule Based Module
	### In this function, Checking of the rule given by user will be done.
	### There are two kinds of model used which are Mamdani and TSK rule model.
	##################
	rule <- rulebase(type.model, rule, func.tsk)

	###################
	### II. Fuzzification Module
	### In this function, we convert crisp value into fuzzy value based on the data and parameter of membership function.
	### There are several membership function can be used such as triangular, trapezoid, gaussian and logistic/sigmoid.
	###################
	MF <- fuzzifier(data, num.varinput, num.fvalinput, varinp.mf)

	###################
	### III. Inference Module
	### In this function, we will calculate the confidence factor on antecedent for each rule. We use AND, OR, NOT operator. 
	################### 
	miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)
	
	alpha.class.all <- miu.rule
	
	for (i in 1 : nrow(miu.rule)){
		alpha.class.all[i, ] <- t(miu.rule[i, ] * grade.certainty[, 2])
	}
	indx.max <- matrix()
	for (i in 1 : nrow(miu.rule)){
		indx.max[i] <- which.max(alpha.class.all[i, ])
	}	
	result <- matrix()
	
	for (i in 1 : length(indx.max)){
		result[i] <- grade.certainty[indx.max[i], 1]
	}	
	res <- matrix(result)
	
	return(res)	
}