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
#' This function is one of the main internal functions of the package. 
#' It determines the values within the prediction phase.  
#'
#' This function involves four different processing steps on fuzzy rule-based systems. 
#' Firstly, the rulebase (see \code{\link{rulebase}}) validates 
#' the consistency of the fuzzy IF-THEN rules form. Then, the fuzzification 
#' (see \code{\link{fuzzifier}}) transforms crisp values 
#' into fuzzy terms. Next, the inference calculates the degree of rule strengths using 
#' the t-norm and the s-norm. 
#' Finally, the defuzzification process calculates the results of the model using the Mamdani 
#' or the Takagi Sugeno Kang model.  
#'
#' @title The prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix(m x n) of data for the prediction process, 
#' where m is the number of instances and n is the number of input variables.
#' @seealso \code{\link{fuzzifier}}, \code{\link{rulebase}}, \code{\link{inference}} 
#' and \code{\link{defuzzifier}}.
#' @return A list with the following items:
#' \item{rule}{the fuzzy IF-THEN rules}
#' \item{varinp.mf}{a matrix to generate the shapes of the membership functions for 
#' the input variables}
#' \item{MF}{a matrix of the degrees of the membership functions}
#' \item{miu.rule}{a matrix of the degrees of the rules}
#' \item{func.tsk}{a matrix of the Takagi Sugeno Kang model for the consequent part of 
#' the fuzzy IF-THEN rules}
#' \item{predicted.val}{a matrix of the predicted values}
#' 
# @export
frbs.eng <- function(object, newdata){
	mod <- object

	data.test <- newdata
	## get all of parameters
	range.input <- mod$range.input
	range.output <- mod$range.output
	num.varinput <- mod$num.varinput	
	num.fvalinput <- mod$num.fvalinput
	names.varinput <- mod$names.varinput
	varout.mf <- mod$varout.mf
	names.varoutput <- mod$names.varoutput
	rule <- mod$rule
	varinp.mf <- mod$varinp.mf	
	names.variable <- c(names.varinput, names.varoutput)
	
		
	if (is.null(mod$method.type) == FALSE && mod$method.type == "ANFIS"){
		type.defuz <- 1
		type.tnorm <- 1
		type.snorm <- 1
		type.model <- 2
		func.tsk <- mod$func.tsk
	}
	else if (is.null(mod$method.type) == FALSE && mod$method.type == "FS.HGD"){
		type.model <- 2
		type.tnorm <- 1
		type.snorm <- 1
		type.defuz <- 1
		func.tsk <- mod$func.tsk
	}	
	else if (is.null(mod$method.type) == FALSE && mod$method.type == "FIR.DM"){
		type.model <- 2
		type.defuz <- 1
		type.tnorm <- 1
		type.snorm <- 1
		func.tsk <- mod$func.tsk
	}
	else if (is.null(mod$method.type) == FALSE && mod$method.type == "HYFIS"){
		type.model = 1
		func.tsk = NULL
		type.defuz = 5
		type.tnorm = 1
		type.snorm = 1
	}
	else {
		type.defuz <- mod$type.defuz
		type.tnorm <- mod$type.tnorm
		type.snorm <- mod$type.snorm	
		type.model <- mod$type.model
		func.tsk <- mod$func.tsk
	}
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
	ncol.MF <- ncol(MF)
	names.var <- names.variable[1 : ncol.MF]
	colnames(MF) <- c(names.var)
	
	miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)
	
	if(is.null(mod$method.type) == FALSE && mod$method.type == "FS.HGD"){
		miu.rule <- miu.rule ^ mod$alpha.heuristic
	}
	
	if(is.null(mod$method.type) == FALSE && mod$method.type == "HYFIS"){
		degree.rule <- mod$degree.rule
		for (i in 1 : nrow(miu.rule)){
			miu.rule[i, ] <- miu.rule[i, ] * t(degree.rule)
		}
	}
	
	###################
	### IV. Defuzzification Module
	### In this function, we calculate and convert fuzzy value back into crisp value. 
	###################
	def <- defuzzifier(data, rule, range.output, names.varoutput, varout.mf, miu.rule, type.defuz, type.model, func.tsk)
	
	res <- list(rule = rule, varinp.mf = varinp.mf, MF = MF, miu.rule = miu.rule, func.tsk = func.tsk, predicted.val = def)
	return(res)
	
}

#' This function is the internal function of the FRBCS method to compute the predicted values.  
#'
#' @title FRBCS: prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix(m x n) of data for the prediction process, 
#' where m is the number of instances and n is the number of input variables.
#' @return A matrix of predicted values.
# @export
FRBCS.eng <- function(object, newdata){
	mod <- object
	data.test <- newdata
	num.varinput <- mod$num.varinput
	num.fvalinput <- mod$num.fvalinput
	names.varinput <- mod$names.varinput
	num.fvaloutput <- mod$num.fvaloutput
	rule <- mod$rule
	varinp.mf <- mod$varinp.mf
	type.model <- 1
	func.tsk <- mod$class
	type.tnorm <- 1
	type.snorm <- 1
	grade.certainty <- mod$grade.cert
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