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
#' This is the internal function that implements the adaptive-network-based 
#' fuzzy inference system (ANFIS). It is used to solve regression tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by Jyh-Shing and Roger Jang. It uses the Takagi Sugeno Kang model
#' on the consequent part of the fuzzy IF-THEN rules.  
#' The ANFIS architecture consists of two processes, the forward and the backward stage.
#' The forward stage has five layers as follows:
#' \itemize{
#' \item Layer 1: The fuzzification process which transforms crisp values into 
#' fuzzy terms using the Gaussian function as the shape of the membership function.
#' \item Layer 2: The inference stage using the t-norm operator (the AND operator).
#' \item Layer 3: Calculating the ratio of the strengths of the rules.
#' \item Layer 4: Calculating the consequent parameters.
#' \item Layer 5: Calculating the overall output as the sum of all incoming signals.
#' }
#' The backward stage is a process of parameter learning. In this step, the least squares method 
#' is used in order to obtain 
#' the parameters, which are coefficients of linear equations on the consequent part, and 
#' mean and variance on the antecedent part.
#' 
#' @title ANFIS model building 
#'
#' @param data.train a matrix(m x n) of normalized data for the training process, 
#'        where m is the number of instances and 
#'         n is the number of variables; the last column is the output variable.
#'        Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix(1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param max.iter the maximal number of iterations.
#' @param step.size a real number between 0 and 1 representing the step size of 
#' the gradient descent. 
#' @param type.tnorm the type of t-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.snorm the type of s-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.implication.func a value representing type of implication function. 
#'        For more detail, please have a look at \code{\link{WM}}. 
#' @seealso \code{\link{ANFIS.update}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return list of the model data. please have a look at \code{\link{frbs.learn}} 
# for looking complete components.
#' @references
#' Jyh-Shing and Roger Jang, "ANFIS: adaptive-network-based fuzzy inference system",
#' IEEE Transactions on Systems, Man, and Cybernetics, vol. 23, no. 3, pp. 665 - 685 (1993).
#'
#' Jyh-Shing Roger Jang, et.al., "Neuro-fuzzy and soft computing: 
#' a computational approach to learning and machine intelligence", Prentice-Hall, Inc (1997).
# @export
ANFIS <- function(data.train, num.labels, max.iter = 10, step.size = 0.01, type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH") {
	
	## fixed value for ANFIS
	type.mf <- "GAUSSIAN"
	range.data <- matrix(nrow = 2, ncol = ncol(data.train))
	range.data[1, ] <- 0
	range.data[2, ] <- 1
	
	## initialize rule and membership function by WM
	mod <- WM(data.train, num.labels, type.mf, type.tnorm, type.implication.func)
	
	## make data test from data training
	data.test <- as.matrix(data.train[, 1 : (ncol(data.train) - 1)])
	
	## get parameters	
	range.input <- range.data[, -ncol(range.data), drop = FALSE]
	range.output <- range.data[, ncol(range.data), drop = FALSE]
	
	num.varinput <- ncol(data.train) - 1
	num.labels.input <- num.labels[, -ncol(num.labels), drop = FALSE]
	names.fvalinput <- colnames(mod$varinp.mf)
	
	rule <- mod$rule

	varinp.mf <- mod$varinp.mf
	rule.data.num <- mod$rule.data.num

	## check a duplication fuzzy linguistics on antecedent part on each rules (should change with duplicated) 
	temp.rule <- matrix(nrow = 1, ncol = (ncol(rule) - 1))
	for (i in 1 : (nrow(rule) - 1)){
		temp.rule[1, ] <- rule[i, 1 : (ncol(rule) - 1)]
		
		for (j in (i + 1) : nrow(rule)) {
			chk <- which(temp.rule[1, ] == rule[j, 1 : (ncol(rule) - 1)])
			if (length(chk) == length(temp.rule)){
				rule[i, ] <- NA
				rule.data.num[i, ] <- NA
			}
		}
	}
	
	## delete NA value
	rule <- na.omit(rule)
	rule.data.num <- na.omit(rule.data.num)
	
	#number of rule
	n.rowrule <- nrow(rule)
	
	#generate func.tsk by random number as consequent part
	num.ele <- n.rowrule * (num.varinput + 1)
	rand <- runif(num.ele, min = 0, max = 1)
	func.tsk<-matrix(rand, nrow = n.rowrule, ncol= num.varinput + 1, byrow=T)
	
	## set into TSK model
	type.model <- "TSK"

	## 	iteration for updating parameters
	## updating uses online learning (it's iterated one by one)
	for (iter in 1 : max.iter){
		for (i in 1 : nrow(data.test)){
			
			## get ith data test
			data <- data.test[i, ,drop = FALSE]
			
			## get ith data training for fitting process
			dt.train.i <- data.train[i, ,drop = FALSE]
		
			### I. Rule Based Module
			rule <- rulebase(type.model, rule, func.tsk)

			### II. Fuzzification Module
			MF <- fuzzifier(data, num.varinput, num.labels.input, varinp.mf)
				
			### III. Inference Module
			miu.rule <- inference(MF, rule, names.fvalinput, type.tnorm, type.snorm)

			## update parameters
			param.new <- ANFIS.update(dt.train.i, range.input, range.output, rule.data.num, miu.rule, func.tsk, varinp.mf, step.size)
			
			## obtain new parameters: antecedent and consequent parts
			func.tsk <- param.new$func.tsk.new
			varinp.mf <- param.new$varinp.mf
		
		}

	}
	
	## change rule format into TSK model 
	length.rule <- length(rule)
	temp <- matrix(rule[[1]], nrow = 1)
	for (i in 2 : length.rule){
		temp.1 <- matrix(rule[[i]], nrow = 1)
		temp <- rbind(temp, temp.1)
	}	
	rule <- temp
	
	num.labels <- num.labels[, -ncol(num.labels), drop = FALSE]
	## collect into mod list
	mod <- list(num.labels = num.labels, rule = rule, rule.data.num = rule.data.num, 
	              varinp.mf = varinp.mf, func.tsk = func.tsk, type.tnorm = type.tnorm, 
				  type.snorm = type.snorm, type.defuz = NULL, type.model = "TSK", type.mf = "GAUSSIAN", type.implication.func = type.implication.func)
	
	return (mod)
}  


#' This is the internal function that implements the hybrid neural fuzzy inference 
#' system (HyFIS). It is used to solve regression tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}
#' 
#' This method was proposed by J. Kim and N. Kasabov. There are two phases in 
#' this method for learning, namely the knowledge acquisition module and the 
#' structure and parameter learning.
#' The knowledge acquition module uses the techniques of Wang and Mendel. 
#' The learning of structure and parameters is a supervised learning method using 
#' gradient descent-based learning algorithms. 
#' This function generates a model which consists of a rule database and parameters 
#' of the membership functions. The rules of HyFIS use the Mamdani model on the antecedent and 
#' consequent parts. Futhermore, 
#' HyFIS uses a Gaussian membership function. So, there are two kinds of 
#' parameters that are optimized, mean and variance of the Gaussian function.  
#' 
#' @title HyFIS model building 
#' 
#' @param data.train a matrix (m x n) of normalized data for the training process, 
#'         where m is the number of instances and 
#'         n is the number of variables; the last column is the output variable.
#'        Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix (1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param max.iter the maximal number of iterations.
#' @param step.size step size of the gradient descent method. 
#' @param type.tnorm the type of t-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.snorm the type of s-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.defuz the type of aggregation function. For more detail, please have a look at \code{\link{defuzzifier}}
#' @param type.implication.func a value representing type of implication function. For more detail, please have a look at \code{\link{WM}}
#' @seealso \code{\link{HyFIS.update}}, \code{\link{frbs.learn}}, and \code{\link{predict}}.
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components.
#' @references 
#' J. Kim and N. Kasabov, "HyFIS: adaptive neuro-fuzzy inference systems and 
#' their application to nonlinear dynamical systems", 
#' Neural Networks, vol. 12, no. 9, pp. 1301 - 1319 (1999).
#'
# @export
HyFIS <- function(data.train, num.labels, max.iter = 10, step.size = 0.01, type.tnorm = "MIN", 
               type.snorm = "MAX", type.defuz = "COG", type.implication.func = "ZADEH") {
	
	type.mf = "GAUSSIAN"
	range.data <- matrix(nrow = 2, ncol = ncol(data.train))
	range.data[1, ] <- 0
	range.data[2, ] <- 1
	
	mod <- WM(data.train, num.labels, type.mf, type.tnorm, type.implication.func)
	
	data.test <- as.matrix(data.train[, 1 : (ncol(data.train) - 1)])
	varout.mf <- mod$varout.mf
	names.varoutput <- colnames(varout.mf)
	rule <- mod$rule
	rule.temp <- mod$rule
	varinp.mf <- mod$varinp.mf
	degree.rule <- mod$degree.rule
	mod$method.type <- "HYFIS"
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$type.defuz <- type.defuz
	mod$type.model <- "MAMDANI"
	mod$func.tsk <- NULL
	
	var.mf <- cbind(varinp.mf, varout.mf)
	var.mf.old <- cbind(varinp.mf, varout.mf)
	
	
	for (iter in 1 : max.iter){
		for (i in 1 : nrow(data.test)){
			dt.i <- data.test[i, ,drop = FALSE]
			dt.train.i <- data.train[i, ,drop = FALSE]
						
			res <- frbs.eng(mod, dt.i)
			def <- res$predicted.val
			
		    miu.rule <- res$miu.rule
			MF <- res$MF

		    ## measure error 
			y.pred <- def
			y.real <- data.train[i, ncol(data.train)]
	
			residuals <- (y.real - y.pred)
			RMSE <- sqrt(mean(residuals^2))		
			error <- RMSE 
					
		    ## stoping criteria by RMSE
			if (error <= 0.001){
				break
			}
			else {
				new.var.mf <- HyFIS.update(dt.train.i, def, rule, names.varoutput, var.mf, miu.rule, num.labels, MF, step.size, degree.rule)		
			}
			
			## update parameters
			mod$varout.mf <- new.var.mf$varout.mf 
			mod$varinp.mf <- new.var.mf$varinp.mf 
			var.mf <- cbind(mod$varinp.mf, mod$varout.mf)
		}	
	}
	varinp.mf <- mod$varinp.mf
	varout.mf <- mod$varout.mf
	
	rule <- rule.temp	
	mod <- list(num.labels = num.labels, varout.mf = varout.mf, rule = rule, varinp.mf = varinp.mf, func.tsk = NULL, 
		 degree.rule = degree.rule, rule.data.num = mod$rule.data.num, method.type = "HYFIS", 
		 type.tnorm = type.tnorm, type.snorm = type.snorm, type.defuz = type.defuz, type.model = "MAMDANI",
		 type.mf = "GAUSSIAN", type.implication.func = type.implication.func)
	
	return (mod)
}  
