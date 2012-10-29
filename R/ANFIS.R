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
#' This is the internal function that implements the adaptive-network-based fuzzy inference system (ANFIS). Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by Jyh-Shing and Roger Jang. It uses the Takagi Sugeno Kang model on the consequent part of the fuzzy IF-THEN rules.  
#' The ANFIS architecture consists of two processes, the forward and the backward stage. The forward stage has five layers as follows:
#' \itemize{
#' \item Layer 1: The fuzzification process which transforms crisp values into fuzzy terms using the Gaussian function as the shape of the membership function.
#' \item Layer 2: The inference stage using the t-norm operator (the AND operator).
#' \item Layer 3: Calculating the ratio of the strengths of the rules.
#' \item Layer 4: Calculating the consequent parameters.
#' \item Layer 5: Calculating the overall output as the sum of all incoming signals.
#' }
#' The backward stage is a process of parameter learning. In this step, the least squares method is used in order to obtain 
#' the parameters, which are coefficients of linear equations on the consequent part, and mean and variance on the antecedent part.
#' 
#' @title ANFIS model building 
#'
#' @param range.data a matrix(2 x n) containing the range of the normalized data, where n is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param num.labels a matrix(1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param max.iter the maximal number of iterations.
#' @param range.data.ori a matrix containing the ranges of the original data. 
#' @param step.size a real number between 0 and 1 representing the step size of the gradient descent. 
#' @seealso \code{\link{ANFIS.update}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return list of the model data. please have a look at \code{\link{frbs.learn}} for looking complete components.
#' @references
#' Andri Riid and Ennu Rustern, "Interpretability versus adaptability in fuzzy systems," 
#' Proceedings of the Estonian Academy of Sciences, Engineering, 6, 2, pp. 76 - 95 (2000).
#'
#' Babuska R. and Verbruggen, H. "Neuro-fuzzy methods for nonlinear system identification," 
#' Annual Reviews in Control, 27 I, pp. 73-85 (2003).
#' 
#' Jyh-Shing and Roger Jang, "ANFIS: adaptive-network-based fuzzy inference system",
#' IEEE Transactions on Systems, Man, and Cybernetics, Vol. 23, No. 3 (1993).
# @export
ANFIS <- function(range.data, data.train, num.labels, max.iter = 100, range.data.ori, step.size = 0.01) {
	
	## fixed value for ANFIS
	type.defuz <- 1
	type.mf <- 3
	type.tnorm <- 1
	type.snorm <- 1
	
	## initialize rule and membership function by WM
	mod <- WM(range.data, data.train, num.labels, type.mf)
	
	## make data test from data training
	data.test <- as.matrix(data.train[, 1 : (ncol(data.train) - 1)])
	
	## get parameters	
	range.input <- mod$range.input
	range.output <- mod$range.output
	num.varinput <- mod$num.varinput
	num.fvalinput <- mod$num.fvalinput
	names.varinput <- mod$names.varinput
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
	type.model <- 2

	## 	iteration for updating parameters
	## updating uses online learning (it's iterated one by one)
	for (iter in 1 : max.iter){
		for (i in 1 : nrow(data.test)){
			
			## get ith data test
			data <- matrix(data.test[i, ], nrow = 1)
			
			## get ith data training for fitting process
			dt.train.i <- matrix(data.train[i, ], nrow = 1)
		
			### I. Rule Based Module
			rule <- rulebase(type.model, rule, func.tsk)

			### II. Fuzzification Module
			MF <- fuzzifier(data, num.varinput, num.fvalinput, varinp.mf)
				
			### III. Inference Module
			miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)

			## update parameters
			param.new <- ANFIS.update(dt.train.i, range.input, range.output, rule.data.num, miu.rule, func.tsk, varinp.mf, step.size)
			
			## obtain new parameters: antecedent and consequent parts
			func.tsk <- param.new$func.tsk.new
			varinp.mf <- param.new$varinp.mf
		
		}

	}

	## collect into mod list
	mod <- list(range.input = range.input, range.output = range.output, num.varinput = num.varinput, num.fvalinput = num.fvalinput, names.varinput = names.varinput,
            	rule = rule, varinp.mf = varinp.mf, func.tsk = func.tsk, range.data.ori = range.data.ori)
	
	return (mod)
}  




