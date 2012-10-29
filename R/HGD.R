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
#' This is the internal function that implements the fuzzy system using heuristics 
#' and gradient descent method (HGD). Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by Ken Nozaki, H. Ishibuchi, and Hideo Tanaka. It 
#' uses fuzzy IF-THEN rules with nonfuzzy singletons (i.e. real numbers) in
#' the consequent parts. The techniques of Wang and Mendel are implemented to 
#' generate the antecedent part, while the initial consequent part of each 
#' rule is determined by the weighted mean value of the given training data. 
#' Then, the gradient descent method updates the value of the consequent part. 
#' Futhermore, the heuristic value given by the user affects the value of weight 
#' of each data. 
#' 
#' @title HGD model building 
#'
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param range.data a matrix(2 x n) containing the range of the normalized data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param num.labels a matrix(1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param max.iter maximal number of iterations.
#' @param step.size step size of the descent method. 
#' @param alpha.heuristic a positive real number which is the heuristic parameter.
#' @seealso \code{\link{frbs.learn}}, \code{\link{predict}}, and \code{\link{HGD.update}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components. 
#' @references
#' H. Ichihashi and T. Watanabe, "Learning control system by a simplified fuzzy reasoning model," Proc. IPMU'90 417 - 419 (1990). 
#'
#' H. Ishibuchi, K. Nozaki, H. Tanaka, Y. Hosaka and M. Matsuda, "Empirical study on learning in fuzzy systems by rice taste analysis,"
#' Fuzzy Set and Systems 64, 129 144 (1994).
# @export

HGD <- function(range.data, data.train, num.labels, max.iter, step.size, alpha.heuristic = 1){
	type.mf = 1
	type.model = 1
	func.tsk = NULL
	type.tnorm = 1
	type.snorm = 1
	
	mod <- WM(range.data, data.train, num.labels, type.mf)
	
	data.test <- as.matrix(data.train[, 1 : (ncol(data.train) - 1)])
	target.dt <- as.matrix(data.train[, ncol(data.train)], ncol = 1)
		
	range.input <- mod$range.input
	range.output <- mod$range.output
	num.varinput <- mod$num.varinput
	num.fvalinput <- mod$num.fvalinput
	names.varinput <- mod$names.varinput
	rule <- mod$rule
	varinp.mf <- mod$varinp.mf
		
	## check a duplication on rule (should change with duplicated)
	temp.rule <- matrix(nrow = 1, ncol = (ncol(rule) - 1))
	for (i in 1 : (nrow(rule) - 1)){
		temp.rule[1, ] <- rule[i, 1 : (ncol(rule) - 1)]
		
		for (j in (i + 1) : nrow(rule)) {
			chk <- which(temp.rule[1, ] == rule[j, 1 : (ncol(rule) - 1)])
			if (length(chk) == length(temp.rule)){
				rule[i, ] <- NA
			}
		}
	}
	rule <- na.omit(rule)
	## for TSK
	rule <- rule[, 1 : (ncol(rule) - 1)]
	mod$rule <- rule
	
	#number of rule
	n.rowrule <- nrow(rule)
	
	gal.iter <- matrix(nrow = max.iter)
	
	#generate func.tsk by heuristics
	#in order to initialize W
	data <- data.test
	rule <- rulebase(type.model, rule, func.tsk)
	MF <- fuzzifier(data, num.varinput, num.fvalinput, varinp.mf)
	
	miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)
	
	miu.rule <- miu.rule ^ alpha.heuristic
	
	func.tsk <- matrix(nrow = n.rowrule, ncol= 1, byrow=T)
	func.tsk <- (t(miu.rule) %*% target.dt) / colSums(miu.rule)
	
		
	type.model <- 2
	type.defuz <- 1
	
	ll <- length(mod)
	
	mod$type.model <- type.model
	mod$func.tsk <- func.tsk
	mod$type.defuz <- type.defuz
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$method.type <- "HGD"
	mod$alpha.heuristic <- alpha.heuristic
	
	galat <- matrix()
	error <- matrix()
		
	for (iter in 1 : max.iter){
		for (i in 1 : nrow(data.test)){
		
			dt.i <- matrix(data.test[i, ], nrow = 1)
			dt.train.i <- matrix(data.train[i, ], nrow = 1)
		
			data <- dt.i
	
			res <- frbs.eng(mod, data)
	
			def <- res$predicted.val
			
			error[i] <- def - dt.train.i[1, ncol(dt.train.i)]
			miu.rule <- res$miu.rule
			MF <- res$MF
			rule <- res$rule
			func.tsk <- res$func.tsk
			varinp.mf <- res$varinp.mf
			
			## update procedure
			param.new <- HGD.update(dt.train.i, miu.rule, func.tsk, varinp.mf, step.size, def)
			
			## get new parameters
			func.tsk <- unlist(param.new[[1]])			
			mod$func.tsk <- func.tsk		
		}
	galat[iter] <- sqrt(mean(error^2))
	}
	
	rule <- mod$rule
		
	mod <- list(range.input = range.input, range.output = range.output, num.varinput = num.varinput, num.fvalinput = num.fvalinput, names.varinput = names.varinput,
            	rule = rule, varinp.mf = varinp.mf, func.tsk = func.tsk, alpha.heuristic = alpha.heuristic)
		
	return (mod)
}