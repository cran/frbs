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
#' This is the internal function that implements the hybrid neural fuzzy inference 
#' system (HyFIS). Users do not need to call it directly,
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
#' @param range.data a matrix(2 x n) containing the range of the normalized data, where n is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param num.labels a matrix(1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param max.iter the maximal number of iterations.
#' @param range.data.ori a matrix containing the ranges of the original data. 
#' @param step.size step size of the gradient descent method. 
#' @seealso \code{\link{HyFIS.update}}, \code{\link{frbs.learn}}, and \code{\link{predict}}.
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components.
#' @references 
#' J. Kim and N. Kasabov, "HyFIS: adaptive neuro-fuzzy inference systems and their application to nonlinear dynamical systems", 
#' Neural Networks 12, 1301 - 1319 (1999).
#'
# @export
HyFIS <- function(range.data, data.train, num.labels, max.iter = 100, range.data.ori, step.size = 0.01) {
	
	type.mf = 3
	
	mod <- WM(range.data, data.train, num.labels, type.mf)
	
	data.test <- as.matrix(data.train[, 1 : (ncol(data.train) - 1)])
	
	range.input <- mod$range.input
	range.output <- mod$range.output
	num.varinput <- mod$num.varinput
	num.fvalinput <- mod$num.fvalinput
	names.varinput <- mod$names.varinput
	varout.mf <- mod$varout.mf
	names.varoutput <- mod$names.varoutput
	rule <- mod$rule
	rule.temp <- mod$rule
	varinp.mf <- mod$varinp.mf
	degree.rule <- mod$degree.rule
	mod$method.type <- "HYFIS"
	
	var.mf <- cbind(varinp.mf, varout.mf)
	var.mf.old <- cbind(varinp.mf, varout.mf)
	
	
	for (iter in 1 : max.iter){
		for (i in 1 : nrow(data.test)){
			dt.i <- matrix(data.test[i, ], nrow = 1)
			dt.train.i <- matrix(data.train[i, ], nrow = 1)					
						
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
			
			MSE <- mean(residuals^2)
					
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
	mod <- list(range.input = range.input, range.output = range.output, num.varinput = num.varinput, num.fvalinput = num.fvalinput,
	     names.varinput = names.varinput, varout.mf = varout.mf, names.varoutput = names.varoutput, rule = rule, varinp.mf = varinp.mf, 
		 range.data.ori = range.data.ori, degree.rule = degree.rule, method.type = "HYFIS")
	
	return (mod)
}  




