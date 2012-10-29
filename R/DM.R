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
#' This is the internal function that implements the fuzzy inference rules by descent method (DM). 
#' Users do not need to call it directly, but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by Hiroyoshi Nomura, Isao Hayashi, and Noboru Wakami. DM uses simplified fuzzy 
#' reasoning where the consequent part is a real number (a particular case within the Takagi Sugeno Kang model),
#' while the membership function on the 
#' antecedent part is expressed by an isosceles triangle. So, in the learning phase, DM updates three parameters 
#' which are center and width of the triangular and a real number on the consequent part using a descent method.
#'
#' @title DM model building 
#'
#' @param range.data a matrix(2 x n) containing the range of the normalized data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param data.train a matrix(m x n) of data for training, where m is the number of instances and 
#' n is the number of variables. The last column is the output variable.
#' @param num.labels a matrix(1 x n) whose elements represent the number of labels (fuzzy terms),
#' where n is the number of variables.
#' @param max.iter the maximal number of iterations.
#' @param step.size the step size of the descent method, between 0 and 1.
#' @seealso \code{\link{DM.update}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking complete components. 
#' @references
#' H. Nomura, I. Hayashi and N. Wakami, "A learning method of fuzzy inference rules by descent method," Proc. FUZZ-IEEE'92
#' pp. 203 - 210 (1992).
# @export

DM <- function(range.data, data.train, num.labels, max.iter, step.size){
	type.mf = 1  
	type.model <- 2
	type.defuz <- 1
	type.tnorm <- 1
	type.snorm <- 1	

	## generate initial model using WM
	mod <- WM(range.data, data.train, num.labels, type.mf)
	
	## get data test from data training
	data.test <- as.matrix(data.train[, 1 : (ncol(data.train) - 1)])

	## get values of model	
	range.input <- mod$range.input
	range.output <- mod$range.output
	num.varinput <- mod$num.varinput
	num.fvalinput <- mod$num.fvalinput
	names.varinput <- mod$names.varinput
	rule <- mod$rule
	rule.data.num <- mod$rule.data.num		
	varinp.mf <- mod$varinp.mf

	## check a duplication on rule 
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
	rule <- na.omit(rule)
	rule.data.num <- na.omit(rule.data.num)
	
	## for TSK
	rule <- rule[, 1 : (ncol(rule) - 1)]
	mod$rule <- rule

	n.rowrule <- nrow(rule)	
	gal.iter <- matrix(nrow = max.iter)
	
	#generate func.tsk by random number
	num.ele <- n.rowrule 
	rand <- runif(num.ele, min = 0, max = 1)
	func.tsk <- matrix(rand, nrow = n.rowrule, ncol= 1, byrow = TRUE)
	
	mod$type.model <- 2
	mod$func.tsk <- func.tsk
	mod$type.defuz <- 1
	mod$type.tnorm <- 1
	mod$type.snorm <- 1
	
	## calculate frbs.eng in order to obtain predicted values (def)
	init.res <- frbs.eng(mod, data.test)		
	init.def <- init.res$predicted.val
	
	init.miu.rule <- init.res$miu.rule
	
	## update initial parameters on consequent part
	for (h in 1 : nrow(data.test)) {
		init.gal <- init.def[h] - data.train[h, ncol(data.train)]
		for (i in 1 : nrow(func.tsk)){
			if (sum(init.miu.rule[h, ]) != 0)
				func.tsk[i, 1] <- func.tsk[i, 1] - step.size * init.miu.rule[h, i]/sum(init.miu.rule[h, ]) * init.gal 
		}
	}
	
	## calcualte error
	init.gal.l <- init.def - data.train[, ncol(data.train)]
	
	## update linear eq.
	mod$func.tsk <- func.tsk

	galat <- matrix()
	error <- matrix()
	
	## update processes
	for (iter in 1 : max.iter){
		for (i in 1 : nrow(data.test)){
		
			data <- matrix(data.test[i, ], nrow = 1)
			dt.train.i <- matrix(data.train[i, ], nrow = 1)			
			res <- frbs.eng(mod, data)			
			def <- res$predicted.val
			
			error[i] <- def - dt.train.i[1, ncol(dt.train.i)]
			
			miu.rule <- res$miu.rule
			MF <- res$MF
			rule <- res$rule
			func.tsk <- res$func.tsk
			varinp.mf <- res$varinp.mf
			
			
			## procedure for getting new parameters
			param.new <- DM.update(dt.train.i, rule.data.num, miu.rule, func.tsk, varinp.mf, step.size, def)

			func.tsk <- param.new$func.tsk.new
			varinp.mf <- param.new$varinp.mf.n
			
			## update with new params
			mod$varinp.mf <- varinp.mf
			mod$func.tsk <- func.tsk
					
		}

		galat[iter] <- sqrt(mean(error^2))
		if (galat[iter] < 0.001) 
			break
	}
	
	rule <- mod$rule
	
	mod <- list(range.input = range.input, range.output = range.output, num.varinput = num.varinput, num.fvalinput = num.fvalinput, names.varinput = names.varinput, 
	            rule = rule, varinp.mf = varinp.mf, func.tsk = func.tsk)
	
	return (mod)
}