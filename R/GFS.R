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
#' This is the internal function that implements the genetic fuzzy systems (GFS) 
#' model. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}
#' 
#' This function is adopted from Jindrich Liska and Stephen S. Melsheimer's paper. 
#' GFS implements a genetic algorithm in order to determine the structure of the fuzzy IF-THEN rules 
#' and the membership function parameters. Genetic algorithms (GA) are a stochastic optimization 
#' technique that mimics natural selection. We conduct selection, uniform crossover, mutation, 
#' and the creep operator for the GA. One important thing when using a GA is to choose the 
#' representation of the solutions within the individuals. 
#' We split every individual to become two matrices which are a matrix of membership 
#' function parameters for all variables and a matrix of the structure of the fuzzy IF-THEN rules. Initially, 
#' one individual is generated using the techniques of Wang and Mendel, and the rest is generated at random. 
#' So, this method is always better or at least the same as the Wang and Mendel method. 
#' The fitness function is computed as the root mean square error for each individual. 
#' 
#' @title GFS model building 
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param range.data a matrix(2 x n) containing the range of the normalized data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param num.labels a matrix(1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param popu.size the size of the population which is generated in each generation.
#' @param persen_cross a real number between 0 and 1 representing the probability of crossover.
#' @param persen_mutant  a real number between 0 and 1 representing the probability of mutation.
#' @param max.iter the maximal number of iterations.
#' @param classification a boolean representing whether it is a classification problem or not.
#' @param range.data.ori a matrix(2 x n) containing the range of the original data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
# @return a list of Genetic Fuzzy System model. please have a look \code{\link{predict}}.
#' @references 
#' Davis L., "Adaptive operator probabilities in genetic algorithm," ICGA'89 (1989).
#'
#' Holland J. H., "Adaptation in natural and artificial systems," University of Michigan Press, Ann Arbor (1975).
#'
#' Jindrich Liska and Stephen S. Melsheimer, "Complete design of fuzzy logic systems using genetic algorithm," 
#' Fuzzy Systems, IEEE World Congress on Computational Intelligence (1994).
#'
#' Jindrich Liska and Stephen S. Melsheimer, "Design of fuzzy logic systems for nonlinear model identification," American Control Conference (1994).
# @export
GFS <- function(data.train, range.data, num.labels, popu.size, persen_cross, persen_mutant, max.iter, classification = FALSE, range.data.ori){
	type.mf = 3
	
	### Wang Mendel for initial individu
	mod <- WM(range.data, data.train, num.labels, type.mf)
    
	num.varinput <- mod$num.varinput
	names.varinput <- mod$names.varinput
	varout.mf <- mod$varout.mf
	names.varoutput <- mod$names.varoutput
	rule <- mod$rule
	varinp.mf <- mod$varinp.mf
	rule.data.num <- mod$rule.data.num
	mod$type.mf <- 3
	mod$type.model <- 1
	mod$func.tsk <- NULL
	mod$type.defuz <- 1
	mod$type.tnorm <- 1
	mod$type.snorm <- 1
	mod$range.data.ori <- range.data.ori
	
	n.labels <- num.labels[1,1]
	names.variable <- c(names.varinput, names.varoutput)
	num.var <- num.varinput + 1
	
	###Build representation of Individi on GA
	############### 1. Representation of MF of Input Variabel and Output Variabel
	
	#combine variance input and output variables
	width.inp <- t(varinp.mf[3, ])
	width.out <- t(varout.mf[3, ])

	width <- cbind(width.inp, width.out)
	
	#combine center input and output variables
	center.inp <- t(varinp.mf[2, ])
	center.out <- t(varout.mf[2, ])
	center <- cbind(center.inp, center.out)	
	
	## it is matrix one dimension: width | center of input and output variables
	## 
	MF.param <- matrix(nrow = 1, ncol = (2 * ncol(width)))
	j <- 1
	k <- 1
	for (i in 1 : (2 * ncol(width))){
		if (i %% 2 == 1){
			MF.param[1, i] <- width[1, j]	
			j <- j + 1
		}
		else {
			MF.param[1, i] <- center[1, k]
			k <- k + 1
		}	
	}
			
	MF.temp <- matrix(nrow = (popu.size - 1), ncol = (2 * ncol(width)))
	
	### representation parameters of MF 
	## MF.param is from wang mendel
	MF.param.popu <- MF.param
	MF.param.init <- MF.param
		
	## generate other individu by random
	for (i in 1 : (popu.size - 1)){
		leap <- 0
		for (j in 1 : ncol(MF.param)){
			if (j %% 2 == 1) {
				eps <- runif(1, min = 0.005, max = 0.3)
				MF.temp[i, j] <- eps	
			}
			else {
				loc.low <- leap * n.labels * 2 + 2
				if (j == loc.low ) {
					min.eps <- 0
					max.eps <- 1
					eps <- runif(1, min = min.eps, max = max.eps)
					MF.temp[i, j] <- eps
					leap <- leap + 1
				}
				else {
					min.eps <- MF.temp[i, j - 2]
					max.eps <- 1
					eps <- runif(1, min = min.eps, max = max.eps)
					MF.temp[i, j] <- eps
				}
			}
		}
	}
	
	## population of MF parameters
	MF.param.popu <- rbind(MF.param.popu, MF.temp)
	
	############### 2. Representation of rule 
	init.popu <- matrix(t(rule.data.num[1, ]), nrow = 1)
	for (i in 2 : nrow(rule.data.num)){
		temp.init <- t(rule.data.num[i, ])
		init.popu <- cbind(init.popu, temp.init)
	}
	
	temp <- matrix(nrow = nrow(rule.data.num), ncol =  (num.varinput + 1)) 
	
	## initial population from wang mendel
	rand.popu <- init.popu
	rule.popu.init <- init.popu
	
	## generate population of rule by random 
	for (h in 1 : (popu.size - 1)) {
		#generate random for input variable
		for (i in 1 : (num.varinput + 1)) {
			min.runif <- (i - 1) * num.labels + 1
			temp[, i] <- as.integer(runif((nrow(rule.data.num)), min = min.runif, max = (min.runif + num.labels - 1)))
		}

		rand.popu.temp <- matrix(t(temp[1, ]), nrow = 1)
		for (i in 2 : nrow(temp)){
			temp.rand <- t(temp[i, ])
			rand.popu.temp <- cbind(rand.popu.temp, temp.rand)
		}
		
		rand.popu <- rbind(rand.popu, rand.popu.temp)
	}
	
	## initial rule population
	rule.popu <- rand.popu
	
	## getting data test from data train
	data.test <- data.train[, 1 : num.varinput]

	stop <- 0
	best.fit <- matrix()
	
	while (stop <= max.iter){
		#####################
		## convert population into rule database
		####################
		res.comp <- matrix(nrow = popu.size + 1, ncol = nrow(data.test))
		for (i in 1 : popu.size){
			
			## convert rule on population into matrix of rule 
			mod$rule <- convert.rule(matrix(rule.popu[i, ], nrow=1), num.varinput, n.labels, names.variable)
			
			## convert var.mf on population into matrix of parameters of MF
			var.mf <- convert.MF(matrix(MF.param.popu[i, ], nrow=1), varinp.mf, varout.mf, n.labels, num.varinput)
			
			## set into mod
			mod$varinp.mf <- var.mf$varinp.mf
			mod$varout.mf <- var.mf$varout.mf 		
						
			## testing using new rule and var.mf 
			res <- frbs.eng(mod, data.test)
			
			res.temp <- matrix(res$predicted.val, nrow = 1)
			
			## collect predicted value
			res.comp <- rbind(res.comp, res.temp)						
		}
		
		res.comp <- na.omit(res.comp)
	
		## check if the purpose is classification
		if (classification == TRUE) {
			res.comp <- round(res.comp)			
		}
		
		real.val <- matrix(data.train[, ncol(data.train)], nrow = 1)
		
		## calculate RMSE as fitness function
		RMSE <- matrix()	
		SMAPE <- matrix()
		for (i in 1 : nrow(res.comp)){
			residuals <- (real.val - matrix(res.comp[i, ], nrow = 1))
			MSE <- mean(residuals^2)
			RMSE[i] <- sqrt(mean(residuals^2))	
		}
		
		fitness <- t(RMSE)
		
		###the best fitness
		best.fit[stop] <- min(fitness)
		
		## get index which has min of fitness
		indx.best <- which.min(fitness)
		
		## save the best individu
		MF.elit.ind <- matrix(MF.param.popu[indx.best, ], nrow = 1)
		rule.elit.ind <- matrix(rule.popu[indx.best, ], nrow = 1)
				
		popu.elit <- list(MF.elit.ind = MF.elit.ind, rule.elit.ind = rule.elit.ind)
		
		## collect population of MF and rule 	
		popu.all <- list(MF.param.popu = MF.param.popu, rule.popu = rule.popu)
		
		### update population using selection, crossover, mutation, and creep operator
		new.popu <-  popubaru(popu.all, popu.elit, persen_cross, persen_mutant, n.labels, num.var)
		
		MF.param.popu <- new.popu$MF.new.popu
		rule.popu <- new.popu$rule.new.popu
		
		MF.param.popu[2, ] <- MF.param.init
		rule.popu[2, ] <- rule.popu.init
				
		stop <- stop + 1
	}
	
	#### The best individu
	
	mod$rule <- convert.rule(rule.elit.ind, num.varinput, n.labels, names.variable)
	var.mf <- convert.MF(MF.elit.ind, varinp.mf, varout.mf, n.labels, num.varinput)
	mod$varinp.mf <- var.mf$varinp.mf
	mod$varout.mf <- var.mf$varout.mf 
	
return(mod)
}
