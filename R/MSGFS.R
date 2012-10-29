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
#' This is the internal function that implements the multi-stage genetic fuzzy systems (MSGFS) based on 
#' iterative rule learning approach. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by Herrera et al. MSGFS implements a genetic algorithm in order 
#' to determine the structure of the fuzzy IF-THEN rules 
#' and the membership function parameters. There are two general types of fuzzy IF-THEN rules, 
#' namely the descriptive and the approximative/free semantic approaches. 
#' A descriptive approach means that the linguistic labels represent a real-world semantic; 
#' the linguistic labels are uniformly defined for all rules. In contrast, in the approximative approach there isn't 
#' any associated linguistic label. This method is based on the latter one. We model a fuzzy IF-THEN rule on a chromosome
#' which consists of the parameter values of the membership function. So, every rule has its own membership function values. 
#' A population contains many such generated chromosomes, based on the iterative rule learning approach (IRL). 
#' IRL means that the chromosomes will be generated one by one, taking into account the fitness value and covering factor, 
#' until there are sufficient chromosomes in the population.
#' After having obtained the population, the genetic algorithm is started, using the genetic operators selection, mutation, and crossover. 
#'  
#' 
#' @title MSGFS model building 
#'
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param popu.size the size of the population which is generated in each generation.
#' @param persen_cross a real number between 0 and 1 representing the probability of crossover.
#' @param persen_mutant  a real number between 0 and 1 representing the probability of mutation.
#' @param max.iter the maximal number of iterations.
#' @param range.data.ori a matrix containing the ranges of the original data. 
#' @param epsilon a real number between 0 and 1 representing the boundary of covering factor. 
#' @seealso \code{\link{MSGFS.test}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return list of the model data. please have a look at \code{\link{frbs.learn}} for looking complete components.
#' @references
#' Antonio Gonzalez and Francisco Herrera, "Multi-stage genetic fuzzy systems based on the iterative rule learning approach,"
#' Mathware & Soft Computing 4, pp. 233 - 249 (1997).
#'
#' F. Herrera, M. Lozano, J.L. Verdegay, "Tuning fuzzy logic controllers by genetic algorithms," 
#' Internat. J. Approx. Reasoning 12, pp. 299 - 315 (1995).
#'
#' F. Herrera, M. Lozano, J.L. Verdegay, "Generating rules from examples using genetic algorithms," 
#' In: B. Bounchon, R. Yager, L. Zadeh (Eds.), Fuzzy Logic and Soft Computing, Word Scientific, pp. 11 - 20 (1995).
#'
#' F. Herrera, M. Lozano, and J.L. Verdegay, "A learning process for fuzzy control rules using genetic algorithms", 
#' Fuzzy Sets and Systems, 100, pp. 143 - 158 (1998).
#'
#' O. Cordon, M.J. del Jesus, F. Herrera, M. Lozano, "MOGUL: A methodology to obtain genetic fuzzy rule-based systems 
#' under the iterative rule learning approach," International Journal of Intelligent Systems, vol. 14, pp. 1123 - 1153 (1999).
# @export
MSGFS <- function(data.train, popu.size, persen_cross, persen_mutant, max.iter, range.data.ori, epsilon){
	## A genetic generation stage
	data.train.ori <- data.train
	num.var <- ncol(data.train)
	rule.gen <- matrix(genetic.gen(matrix(data.train[1, ], nrow =1), epsilon), nrow= 1)
	data.train <- matrix(data.train[-1, ], ncol = ncol(data.train))
	
	i <- nrow(data.train)
	while (i > 0){
		
		dt.i <- matrix(data.train[1, ], nrow =1)
		temp.rule <- matrix(genetic.gen(dt.i, epsilon), nrow = 1)
		data.train <- matrix(data.train[-1, ], ncol = ncol(data.train))
		i <- i - 1
		
		for (j in 1 : nrow(data.train)){
			if (j > nrow(data.train)){
				break
			}
			## check completeness and covering ==> ch.cover
			comp.degree <- ch.cover(matrix(data.train[j, ], nrow = 1), temp.rule)	

			if (comp.degree >= epsilon){
				data.train <- matrix(data.train[-j, ], ncol = ncol(data.train))
				i <- i - 1
			}
		}			
		rule.gen <- rbind(rule.gen, temp.rule)	
		
	}
	
	rule.gen.mean <- rule.gen[, seq(2, ncol(rule.gen), by = 3)]
	for (i in 1 : (nrow(rule.gen.mean) - 1)){	
		temp <- matrix(rule.gen.mean[i, ], nrow = 1)
	
		for (j in (i + 1) : nrow(rule.gen.mean)){
			margin.rule <- sum(abs(temp - matrix(rule.gen.mean[j, ], nrow = 1)), na.rm = TRUE) / ncol(rule.gen.mean)
			if (margin.rule <= 0.01) {
				rule.gen[j, ] <- NA
			}
		}
	}
	
	rule.gen <- na.omit(rule.gen)
	
	if (nrow(rule.gen) > popu.size){
		exceeded.rule.gen <- nrow(rule.gen) - popu.size
		indx <- sample(1 : nrow(rule.gen), size = exceeded.rule.gen) 
		rule.gen <- rule.gen[-indx, ]
	}
	else {
		added.rule.gen <- popu.size - nrow(rule.gen) 
		num.element <- added.rule.gen * ncol(rule.gen)
		added.rule <- matrix(runif(num.element, min = 0, max = 1), nrow = added.rule.gen, ncol = ncol(rule.gen))
		rule.gen <- rbind(rule.gen, added.rule)
	}
	
	##calculate individual fitness
	ind.fit <- calc.ind.fitness(data.train.ori, rule.gen)
	
	## combine rule.gen and its fitness
	rule.ind.fitness <- cbind(rule.gen, ind.fit)
	by.col <- ncol(rule.ind.fitness)
	
	rule.ind.fitness <- rule.ind.fitness[order(rule.ind.fitness[, by.col]), ]

	## calculate fitness 
	data.input <- data.train.ori[, -num.var]
	target.val <- matrix(data.train.ori[, ncol(data.train.ori)], nrow = nrow(data.train.ori))
	fit <- matrix()
	residual <- matrix()
	new.fitness <- 1000000000
	iter <- 1
	fit <- matrix(nrow = max.iter, ncol = 1)
	
	while (iter <= max.iter){
		
		## split in order to get rule.gen
		rule.gen <- matrix(rule.ind.fitness[, 1 : (ncol(rule.ind.fitness) - 1)], nrow = nrow(rule.ind.fitness))
		
		## calculate fitness
		for (i in 1 : nrow(data.input)){
			residual[i] <- calc.fitness(matrix(data.input[i, ], nrow = 1), target.val[i, 1], rule.gen)		
		}
		fitness <- sqrt(mean(residual))
				
		if (new.fitness >= fitness){
			## A postprocessing stage
			best.chromosome <- rule.gen
			new.fitness <- fitness
			fit[iter, 1] <- c(new.fitness)
		}
	
		##mutation
		rule.gen.mut <- rule.gen
		for (i in 1 : nrow(rule.gen.mut)){
			loc.mutation <- runif(ncol(rule.gen.mut), min=0, max = 1)		
			for (j in 1 : ncol(rule.gen.mut)){
				if (loc.mutation[j] <= persen_mutant){
					cond <- round(runif(1, min = 0, max = 1))
					r <- runif(1, min = 0, max = 1)
					b = 1
					if (cond == 0){
						y <- 1 - rule.gen.mut[i, j]
						delta <- y * (1 -  r ^ ((1 - iter/max.iter) ^ 1))
						rule.gen.mut[i, j] <- rule.gen.mut[i, j] + delta
						
						## check whether the new chromosome is better with the old one
						ind.fit.mut <- calc.ind.fitness(data.train.ori, matrix(rule.gen.mut[i, ], nrow = 1))
						if (ind.fit.mut < rule.ind.fitness[i, ncol(rule.ind.fitness)]){
							rule.ind.fitness[i, j] <- rule.gen.mut[i, j]
							rule.ind.fitness[i, ncol(rule.ind.fitness)] <- ind.fit.mut
						}	
					}
					else {
						y <- rule.gen.mut[i, j]
					 	delta <- y * (1 -  r ^ ((1 - iter/max.iter) ^ 1))
						rule.gen.mut[i, j] <- rule.gen.mut[i, j] - delta
						
						## check whether the new chromosome is better with the old one
						ind.fit.mut <- calc.ind.fitness(data.train.ori, matrix(rule.gen.mut[i, ], nrow = 1))
						if (ind.fit.mut < rule.ind.fitness[i, ncol(rule.ind.fitness)]){
							rule.ind.fitness[i, j] <- rule.gen.mut[i, j]
							rule.ind.fitness[i, ncol(rule.ind.fitness)] <- ind.fit.mut
						}
					}
				}
			}
		}
		
		### crossover 
		## determine the number of parents
		num.parents <- as.integer(runif(1, min = 1, max = nrow(rule.gen)/2))
		parent.seq <- sample(1 : nrow(rule.gen), size = (num.parents * 2))
		
		for (i in seq(1, length(parent.seq), 2)){
			gen.rand <- runif(1, min = 0, max = 1)
			if (gen.rand <= persen_cross){
				cross.point.rule <- as.integer(runif(1, min = 1, max = ncol(rule.gen)))
				
				child.rule.gen.a1 <- matrix(rule.gen[parent.seq[i], 1 : cross.point.rule], nrow = 1)
				child.rule.gen.a2 <- matrix(rule.gen[parent.seq[i + 1], (cross.point.rule + 1) : ncol(rule.gen)], nrow = 1)
				child.rule.gen.a <- cbind(child.rule.gen.a1, child.rule.gen.a2)
				
				child.rule.gen.b1 <- matrix(rule.gen[parent.seq[i], (cross.point.rule + 1) : ncol(rule.gen)], nrow = 1)
				child.rule.gen.b2 <- matrix(rule.gen[parent.seq[i + 1], 1 : cross.point.rule], nrow = 1)
				child.rule.gen.b <- cbind(child.rule.gen.b1, child.rule.gen.b2)			
				
				## check fitness of child
				ind.fit.cr.a <- calc.ind.fitness(data.train.ori, child.rule.gen.a)
				ind.fit.cr.b <- calc.ind.fitness(data.train.ori, child.rule.gen.b)
				
				if (ind.fit.cr.a < rule.ind.fitness[parent.seq[i], ncol(rule.ind.fitness)]){
					rule.ind.fitness[parent.seq[i], (1 : ncol(rule.ind.fitness) - 1)] <- child.rule.gen.a
					rule.ind.fitness[parent.seq[i], ncol(rule.ind.fitness)] <- ind.fit.cr.a
				}
				
				if (ind.fit.cr.b < rule.ind.fitness[parent.seq[i + 1], ncol(rule.ind.fitness)]){
					rule.ind.fitness[parent.seq[i + 1], (1 : ncol(rule.ind.fitness) - 1)] <- child.rule.gen.b
					rule.ind.fitness[parent.seq[i + 1], ncol(rule.ind.fitness)] <- ind.fit.cr.b
				}
			}
		}
							
	iter <- iter + 1
	}

	## A genetic tuning stage
#plot(fit)
rule <- list(rule = best.chromosome)
mod <- list(rule = best.chromosome, range.data.ori = range.data.ori)	
return(mod)
}


# This function is to generate chromosome representing fuzzy IF-THEN rules. 
#
# @title The chromosome generating function
# @param data.i a matrix(1 x m) of the normalized data.
# @param epsilon a parameter of covering factor.  
# @return rules
# @export
genetic.gen <- function(data.i, epsilon){
	
	## now, this method just implement triangular membership function
	rule.gen <- matrix()
	for (i in 1 : ncol(data.i)){					
		dt.i <- data.i[1, i]
		delta.data <- min(dt.i, 1 - dt.i)
		rand.delta <- runif(1, min = 0, max = delta.data)
		
		temp.rule.i <- matrix(c((dt.i - rand.delta), dt.i, (dt.i + rand.delta)), nrow = 1)
		
		rule.gen <- cbind(rule.gen, temp.rule.i)
	}

	rule.gen <- rule.gen[1, 2 : ncol(rule.gen)]
	
	return(rule.gen)
}

# This function is to calculate compatibility degrees. 
#
# @title The compability degree
# @param data.train.i a matrix(1 x m) of the normalized data.
# @param rule.gen chromosome which represents fuzzy IF-THEN rules.  
# @return The compability degree
# @export
ch.cover <- function(data.train.i, rule.gen){
## calculate degree of MF using chromosome ==> calc.compDegree
comp.degree <- matrix(nrow = nrow(rule.gen))

for (i in 1 : nrow(rule.gen)){
	# do t-norm over all variables
	degree.var <- calc.compDegree(data.train.i, matrix(rule.gen[i, ], nrow = 1))
	comp.degree[i, 1] <- min(degree.var)
}		
return(comp.degree)	
}

# This function is to calculate membership function (MF) degrees. 
#
# @title The MF degree
# @param data.i a matrix(1 x m) of the normalized data.
# @param rule.gen.i a chromosome which represents fuzzy IF-THEN rules.  
# @return The degree
# @export
calc.compDegree <- function(data.i, rule.gen.i){
	## calculate degree of MF

	ncol.data.i <- ncol(data.i)
	degree.R <- matrix(nrow = 1, ncol = ncol.data.i)
	for (j in 1 : ncol.data.i){
		xx <- data.i[1, j]
		aa <- rule.gen.i[1, (j * 3 - 2)]
		bb <- rule.gen.i[1, (j * 3 - 1)]
		cc <- rule.gen.i[1, (j * 3)]
		
		if (xx <= aa){
			degree.R[1, j] <- 0
		}
		else if (xx <= bb) {
			degree.R[1, j] <- (xx - aa) / (bb - aa)
		}
		else if (xx <= cc) {
			degree.R[1, j] <- (cc - xx) / (cc - bb)
		}
		else {
			degree.R[1, j] <- 0
		}
	}
return(degree.R)
}

# This function is to calculate fitness. 
#
# @title The fitness function
# @param data.test.i a matrix(1 x m) of the normalized data.
# @param target.val.i a target value
# @param rule.gen the chromosomes which represent fuzzy IF-THEN rules.  
# @return The residual
# @export
calc.fitness <- function(data.test.i, target.val.i, rule.gen){
	output.val.pred <- calc.pred.val(data.test.i, rule.gen)
	benchmarks <- cbind(output.val.pred, target.val.i)
	residuals <- (benchmarks[, 1] - benchmarks[, 2])^2
	fit <- residuals
	return(fit)
}

# This function is to calculate fitness of each chromosome. 
#
# @title The individual fitness
# @param data.train a matrix(n x m) of the normalized data.
# @param rule.gen the chromosomes which represent fuzzy IF-THEN rules.  
# @return The individual fitness
# @export
calc.ind.fitness <- function(data.train, rule.gen){
	data.test <- matrix(data.train[, 1 : (ncol(data.train) - 1)], nrow = nrow(data.train))
	real.val <- matrix(data.train[, ncol(data.train)], nrow = nrow(data.train))
	ind.fit <- matrix(nrow = nrow(rule.gen), ncol = 1)
	cov.deg <- matrix(nrow = nrow(rule.gen), ncol = 1)
	R <- matrix(nrow = nrow(rule.gen), ncol = nrow(data.train))
	for (i in 1 : nrow(rule.gen)){
		mod <- list(rule = matrix(rule.gen[i, ], nrow =1))
		res <- MSGFS.test(mod, data.test)
		y.pred <- res
		y.real <- real.val

		bench <- cbind(y.pred, y.real)

		#### Measure error for classification
		residuals <- (y.real - y.pred)
		#MSE <- mean(residuals^2)
		RMSE <- sqrt(mean(residuals^2))
		
		ind.fit[i, 1] <- RMSE	
	}
	
	### check high-frequency value
	for(i in 1 : nrow(rule.gen)){
		for (j in 1 : nrow(data.train)){
			R[i, j] <- ch.cover(matrix(data.train[j, ], nrow = 1), matrix(rule.gen[i, ], nrow = 1))
		}
	}
	
	cov.deg <- matrix(rowSums(R)/nrow(data.train), ncol = 1)
	
	res.ind.fit <- ind.fit / (1 + cov.deg)
	return(res.ind.fit)
}