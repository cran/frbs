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
# This function is to generation a new population.
#
# @title the population generating function
# @param popu.all the population on GA representing the rule database and membership function 
# @param popu.elit the best individu on previous generation 
# @param persen_cross  a real number representing a crossover percentage
# @param persen_mutant  a real number representing a mutation percentage
# @param n.labels number of label used as linguistic fuzzy term
# @param num.var a number of variables on data training
# @return a list of the population (MF parameters and rule database)
# @export
popubaru <- function(popu.all, popu.elit, persen_cross, persen_mutant, n.labels, num.var){
	
	MF.param.popu <- popu.all$MF.param.popu
	rule.popu <- popu.all$rule.popu
	popu.size <- nrow(MF.param.popu)
	
	length.MF <- ncol(MF.param.popu)
	length.rule <- ncol(rule.popu)
	
	MF.elit <- popu.elit$MF.elit.ind
	rule.elit <- popu.elit$rule.elit.ind
	
	num.popu <- nrow(MF.param.popu)
	
	num.parent <- floor(num.popu / 3)
	
	rand.indx <- matrix(as.integer(runif(num.parent, min = 1, max = num.popu)), nrow = 1)
	
	child.MF.a1 <- matrix(nrow = 1)
	child.MF.a2 <- matrix(nrow = 1)
	child.MF.b1 <- matrix(nrow = 1)
	child.MF.b2 <- matrix(nrow = 1)
	child.MF.a <- matrix(nrow = 1)
	child.MF.b <- matrix(nrow = 1)
	MF.new.popu <- matrix(ncol = length.MF)
	rule.new.popu <- matrix(ncol = length.rule)
	
	##selection parent and crossover
	if (length(rand.indx) != 0){	
		for (i in 1 : length(rand.indx)){
			gen.rand <- runif(1, min = 0, max = 1)
			if (gen.rand <= persen_cross){
				cross.point.MF <- as.integer(runif(1, min = 1, max = length.MF))
				cross.point.rule <- as.integer(runif(1, min = 1, max = length.rule))
				
				child.MF.a1 <- matrix(MF.elit[1, 1 : cross.point.MF], nrow = 1)
				child.MF.a2 <- matrix(MF.param.popu[rand.indx[1, i], (cross.point.MF + 1) : length.MF], nrow = 1)
				child.MF.a <- cbind(child.MF.a1, child.MF.a2)
			
				child.rule.a1 <- matrix(rule.elit[1, 1 : cross.point.rule], nrow = 1)
				child.rule.a2 <- matrix(rule.popu[rand.indx[1, i], (cross.point.rule + 1) : length.rule], nrow = 1)
				child.rule.a <- cbind(child.rule.a1, child.rule.a2)

				MF.new.popu <- rbind(MF.new.popu, child.MF.a)
				rule.new.popu <- rbind(rule.new.popu, child.rule.a)
	
				child.MF.b1 <- matrix(MF.elit[1, (cross.point.MF + 1) : length.MF], nrow = 1)
				child.MF.b2 <- matrix(MF.param.popu[rand.indx[1, i], 1 : cross.point.MF], nrow = 1)
				child.MF.b <- cbind(child.MF.b2, child.MF.b1)
				
				child.rule.b1 <- matrix(rule.elit[1, (cross.point.rule + 1) : length.rule], nrow = 1)
				child.rule.b2 <- matrix(rule.popu[rand.indx[1, i], 1 : cross.point.rule], nrow = 1)
				child.rule.b <- cbind(child.rule.b2, child.rule.b1)
		
				MF.new.popu <- rbind(MF.new.popu, child.MF.b)
				rule.new.popu <- rbind(rule.new.popu, child.rule.b)
			}
		}
	}
	
	MF.new.popu <- na.omit(MF.new.popu)
	rule.new.popu <- na.omit(rule.new.popu)
	
	MF.new.popu <- rbind(MF.new.popu, MF.param.popu)
	rule.new.popu <- rbind(rule.new.popu, rule.popu)
	
	length.MF.new.popu <- nrow(MF.new.popu)
	length.rule.new.popu <- nrow(rule.new.popu)

	## keep to be the same size
	if (length.MF.new.popu > popu.size) {
		num.del <- (length.MF.new.popu - popu.size)
		seqq <- seq(1, num.del, 1)
		MF.new.popu <- MF.new.popu[-seqq, ]
		rule.new.popu <- rule.new.popu[-seqq, ]
	}
	#####Mutation of MF
	for (i in 1 : popu.size){
		loc.mutation <- runif(length.MF, min=0, max = 1)
		leap <- 0
		for (j in 1 : length(loc.mutation)){
			if (loc.mutation[j] <= persen_mutant){
				if (j %% 2 == 1) {
					eps <- runif(1, min = 0.005, max = 0.3)
					MF.new.popu[i, j] <- eps	
				}
				else {
					loc.low <- leap * n.labels * 2 + 2
					if (j == loc.low ) {
						min.eps <- 0
						max.eps <- 1
						eps <- runif(1, min = min.eps, max = max.eps)
						MF.new.popu[i, j] <- eps
						leap <- leap + 1
					}
					else {
						min.eps <- MF.new.popu[i, j - 2]
						max.eps <- 1
						eps <- runif(1, min = min.eps, max = max.eps)
						MF.new.popu[i, j] <- eps
					}
				}
			}
		}
		
		### mutation of rule
		loc.mutation.rule <- as.integer(runif(1, min=1, max = length.rule))
		loc.var <- loc.mutation.rule %% num.var
		if (loc.var == 0){
			end.mut <- num.var * n.labels
			start.mut <- end.mut - n.labels	+ 1		
		}
		else {
			end.mut <- loc.var * n.labels
			start.mut <- end.mut - n.labels + 1
		}
			label.rand <- as.integer(runif(1, min = start.mut, max = end.mut))
			rule.new.popu[i, loc.mutation.rule] <- label.rand
	}
	
	#### creep operator
	for (i in 1 : popu.size){
		creep.point <- matrix(round(runif(length.MF, min = 0, max = 1)), nrow = 1)
		eps <- matrix(runif(length.MF, min = -0.05, max = 0.05), nrow = 1)		
		creep.eps <- creep.point * eps
		
		for (j in 1 : length.MF){
			temp = MF.new.popu[i, j] + creep.eps[1, j] 
			if (temp > 0 && temp < 1)
				MF.new.popu[i, j] <- MF.new.popu[i, j] + creep.eps[1, j]
		}
	}
	
	MF.new.popu[1, ] <- MF.elit[1, ]
	rule.new.popu[1, ] <- rule.elit[1, ]
	
	new.popu <- list(MF.new.popu = MF.new.popu, rule.new.popu = rule.new.popu)
	
return(new.popu)
}