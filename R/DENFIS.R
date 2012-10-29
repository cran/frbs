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
#' This is the internal function that implements the DENFIS method. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by Nikola K. Kasabov and Qun Song. There are several steps in this method that 
#' are to determine the cluster centers using the envolving clustering method (ECM) to partition the input space 
#' and to find optimal parameters on the consequent part (Takagi Sugeno Kang model) for the IF-THEN rule using a least 
#' squares estimator. 
#'
#' ECM is a distance-based clustering method which is determined by a threshold value, \code{Dthr}. This parameter 
#' influences how many clusters are created. In the beginning of the clustering process, the first instance from the 
#' training data is chosen to be a cluster center, and the determining radius is set to zero. Afterwards, using the 
#' next instance, cluster centers and radius are changed based on certain mechanisms of ECM (please see \code{\link{ECM}}). 
#' All of the cluster centers are then obtained after evaluating all the training data. 
#' The next step is to update the parameters on the consequent part with the assumption that the antecedent part which we got from ECM is fixed. 
#' Actually, ECM can perform well as an online clustering method, but in this package it is used in an offline mode. 
#' 
#' @title DENFIS model building
#'
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables (input and output variables).
#' @param range.data.ori a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param Dthr the threshold value for the envolving clustering method (ECM), between 0 and 1. 
#' @param max.iter the maximal number of iterations.
#' @param step.size the step size of the least squares method, between 0 and 1.
#' @param d a parameter for the width of the triangular membership function.
#' @seealso \code{\link{DENFIS.eng}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return list of the model data. Please have a look at \code{\link{frbs.learn}} for looking complete components. 
#' @references
#' Nikola K. Kasabov and Qun Song, "DENFIS: dynamic evolving neural-fuzzy inference system and its Application for time-series prediction", 
#' IEEE Transactions on Fuzzy Systems, Vol. 10, No. 2, pp. 144 - 154 (2002).
# @export
DENFIS <- function(data.train, range.data.ori, Dthr = 0.1, max.iter = 100, step.size = 0.01, d = 2){

min.scale <- 0
max.scale <- 1

data.train <- norm.data(data.train, range.data.ori, min.scale, max.scale)
alpha <- step.size
cluster.c <- ECM(data.train, Dthr)
MSE <- matrix()
num.cls <- nrow(cluster.c)
num.dt <- nrow(data.train) 
num.inpvar <- (ncol(data.train) - 1)
temp <- matrix(nrow = num.cls, ncol=num.inpvar)
miu.rule <- matrix(nrow = num.dt, ncol = num.cls)

def <- matrix(nrow=num.dt, ncol = 1)

for (i in 1 : num.dt){
	for (j in 1 : num.cls){
		for (k in 1 : num.inpvar){
			b <- cluster.c[j, k]
			a <- b - d * Dthr
			cc <- b + d * Dthr
			x <- data.train[i, k]
			left <- (x - a)/(b - a)
			right <- (cc - x)/(cc - b)
			temp[j, k] <- max(min(left, right), 0)
		}
		miu.rule[i, j] <- prod(temp[j, ])	
	}	
}

num.ele <- num.cls * (num.inpvar + 1)
rand <- runif(num.ele, min = 0, max = 1)
func.tsk <- matrix(rand, nrow = num.cls, ncol= (num.inpvar + 1), byrow=T)

dt.input <- matrix(data.train[, 1 : (ncol(data.train) - 1)], nrow = num.dt)

for (i in 1 : max.iter) {
### Calculate defuzzification
	for (j in 1 : num.dt) {
		data.k <- (dt.input[j, ])
		data.m <- as.matrix(data.k)
		
		func.tsk.var <- as.matrix(func.tsk[, 1: (ncol(func.tsk) - 1)])
		func.tsk.cont <- matrix(func.tsk[, ncol(func.tsk)], ncol = 1)
		
		ff <- func.tsk.var %*% data.m + func.tsk.cont
		
		miu.rule.t <- as.matrix(miu.rule[j, ])
		cum <- sum(ff * miu.rule.t)
		div <- sum(miu.rule.t)
		
		if (div == 0)
			def[j] <- 0
		else {
			def[j] <- cum / div
		}
	}

gal1 <- def - data.train[, ncol(data.train)]

for (ii in 1 : num.dt){
	gal <- (def[ii] - data.train[ii, ncol(data.train)])
	
	for (m in 1 : nrow(func.tsk.var)){
		for (n in 1 : ncol(func.tsk.var)){
			sum.miu <- sum(miu.rule[ii, ], na.rm = TRUE)
			if (sum.miu != 0)
				func.tsk.var[m, n] <- func.tsk.var[m, n] - alpha * gal * (miu.rule[ii, m]/sum.miu) *  data.train[ii, n]	
			else
				func.tsk.var[m, n] <- func.tsk.var[m, n] - alpha * gal * (miu.rule[ii, m]) *  data.train[ii, n]	
		}	
	}
	
	for (mm in 1 : nrow(func.tsk.cont)){
		sum.miu <- sum(miu.rule[ii, ], na.rm = TRUE)
		if (sum.miu != 0)
			func.tsk.cont[mm, 1] <- func.tsk.cont[mm, 1] - alpha * gal * (miu.rule[ii, m]/sum.miu)
		else
			func.tsk.cont[mm, 1] <- func.tsk.cont[mm, 1] - alpha * gal * miu.rule[ii, m]
 	}
}

func.tsk.new <- cbind(func.tsk.var, func.tsk.cont)
func.tsk <- func.tsk.new

	
}	

cluster.c <- denorm.data(cluster.c, range.data.ori, min.scale, max.scale)

mod <- list(cls = cluster.c, func.tsk = func.tsk, range.data.ori = range.data.ori, Dthr = Dthr, d = d)

return(mod)
}