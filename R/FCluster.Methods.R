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
#' This is the internal function that implements the dynamic evolving neural-fuzzy inference system (DENFIS). 
#' It is used to handle regression tasks. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by Nikola K. Kasabov and Qun Song. There are several steps in this method that 
#' are to determine the cluster centers using the evolving clustering method (ECM), to partition the input space 
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
#' @param data.train a matrix (m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables (input and output variables).
#' @param range.data.ori a matrix (2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param Dthr the threshold value for the evolving clustering method (ECM), between 0 and 1. 
#' @param max.iter the maximal number of iterations.
#' @param step.size the step size of the least squares method, between 0 and 1.
#' @param d a parameter for the width of the triangular membership function.
#' @seealso \code{\link{DENFIS.eng}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return list of the model data. Please have a look at \code{\link{frbs.learn}} for looking complete components. 
#' @references
#' N. K. Kasabov and Q. Song, "DENFIS: dynamic evolving neural-fuzzy inference system and its Application for time-series prediction", 
#' IEEE Transactions on Fuzzy Systems, vol. 10, no. 2, pp. 144 - 154 (2002).
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

dt.input <- data.train[, 1 : (ncol(data.train) - 1), drop = FALSE]

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
i <- seq(from = 1, to = ncol(data.train))
colnames(cluster.c) <- paste("var", i, sep = ".")
colnames(range.data.ori) <- paste("var", i, sep = ".")
mod <- list(cls = cluster.c, func.tsk = func.tsk, range.data.ori = range.data.ori, Dthr = Dthr, d = d, type.model = "CLUSTERING")

return(mod)
}

#' This is the internal function that implements a combination of the subtractive 
#' clustering method and fuzzy c-means. It is used to solve regression tasks. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}
#' 
#' This method was proposed by Stephen Chiu. For generating the rules 
#' in the learning phase, the subtractive clustering method is used to obtain the cluster 
#' centers. Subtractive clustering (SBC) is an extension of Yager and Filev's 
#' mountain method.
#' SBC considers each data point as a potential cluster center by determining 
#' the potential of a data point as a function of its distances to all the other 
#' data points. A data point has a high potential value if that data 
#' point has many nearby neighbors.
#' The highest potential is chosen as the cluster center and then the potential of 
#' each data point will be updated. The process of determining new clusters and updating 
#' potentials repeats until the remaining potential of all data points falls below 
#' some fraction  
#' of the potential of the first cluster center. After getting all the cluster 
#' centers from subtractive clustering,
#' the cluster centers are optimized by fuzzy c-means. 
#' 
#' @title The subtractive clustering and fuzzy c-means (SBC) model building 
#' @param data.train a matrix (m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param range.data.ori a matrix (2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param r.a the radius defining a neighborhood. 
#' @param eps.high an upper threshold value. 
#' @param eps.low a lower threshold value. 
#' @seealso \code{\link{SBC.test}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components.
#' @references
#' R. Yager and D. Filev, "Generation of fuzzy rules by mountain clustering," 
#' J. of Intelligent and Fuzzy Systems, vol. 2, no. 3, pp. 209 - 219 (1994).
#' 
#' Stephen Chiu, "Method and software for extracting fuzzy classification rules by subtractive clustering", 
#' Fuzzy Information Processing Society, NAFIPS, pp. 461 - 465 (1996).
# @export
SBC <- function(data.train, range.data.ori, r.a = 0.5, eps.high = 0.5, eps.low = 0.15){

	req.suc <- require("e1071", quietly=TRUE)
	if(!req.suc) stop("In order to use this function, you need to install the package e1071.")
	
	num.outvar = 1
	
	#c.center <- matrix(ncol, 1) 
	num.inpvar <- (ncol(data.train) - num.outvar)
	cluster.ctr <- matrix(NA, nrow = nrow(data.train), ncol = ncol(data.train))

	## normalization of data into [-1, 1]
	data.norm <- norm.data(data.train, range.data.ori, min.scale = -1, max.scale = 1)
	
	## calculate potential of data points
	alpha <- 4 / ((r.a)^2)
	pot <- potential(data.norm, alpha) 
	
	## the first cluster center
	P.star <- max(pot)
	indx.P.star <- which.max(pot)
	x.star <- data.norm[indx.P.star, ]
	
	cluster.ctr[1, ] <- x.star
	
	Beta <- 4 / (1.5 * r.a)^2
	new.pot <- pot
	
	for(iii in 2 : nrow(data.norm)){
		
		# calculate revision of potential on each point		
		rev.pot <- revise.pot(P.star, x.star, data.norm, Beta)
		
		new.pot <- (new.pot - rev.pot)
		
		PP.star <- max(new.pot)
		indx.PP.star <- which.max(new.pot)
		xx.star <- data.norm[indx.PP.star, ]
		
		cluster.ctr[iii, ] <- xx.star
		
		st.cr <- stop.criteria(r.a, eps.high = 0.5, eps.low = 0.15, PP.star, P.star, cluster.ctr)
		
		if (st.cr == 1){	
			P.star <- PP.star
			x.star <- cluster.ctr[iii, ]
		}
		else if (st.cr == 2){
			cluster.ctr[iii, ] <- NA	
			break		
		}
		else if (st.cr == 3) {
			cluster.ctr[iii, ] <- NA
			new.pot[indx.PP.star] <- 0
			PP.star.2 <- max(new.pot)
			indx.PP.star <- which.max(new.pot)
			xx.star <- data.norm[indx.PP.star, ]		
			cluster.ctr[iii, ] <- xx.star
			st.cr <- stop.criteria(r.a, eps.high = 0.5, eps.low = 0.15, PP.star.2, P.star, cluster.ctr)
			if (st.cr == 1){
				P.star <- PP.star.2
				x.star <- cluster.ctr[iii, ]
			}
			else if (st.cr == 2){
				cluster.ctr[iii, ] <- NA
				break
			}
			else if (st.cr == 3) {
				cluster.ctr[iii, ] <- NA
				break
			}
		}
		
	}
	
	
	res.c.ctr <- na.omit(cluster.ctr)	
	cl.fcm <- cmeans(data.norm, res.c.ctr, 300, method="cmeans")
	cls <- denorm.data(cl.fcm$centers, range.data.ori, min.scale = -1, max.scale = 1)	
	i <- seq(from = 1, to = ncol(data.train))
	colnames(cls) <- paste("var", i, sep = ".")
	colnames(range.data.ori) <- paste("var", i, sep = ".")
	mod <- list(cls = cls, range.data.ori = range.data.ori, r.a = r.a, type.model = "CLUSTERING")
	return(mod)
} 

