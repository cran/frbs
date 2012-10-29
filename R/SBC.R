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
#' This is the internal function that implements a combination of the subtractive 
#' clustering method and fuzzy c-means. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}
#' 
#' This method was proposed by Stephen Chiu and J.C. Bezdek. For generating the rules 
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
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param range.data.ori a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param r.a the radius defining a neighborhood. 
#' @param eps.high an upper threshold value. 
#' @param eps.low a lower threshold value. 
#' @seealso \code{\link{SBC.test}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components.
#' @references
#' Bezdek, J.C., "Pattern recognition with fuzzy objective function algorithms", Plenum Press, New York (1981).
#'
#' Nikhil R. Pal, James C. Bezdek, and Richard J. Hathaway, "Sequential competitive learning and the fuzzy c-means clustering algorithms," 
#' Neural Networks, Vol. 9, No. 5, pp. 787-796 (1996).
#'
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
	mod <- list(cls = cls, range.data.ori = range.data.ori, r.a = r.a)
	return(mod)
} 

# This function is to calculate the potential of data point that will be used on subtractive clustering (SBC) method. 
#
# @title The potential of data
# @param dt.norm a matrix(n x m) of the normalized data.
# @param alpha a parameter of distance  
# @seealso \code{\link{SBC}} and \code{\link{SBC.test}}
# @return the potential on each data
# @export
potential <- function(dt.norm, alpha){
	counter <- 0
	n <- nrow(dt.norm)
	potential <- matrix(nrow = n)
	
	for (i in 1:n){
		temp.x <- dt.norm[i,]
		for (j in 1:n){
			euc <- dist(rbind(temp.x, dt.norm[j,]))
			cum <- exp(-alpha * (euc)^2)
			counter <- cum + counter
		}
		potential[i] <- counter
		counter <- 0
	}	
	return(potential)
}

# This function is to revise of the potential. 
#
# @title The revision of the potential
# @param P.star a value of potential of cluster center
# @param x.star a vector of cluster center
# @param dt.norm a matrix(n x m) of the normalized data.
# @param Beta a parameter of distance  
# @seealso \code{\link{SBC}} and \code{\link{SBC.test}}
# @return the revision of potential on each data
# @export
revise.pot <- function(P.star, x.star, dt.norm, Beta){
	counter <- 0
	n <- nrow(dt.norm)
	potential <- matrix(nrow = n)
	
	temp.x <- x.star
	for (i in 1:n){
		euc <- dist(rbind(temp.x, dt.norm[i,]))
		potential[i] <- exp(-Beta * (euc)^2)		
	}
	
	rev.pot <- P.star * potential 
	return(rev.pot)
}

# This function is to get stopping criteria on SBC
#
# @title The stopping criteria on Subtractive Clustering (SBC)
# @param r.a a positive constant which is effectively the radius defining a neighborhood. (default = 0.5)
# @param eps.high a parameter which is upper threshold value. (default = 0.5)
# @param eps.low a parameter which is lower threshold value. (default = 0.15)
# @param PP.star a value of potential of cluster center on k-th iteration
# @param P.star a value of potential of cluster center on 1-st iteration
# @param cluster.ctr a matrix of cluster center
# seealso \code{\link{SBC.test}} and \code{\link{SBC}}
# @return stopping criteria
# @export
stop.criteria <- function(r.a=0.5, eps.high = 0.5, eps.low = 0.15, PP.star, P.star,  cluster.ctr){
	temp.c.ctr <- na.omit(cluster.ctr)
	
	if (PP.star > (eps.high * P.star)){
		st.cr <- c(1)
	}
	else if (PP.star < (eps.low * P.star)){
		st.cr <- c(2)
	}
	else {
		x.k <- temp.c.ctr[nrow(temp.c.ctr), ]
		for (i in 1: (nrow(temp.c.ctr) - 1)){
			dmin <- dist(rbind(x.k, temp.c.ctr[i, ]))
			crt <- ((dmin/r.a) + (PP.star/P.star))
			if (crt >= 1){
				st.cr <- c(1)
			}
			else {
				st.cr <- c(3)
			}
		}
	}
	return(st.cr)
}