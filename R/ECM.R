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
#' This function is a part of the DENFIS method to generate cluster centers. 
#'
#' @title Envolving Clustering Method
#'
#' @param data.train a matrix(m x n) of data for training, where m is the number of instances and 
#' n is the number of variables where the last column is the output variable.
#' @param Dthr the threshold value for the envolving clustering method (ECM), between 0 and 1.
#' @seealso \code{\link{DENFIS}} and \code{\link{DENFIS.eng}}
#' @return a matrix of cluster centers
#' @export
ECM <- function(data.train, Dthr){

Cc1 <- matrix(data.train[1, ], nrow = 1)

Ru1 <- 0
nrow.dt <- nrow(data.train)
ncol.dt <- ncol(data.train)
Cc.j <- Cc1
Ru.j <- matrix(Ru1)

D.ij <- matrix()
temp <- matrix()
temp <- c(1)
n.cc <- 1
D.ij <- matrix()

num.cls <- 1
num.mem <- 1
for (i in 2 : nrow.dt){
	for (j in 1 : nrow(Cc.j)){
		norm.euc <- dist(rbind(data.train[i,], Cc.j[j, ]))
		D.ij[j] <- (norm.euc)/sqrt(ncol.dt)
	}
	
	D.ij <- matrix(D.ij)
	indx <- which.min(D.ij)	
	
	if (any(D.ij[indx] <= Ru.j)){
		#do nothing
	}
	else{
		S.ij <- Ru.j + D.ij
		indx.s <- which.min(S.ij)
		conts <- (2 * Dthr)
		if (S.ij[indx.s] > conts){
			Cc.j <- rbind(Cc.j, data.train[i, ])
			Ru.j <- rbind(Ru.j, 0)
		}
		else {		
			Ru.j.temp <- 0.5 * S.ij[indx.s]	
			
			if (Ru.j.temp > Dthr){
				Cc.j <- rbind(Cc.j, data.train[i, ])
				Ru.j <- rbind(Ru.j, 0)
			}
			else {
				temp <- data.train[i, ] - Cc.j[indx.s, ]
				d.temp <- sqrt(sum(temp^2))
				ratio <- abs((d.temp - Ru.j.temp))/d.temp
				new.vec <- ratio * temp
				new.cls <- Cc.j[indx.s, ] + new.vec
				Cc.j[indx.s, ] <- new.cls
				Ru.j[indx.s] <- Ru.j.temp
			}
		}
	}
}

res <- Cc.j

return(res)
}